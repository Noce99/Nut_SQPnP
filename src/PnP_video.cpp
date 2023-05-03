#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <opencv2/aruco.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/objdetect/aruco_detector.hpp>
#include <opencv2/calib3d.hpp>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "pnp_solver.cpp"

#define DATA_PATH "/home/enrico/Progetti/Nut_SQPnP/resources/3D_points.data"
#define VIDEO_PATH "/home/enrico/Progetti/Nut_SQPnP/resources/test_video.mp4"
#define OUTPUT_PATH "/home/enrico/Progetti/Nut_SQPnP/resources/camera_positions.txt"
#define OUTPUT_VIDEO_PATH "/home/enrico/Progetti/Nut_SQPnP/resources/video_with_BBs.avi"
int main () {
	/*
		START PART 1
		There I read the file where I have putted on each line the id of a marker
		followed by 4 3D coordinates that are the 4 corners coordinates mesured in
		my room with a meter.
	*/
	std::string line;
  std::ifstream myfile(DATA_PATH);
  std::vector<std::vector<Eigen::Matrix<double, 3, 1>>> marker_points;
  std::vector<int> marker_id;
  if (myfile.is_open()){
  	while (getline(myfile,line)){
    	int start = 0;
    	std::vector<double> line_of_data;
    	for (int i=0; i<line.length(); i++){
    		if (line[i] == ' '){
    			line_of_data.push_back(std::stod(line.substr(start, i-start)));
    			start = i+1;
      	}
    	}
    	line_of_data.push_back(std::stod(line.substr(start, line.length())));
    	int id = int(line_of_data[0]);
    	std::vector<Eigen::Matrix<double, 3, 1>> tmp;
    	for (int i=0; i<4; i++){
    		tmp.push_back((Eigen::Matrix<double, 3, 1>() <<
											 line_of_data[1 + 3*i],
    									 line_of_data[1 + 3*i + 1],
    									 line_of_data[1 + 3*i + 2]).finished());
    	}
    	marker_points.push_back(tmp);
    	marker_id.push_back(id);
    }
    myfile.close();
  }else std::cout << "Unable to open file" << std::endl;
	/*
		FINISHED PART 1
	*/

	/*
		START PART 2
		There I undistort the image given the Camera Matrix and the Distortion Coefficients
		obtained with the camera calibration process
	*/
  cv::aruco::DetectorParameters detectorParams = cv::aruco::DetectorParameters();
	cv::aruco::Dictionary dictionary = cv::aruco::getPredefinedDictionary(cv::aruco::DICT_APRILTAG_16h5);
	cv::aruco::ArucoDetector detector(dictionary, detectorParams);

	Eigen::Matrix<double, 3, 3> K{(Eigen::Matrix<double, 3, 3>()  << 9.71143e+02, 0.00000e+00, 6.5138910e+02,
										   							 0.00000e+00, 9.63718e+02, 3.63198e+02,
										   							 0.00000e+00, 0.00000e+00, 1).finished()
							 	 };

	cv::Mat Camera_Matrix = cv::Mat(3, 3, CV_64F);
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			Camera_Matrix.at<double>(i, j) = K(i, j);
		}
	}
	std::vector<double> Distortion_Coefficients;
	Distortion_Coefficients.push_back(-1.417623607107070319e-02);
	Distortion_Coefficients.push_back(4.599614520141423468e-01);
	Distortion_Coefficients.push_back(1.775918531563122323e-03);
	Distortion_Coefficients.push_back(2.065892344481764023e-03);
	Distortion_Coefficients.push_back(-1.038305182066649568e+00);
	/*
		FINISH PART 2
	*/

	// Create a VideoCapture object and open the input file
  // If the input is the web camera, pass 0 instead of the video file name
  cv::VideoCapture cap(VIDEO_PATH);
	std::vector<Eigen::Matrix<double, 4, 1>> camera_positions;
  // Check if camera opened successfully
  if(!cap.isOpened()){
    std::cout << "Error opening video stream or file" << std::endl;
    return -1;
  }

	Eigen::Matrix<double, 4, 1> last_point = (Eigen::Matrix<double, 4, 1>()  << 0, 0, 0, 0).finished();

	cv::VideoWriter video(OUTPUT_VIDEO_PATH, cv::VideoWriter::fourcc('M','J','P','G'), 10, cv::Size(720,1280));

  while(1){
    cv::Mat frame;
    cap >> frame;

    if (frame.empty())
      break;

		auto start = std::chrono::high_resolution_clock::now();
		cv::Mat frameUndistorted;
		cv::undistort(frame, frameUndistorted, Camera_Matrix, Distortion_Coefficients);

		std::vector<int> ids;
		std::vector<std::vector<cv::Point2f>> corners, rejected;
		detector.detectMarkers(frame, corners, ids, rejected);

		std::vector<Eigen::Matrix<double, 3, 1>> points_3D;
		std::vector<Eigen::Matrix<double, 2, 1>> points_2D;


		for (int i=0; i<ids.size(); i++){
			int id = ids[i];
			for (int j=0; j<marker_id.size(); j++){
				if (marker_id[j] == id){
					for (int ii=0; ii<4; ii++){
						points_3D.push_back(marker_points[j][ii]);
						points_2D.push_back((Eigen::Matrix<double, 2, 1>() << corners[i][ii].x, corners[i][ii].y).finished());
					}
				}
			}
		}

		if (points_3D.size() > 0){
			cv::aruco::drawDetectedMarkers(frame, corners, ids);
			std::cout << "I have [" << points_3D.size() << "] 3D points, [" << points_2D.size() << "] 2D points and a (" << K.rows() << ", " << K.cols() << ") K" << std::endl;
			PnPSolver PnP = PnPSolver(points_3D, points_2D, K);
			PnP.Solve();
			if (last_point[3] == 0){
				Eigen::Matrix<double, 4, 4> sol = PnP.get_best_solution();
				Eigen::Matrix<double, 4, 1> camera_position = sol.inverse() * (Eigen::Matrix<double, 4, 1>()  << 0, 0, 0, 1).finished();
				camera_positions.push_back(camera_position);
				last_point = camera_position;
			}else{
				std::vector<Eigen::Matrix<double, 4, 4>> solutions = PnP.get_all_solutions();
				Eigen::Matrix<double, 4, 1> min_camera_position = solutions[0].inverse() * (Eigen::Matrix<double, 4, 1>()  << 0, 0, 0, 1).finished();
				double min_distance = pow(min_camera_position[0] - last_point[0], 2) +
															pow(min_camera_position[1] - last_point[1], 2) +
															pow(min_camera_position[2] - last_point[2], 2);
				for (int i=1; i<solutions.size(); i++){
					Eigen::Matrix<double, 4, 1> camera_position = solutions[i].inverse() * (Eigen::Matrix<double, 4, 1>()  << 0, 0, 0, 1).finished();
					double distance = pow(camera_position[0] - last_point[0], 2) +
														pow(camera_position[1] - last_point[1], 2) +
														pow(camera_position[2] - last_point[2], 2);
					if (distance < min_distance){
						min_camera_position = camera_position;
					}
				}
				std::cout << "Distance: " << min_distance << std::endl;
				camera_positions.push_back(min_camera_position);
				last_point = min_camera_position;
			}
			auto stop = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
			std::cout << "Time taken: " << duration.count()*10e-03 << " milliseconds" << std::endl << std::endl;
			std::cout << "--------------------------" << std::endl;
		}

    // Display the resulting frame
    cv::imshow("Frame", frame);
		video.write(frame);

    // Press  ESC on keyboard to exit
    char c=(char)cv::waitKey(25);
    if(c==27)
      break;
  }

  // When everything done, release the video capture object
  cap.release();
	video.release();

  // Closes all the frames
  cv::destroyAllWindows();

	std::ofstream output_file;
	output_file.open (OUTPUT_PATH);

	for (int i=0; i<camera_positions.size(); i++){
		output_file << camera_positions[i][0] << " " << camera_positions[i][1] << " " << camera_positions[i][2] << std::endl;
	}
	output_file.close();
	return 0;
	/*
		START PART 3
		There I detect the markers in the image getting all the corners 2D coordinates
		and I create a correspondencies with the 3D coorinate. The result are two vectors
		of equal size: points_3D, points_2D containing the corresponding 3D point and 2D
		point for the same marker (identifed by the id).
	*/

		/*
	for (int i=0; i<points_3D.size(); i++){
		std::cout << points_3D[i] << std::endl;
		std::cout << "--------------------------------" << std::endl;
	}
	for (int i=0; i<points_2D.size(); i++){
		std::cout << points_2D[i] << std::endl;
		std::cout << "--------------------------------" << std::endl;
	}
	*/
	/*
		FINISH PART 3
	*/

	/*
		START PART 4
		There I use the PnP solver that I have implemented and I print the results
	*/
	/*
		FINISH PART 4
	*/
}
