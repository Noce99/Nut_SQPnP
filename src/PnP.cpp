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
#define IMAGE_PATH "/home/enrico/Progetti/Nut_SQPnP/resources/test_image.jpg"

int main () {

	auto start = std::chrono::high_resolution_clock::now();

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
  }else std::cout << "Unable to open file";
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
	cv::Mat image = cv::imread(IMAGE_PATH);

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

	cv::Mat imageUndistorted;
	cv::undistort(image, imageUndistorted, Camera_Matrix, Distortion_Coefficients);
	// cv::imwrite("./images/real_really_small_marker_1_u.jpg", imageUndistorted);
	/*
		FINISH PART 2
	*/

	/*
		START PART 3
		There I detect the markers in the image getting all the corners 2D coordinates
		and I create a correspondencies with the 3D coorinate. The result are two vectors
		of equal size: points_3D, points_2D containing the corresponding 3D point and 2D
		point for the same marker (identifed by the id).
	*/
	std::vector<int> ids;
	std::vector<std::vector<cv::Point2f>> corners, rejected;
	detector.detectMarkers(imageUndistorted, corners, ids, rejected);

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

	std::cout << "I have [" << points_3D.size() << "] 3D points, [" << points_2D.size() << "] 2D points and a (" << K.rows() << ", " << K.cols() << ") K" << std::endl;
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
	PnPSolver PnP = PnPSolver(points_3D, points_2D, K);
	PnP.Solve();
	Eigen::Matrix<double, 4, 4> sol = PnP.get_best_solution();
	std::cout << "Best solution found:" << std::endl;
	std::cout << sol << std::endl;
	std::cout << "--------------------------" << std::endl;
	std::cout << sol.inverse() << std::endl;
	std::cout << "--------------------------" << std::endl;
	std::cout << sol.inverse() * (Eigen::Matrix<double, 4, 1>()  << 0, 0, 0, 1).finished() << std::endl;
	std::cout << "--------------------------" << std::endl;
	std::cout << sol.inverse() * (Eigen::Matrix<double, 4, 1>()  << 0, 0, 10, 1).finished() << std::endl;

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "Time taken: " << duration.count()*10e-03 << " milliseconds" << std::endl << std::endl;
  return 0;
	/*
		FINISH PART 4
	*/
}
