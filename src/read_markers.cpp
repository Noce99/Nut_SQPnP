#include <opencv2/aruco.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/objdetect/aruco_detector.hpp>
#include <iostream>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string> 

int main(){
	cv::aruco::DetectorParameters detectorParams = cv::aruco::DetectorParameters();
	cv::aruco::Dictionary dictionary = cv::aruco::getPredefinedDictionary(cv::aruco::DICT_APRILTAG_16h5);
	cv::aruco::ArucoDetector detector(dictionary, detectorParams);
	cv::Mat image = cv::imread("./images/real_really_small_marker.jpg");
	cv::Mat imageCopy;
	image.copyTo(imageCopy);
	
	std::vector<int> ids;
	std::vector<std::vector<cv::Point2f>> corners, rejected;
	detector.detectMarkers(image, corners, ids, rejected);
	
	for (int i=0; i<ids.size(); i++){
		for (int ii=0; ii<4; ii++){
			std::cout << corners[i][ii] << " ";
		}
		std::cout << std::endl;
	}
	
	if (ids.size() > 0)
		cv::aruco::drawDetectedMarkers(imageCopy, corners, ids);
	cv::imwrite("trovato.jpg", imageCopy);
}
