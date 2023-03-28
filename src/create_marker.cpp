#include <opencv2/aruco.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/objdetect/aruco_detector.hpp>
#include <iostream>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string> 

int main(){
	cv::Mat markerImage;
	cv::aruco::Dictionary dictionary = cv::aruco::getPredefinedDictionary(cv::aruco::DICT_APRILTAG_16h5);
	mkdir("markers",0777);
	for (int i=0; i<20; i++){
		cv::aruco::generateImageMarker(dictionary, i, 500, markerImage, 1);
		cv::imwrite("./markers/marker"+std::to_string(i)+".png", markerImage);
	}
}
