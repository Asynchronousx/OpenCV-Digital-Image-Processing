#include <opencv2/opencv.hpp>
#include <iostream>

using namespace std;
using namespace cv;

int main(int argc, char** argv) {

	if(argc != 2) {
		cerr << "Usage: ./<programname> <imagename>.<format>" << endl;
		exit(EXIT_FAILURE);
	}
	
	Mat raw_image = imread(argv[1], IMREAD_GRAYSCALE);
	Mat dest_image = Mat::zeros(raw_image.rows, raw_image.cols, raw_image.type());
	
	if(raw_image.empty()) {
		cerr << "Image must be valid." << endl;
		exit(EXIT_FAILURE);
	}
	
	//What the normalize function do, is called "Constrast Stretching. 
	//Given an Alpha, and a Beta, Where:
	//Alpha is the minimum value of our interval range,
	//Beta is the maximum value of our interval range,
	//It stretches to left and right the histogram values to cover all the histogram range.
	//This increase the constrast. 
	//For example, if the minimum value of the histogram of the raw images was 20, and the maximum 200
	//with the normalization we say that 20 (Alpha) should be 0, and 200 should be 255. ALl the value 
	//in the middle got normalized to this new range scale.
	//This function takes in input the:
	//SOURCE image, the DESTINATION image, ALPHA value, BETA value, NORMALIZATION TYPE and the RESULTANT IMAGE TYPE 
	//(if not specified, default is the src image type).
	//The NORM_MINMAX is the normalization type and tell us that the minimum value is alpha, and the maximum is beta. 
	normalize(raw_image, dest_image, 0, 255, NORM_MINMAX, CV_8UC1);
	
	imshow("Original Image", raw_image);
	imshow("Normalized Image", dest_image); 
	waitKey(0);
	return(0);
	
}
