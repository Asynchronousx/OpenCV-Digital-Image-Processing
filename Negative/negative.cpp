#include <opencv2/opencv.hpp>
#include <iostream>
#define L 256

using namespace cv;
using namespace std;

int main(int argc, char** argv) {

	if(argc!=2) {
		cerr << "Usage: ./programname <imagename>.<format>" << endl;
		exit(EXIT_FAILURE);
	}
	
	Mat image = imread(argv[1], IMREAD_GRAYSCALE);
	
	if(image.empty()) {
		cerr << "Invalid image. Try again with another one." << endl;
		exit(EXIT_FAILURE);
	}
	
	for(int i=0; i<image.rows; i++) {
		for(int j=0; j<image.cols; j++) {
			image.at<uchar>(i,j) = L - 1 - image.at<uchar>(i,j);		
		}
	}
	
	imshow("Negative Image", image);
	waitKey(0);
	return(0);

}

/* RGB VERSION 

Mat image = imread(argv[1], IMREAD_COLOR);
	
	for(int i=0; i<image.rows; i++) {
		for(int j=0; j<image.cols; j++) {
			image.at<Vec3b>(i,j)[0] = L - 1 - image.at<Vec3b>(i,j)[0];	
			image.at<Vec3b>(i,j)[1] = L - 1 - image.at<Vec3b>(i,j)[1];
			image.at<Vec3b>(i,j)[2] = L - 1 - image.at<Vec3b>(i,j)[2];	
		}
	}
	
*/
