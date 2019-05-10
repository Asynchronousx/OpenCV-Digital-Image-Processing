#include <opencv2/opencv.hpp>
#include <iostream>
#include <string>
#include <math.h>
#define L 256
#define Ksize 3

using namespace std;
using namespace cv;

Mat Robert_Horizontal = (Mat_<char>(3,3) << 0, 0, 0,
									     0, 0, 1,
									     0, -1, 0);
									   
Mat Robert_Vertical = (Mat_<char>(3,3) << 0, 0, 0,
									   0, 1, 0,
									   0, 0, -1);
									  
									  
void padding(Mat, Mat&);
void RobertOperator(Mat, Mat, Mat&);

int main(int argc, char** argv) {

	if(argc!=2) {
		cerr << "Usage: ./<programname> <imagename>.<format>" << endl;
		exit(EXIT_FAILURE);
	}
	
	Mat original_image = imread(argv[1], IMREAD_COLOR);
	Mat raw_image = imread(argv[1], IMREAD_GRAYSCALE);
	
	if(raw_image.empty()) {
		cerr << "Usage ./<programname> <imagename>.<format>" << endl;
		exit(EXIT_FAILURE);
	}
	
	Mat dest_image = Mat::zeros(raw_image.rows, raw_image.cols, raw_image.type());
	Mat padded_image;
	
	padding(raw_image, padded_image);
	RobertOperator(raw_image, padded_image, dest_image);
	
	string name(argv[1]);
	
	imshow("Original Image", original_image);
	imshow("Raw Image", raw_image);
	imshow("Robert Image", dest_image);
	imwrite("RobertConv"+name, dest_image);
	waitKey(0);
	return(0);

}

void padding(Mat raw_image, Mat& padded_image) {
	
	int pad = floor(Ksize/2);
	padded_image = Mat::zeros(raw_image.rows+2*pad, raw_image.cols+2*pad, raw_image.type());
	
	for(int i=0; i<raw_image.rows; i++) {
		for(int j=0; j<raw_image.cols; j++) {
			padded_image.at<uchar>(i+pad, j+pad) = raw_image.at<uchar>(i,j);
		}
		
	}
 
}

void RobertOperator(Mat raw_image, Mat padded_image, Mat& dest_image) {

	int pixel_intensity = 0;
	int gradientX = 0;
	int gradientY = 0;
	int pad = floor(Ksize/2);
	
	for(int i=pad; i<padded_image.rows-pad; i++) {
		for(int j=pad; j<padded_image.cols-pad; j++) {
			for(int u=0; u<Ksize; u++) {
				for(int v=0; v<Ksize; v++) {
					gradientX += padded_image.at<uchar>(i-pad+u, j-pad+v)*Robert_Horizontal.at<char>(u,v);
					gradientY += padded_image.at<uchar>(i-pad+u, j-pad+v)*Robert_Vertical.at<char>(u,v);
				}
			}
			
			//do the math
			pixel_intensity = abs(gradientX) + abs(gradientY);
			pixel_intensity = pixel_intensity > 255? 255: pixel_intensity < 0? 0: pixel_intensity;
			dest_image.at<uchar>(i-pad,j-pad) = pixel_intensity;
			gradientX = 0;
			gradientY = 0;
			
		}
	}
	
}
