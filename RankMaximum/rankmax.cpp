#include <opencv2/opencv.hpp>
#include <iostream>
#include <algorithm>
#include <math.h>
#define K 3

//Rank Max Operator is one of the operator of the Rank operators family.
//The Rank Operator Family consisnt in a non-linear filtering of an image, involing
//a "sorting" of the calculated value of a certain local area. With different values taken,
//there are different result obtained. That's why we use a ones filter, because we want the exacts
//value from the original image to put into a vector to sort after the convolution is done.
//The rank max filter, produces a phenomenon called: DILATATION: it is the corrispective of the espansion
//on an instogram, with the result of "shifting" the grayscale value to the right of it, that means to the
//bright side of the image. 
//we can use the dilate operation for:
//°Removing noise
//°Isolation of individual elements and joining disparate elements in an image.
//°Finding of intensity bumps or holes in an image.
//The process of dilatation is given by, having a vector of values of size N, sort this vector
//and take the right-most element of it (at size-1). In this way we're increasing the luminosity
//into an image.


using namespace cv;
using namespace std;

Mat rank_max_filter = (Mat_<uchar>(3,3) << 1,1,1,
										1,1,1,
										1,1,1);
										
void padding(Mat, Mat&);
void dilate(Mat, Mat&);

int main(int argc, char** argv) {
	
	if(argc!=2) {
		cerr << "Usage: ./<programname> <imagename>.<filter>" << endl;
		exit(EXIT_FAILURE);
	}
	
	Mat raw_image = imread(argv[1], IMREAD_GRAYSCALE);
	
	if(raw_image.empty()) {
		cerr << "Image format not valid." << endl;
		exit(EXIT_FAILURE);
	}
	
	Mat padded_image, dest_image;
	
	padding(raw_image, padded_image);
	
	dilate(padded_image, dest_image);
	
	imshow("Original Image", raw_image);
	imshow("Dilated Image", dest_image);
	waitKey(0);
	return(0);

}

void padding(Mat raw_image, Mat& padded_image) {

	int pad = floor(K/2);
	padded_image = Mat::zeros(raw_image.rows+2*pad, raw_image.cols+2*pad, raw_image.type());
	
	for(int i=0; i<raw_image.rows; i++) {
		for(int j=0; j<raw_image.cols; j++) {
			padded_image.at<uchar>(i+pad, j+pad) = raw_image.at<uchar>(i,j);		
		}
	}

}

void dilate(Mat padded_image, Mat& dest_image) {

	int pad = floor(K/2);
	vector<int> dilated_values;
	dest_image = Mat::zeros(padded_image.rows-2*pad, padded_image.cols-2*pad, padded_image.type());
	
	for(int i=pad; i<padded_image.rows-2*pad; i++) {
		for(int j=pad; j<padded_image.cols-2*pad; j++) {
			for(int u=0; u<rank_max_filter.rows; u++) {
				for(int v=0; v<rank_max_filter.cols; v++) {
					dilated_values.push_back((int)padded_image.at<uchar>(i-pad+u, j-pad+v)*rank_max_filter.at<uchar>(u,v));
				}
			}
			
			sort(dilated_values.begin(), dilated_values.end());
			dest_image.at<uchar>(i-pad, j-pad) = dilated_values.at(dilated_values.size()-1);
			dilated_values.clear();
			
		}	
	}

}
