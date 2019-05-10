#include <opencv2/opencv.hpp>
#include <iostream>
#include <math.h>
#define K 3

//Note: if you was looking for theory, just look at the folder AverageBlur.
//Weighblur is the same of the Average Blur, but with "weighted values", that means
//that the mask is composed not by only ones but from different values. 
//In this special case, the weight_average_filter mask correspond to the following mask filter:
// |1 2 1| 
// |2 4 2| -> The sum of this filter is 16, then every value in this filter is divided by 16.
// |1 2 1|
//For example, 1/16 is 1/16. 2/16 is 1/8. 4/16 is 1/4. that's how we've calculated this filter.
//In this special case we're giving more importance to the center of the filter. That means if a noise pixel
//pass under the center, the noise will be amplificated. Consider a Weighted filter with a 0 center to reduce this problem.

using namespace std;
using namespace cv;

Mat weight_average_filter = (Mat_<float>(3,3) << 1/16.0, 1/8.0, 1/16.0,
												 1/8.0,  1/4.0, 1/8.0,
												 1/16.0, 1/8.0, 1/16.0);
												 
float calculate_total_value(Mat);
void padding(Mat, Mat&);
void blur(Mat, Mat&);

int main(int argc, char** argv) {

	if(argc!=2) {
		cerr << "Usage: ./<programname> <imagename>.<format>" << endl;
		exit(EXIT_FAILURE);
	}
	
	Mat raw_image = imread(argv[1], IMREAD_GRAYSCALE);
	
	if(raw_image.empty()) {
		cerr << "Image format not valid." << endl;
		exit(EXIT_FAILURE);
	}	

	Mat padded_image, dest_image;
	
	padding(raw_image, padded_image);
	blur(padded_image, dest_image);
	
	imshow("Original Image", raw_image);
	imshow("Blurred Image", dest_image);
	waitKey(0);
	return(0);
	
}

float calculate_total_value(Mat filter) {

	float value = 0;
	
	for(int i=0; i<filter.rows; i++) {
		for(int j=0; j<filter.cols; j++) {
			value += filter.at<float>(i,j);
		}
	}
	
	return value;
	
}

void padding(Mat raw_image, Mat& padded_image) {

	int pad = floor(K/2);
	padded_image = Mat::zeros(raw_image.rows+2*pad, raw_image.cols+2*pad, raw_image.type());
	
	for(int i=0; i<raw_image.rows; i++) {
		for(int j=0; j<raw_image.rows; j++) {
			padded_image.at<uchar>(i+pad, j+pad) = raw_image.at<uchar>(i,j);
		}
	}

}

void blur(Mat padded_image, Mat& dest_image) {
	
	float pixel_intensity = 0;
	int pad = floor(K/2);
	int total_value = calculate_total_value(weight_average_filter);
	dest_image = Mat::zeros(padded_image.rows-2*pad, padded_image.cols-2*pad, padded_image.type());
	
	for(int i=pad; i<padded_image.rows-2*pad; i++) {
		for(int j=pad; j<padded_image.cols-2*pad; j++) {
			for(int u=0; u<weight_average_filter.rows; u++) {
				for(int v=0; v<weight_average_filter.cols; v++) {
					pixel_intensity += padded_image.at<uchar>(i-pad+u, j-pad+v)*weight_average_filter.at<float>(u,v);
				}
			}
					
			dest_image.at<uchar>(i-pad, j-pad) = floor(pixel_intensity/total_value);
			pixel_intensity = 0;
		
		}
	}
	
}
