#include <opencv2/opencv.hpp>
#include <iostream>
#include <math.h>
#define L 256

using namespace std;
using namespace cv;

void PrewittOperator(Mat, Mat&);

int main(int argc, char** argv) {

	if(argc != 2) {
		cerr << "Usage: ./<programname> <imagename>.<format>" << endl;
		exit(EXIT_FAILURE);
	}
	
	Mat raw_image = imread(argv[1], IMREAD_GRAYSCALE);
	
	if(raw_image.empty()) {
		cerr << "Image data invalid. Try again." << endl;
		exit(EXIT_FAILURE);	
	}
	
	Mat dest_image = Mat::zeros(raw_image.rows, raw_image.cols, raw_image.type());
	
	PrewittOperator(raw_image, dest_image);
	
	imshow("Original Image", raw_image);
	imshow("Prewitt Image", dest_image);
	waitKey(0);
	return(0);

}

//Sobel Operator is one of the Edge Calculator for an image given a convolution Matrix 3x3. Sobel is similar to Prewitt and Isotropic operator.
//We know that the Gradient is a "substitute" for the derivates of a function in a bidimensional space.
//The Gradient is composed by two element, Partial Derivate of X, and Partial Derivate of Y.
//The derivates calculates how the function is going to change and so do the Gradient: it POINT to the function changing verse.
//That said, the gradient is composed by two element because we're in two dimension: partial derivate of x and y.
//and we also knwo that the edges pixel are ORTOGONAL to the gradient.
//Sobel matrix is composed like that:
//	-1  0  1	       					    1  2  1
//  -2  0  2  identificable as Gx(x,y) and  0  0  0 identificable as Gy(x,y). 
//  -1  0  1							   -1 -2 -1
//Given those two 3x3 matrix, we need to convolute them on the original image.
//We can also avoid convolution by simply adding and subtracting values between them, avoiding the two additional for needed to convolute
//the matrix on the image.
void PrewittOperator(Mat raw_image, Mat& dest_image) {

	//Here, gradient x and y represent the two component of our Gradient. 
	//pixel_intensity represent how the pixel vary on the grayscale (0-255).
	int pixel_intensity = 0;
	int gradientX, gradientY;
	
	//for the size of the raw image
	for(int i=0; i<raw_image.rows; i++) {
		for(int j=0; j<raw_image.cols; j++) {
			gradientX = ((-1*raw_image.at<uchar>(i,j)) + (-2*raw_image.at<uchar>(i+1,j)) + (-1*raw_image.at<uchar>(i+2,j))) + 
						((1*raw_image.at<uchar>(i,j+2)) + (2*raw_image.at<uchar>(i+1,j+2)) + (1*raw_image.at<uchar>(i+2,j+2)));
			gradientY = ((1*raw_image.at<uchar>(i,j)) + (1*raw_image.at<uchar>(i,j+1)) + (1*raw_image.at<uchar>(i,j+2))) + 
						((-1*raw_image.at<uchar>(i+2,j)) + (-1*raw_image.at<uchar>(i+2,j+1)) + (-1*raw_image.at<uchar>(i+2,j+2)));
			pixel_intensity = abs(gradientX + gradientY);
			pixel_intensity = pixel_intensity > 255? 255: pixel_intensity < 0? 0: pixel_intensity;
			dest_image.at<uchar>(i,j) = pixel_intensity;
		
		}
	}

}

