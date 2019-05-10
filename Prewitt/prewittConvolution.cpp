#include <opencv2/opencv.hpp>
#include <iostream>
#include <math.h>
#define L 256
#define Ksize 3

using namespace std;
using namespace cv;

//Defining a Prewitt Matrix
Mat prewitt_horizontal = (Mat_<char>(3,3) << -1, 0, 1,
											 -1, 0, 1,
											 -1, 0, 1);

Mat prewitt_vertical = (Mat_<char>(3,3) <<  1,  1,  1, 
										    0,  0,  0,
										   -1, -1, -1);

void padding(Mat, Mat&);
void PrewittOperator(Mat, Mat, Mat&);

int main(int argc, char** argv) {

	if(argc != 2) {
		cerr << "Usage: ./<programname> <imagename>.<format>" << endl;
		exit(EXIT_FAILURE);
	}
	
	Mat original_image = imread(argv[1], IMREAD_COLOR);
	Mat raw_image = imread(argv[1], IMREAD_GRAYSCALE);
	
	if(raw_image.empty()) {
		cerr << "Image data invalid. Try again." << endl;
		exit(EXIT_FAILURE);	
	}
	
	Mat padded_image;
	Mat dest_image = Mat::zeros(raw_image.rows, raw_image.cols, raw_image.type());
	
	padding(raw_image, padded_image);
	PrewittOperator(raw_image, padded_image, dest_image);

	string name(argv[1]);
	
	imshow("Colored Image", original_image);
	imshow("Original Image", raw_image);
	imshow("Prewitt Image", dest_image);
	imwrite("PrewittConv"+name, dest_image);
	waitKey(0);
	return(0);
	
}

//Function to create a zero padding layer to the dest image
void padding(Mat raw_image, Mat& padded_image) {
	
	int pad = floor(Ksize/2);
	
	//fill the matrix 
	padded_image = Mat::zeros(raw_image.rows+2*pad, raw_image.cols+2*pad, raw_image.type());
	
	//copy the image into the padded one without the borders
	for(int i=0; i<raw_image.rows; i++) {
		for(int j=0; j<raw_image.cols; j++) {
			padded_image.at<uchar>(i+pad, j+pad) = raw_image.at<uchar>(i,j);
		}
	}

}

//Prewitt Operator is one of the Edge Calculator for an image given a convolution Matrix 3x3. Prewitt is similar to Sobel and Isotropic operator.
//We know that the Gradient is a "substitute" for the derivates of a function in a bidimensional space.
//The Gradient is composed by two element, Partial Derivate of X, and Partial Derivate of Y.
//The derivates calculates how the function is going to change and so do the Gradient: it POINT to the function changing verse.
//That said, the gradient is composed by two element because we're in two dimension: partial derivate of x and y.
//and we also knwo that the edges pixel are ORTOGONAL to the gradient.
//Prewitt matrix is composed like that:
//	-1  0  1	       					    1  1  1
//  -1  0  1  identificable as Gx(x,y) and  0  0  0 identificable as Gy(x,y). 
//  -1  0  1							   -1 -1 -1
//Given those two 3x3 matrix, we need to convolute them on the original image.
void PrewittOperator(Mat raw_image, Mat padded_image, Mat& dest_image) {

	//Here, gradient x and y represent the two component of our Gradient. 
	//pixel_intensity represent how the pixel vary on the grayscale (0-255).
	int pixel_intensity = 0;
	int gradientX = 0;
	int gradientY = 0;
	int pad = floor(Ksize/2);
	
	Mat dx = (Mat_<uchar>(raw_image.rows, raw_image.cols));
	Mat dy = (Mat_<uchar>(raw_image.rows, raw_image.cols));

	for(int i=pad; i<padded_image.rows-pad; i++) {
		for(int j=pad; j<padded_image.cols-pad; j++) {
			for(int u=0; u<Ksize; u++) {
				for(int v=0; v<Ksize; v++) {
					gradientX += padded_image.at<uchar>(i-pad+u, j-pad+v)*prewitt_horizontal.at<char>(u,v);
					gradientY += padded_image.at<uchar>(i-pad+u, j-pad+v)*prewitt_vertical.at<char>(u,v);
				}
			}
			
			//do the math
			pixel_intensity = abs(gradientX) + abs(gradientY);
			pixel_intensity = pixel_intensity > 255? 255: pixel_intensity < 0? 0: pixel_intensity;
			dest_image.at<uchar>(i-pad,j-pad) = pixel_intensity;
			double theta = atan2(gradientX, gradientY);
			double degree = theta*180/M_PI;
			
			dx.at<uchar>(i-pad, j-pad) = gradientX;
			dy.at<uchar>(i-pad, j-pad) = gradientY;
		
			//Once here, reset the gradientx and y because we're considering another pixel.
			gradientX = 0;
			gradientY = 0;
			 
		}
	
	}
	
	
	imshow("DX", dx);
	imshow("DY", dy);
	waitKey(0);
	

}

