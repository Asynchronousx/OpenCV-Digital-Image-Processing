#include <opencv2/opencv.hpp>
#include <iostream>
#include <math.h>
#define L 256

using namespace std;
using namespace cv;

void RobertOperator(Mat, Mat&);

int main(int argc, char** argv) {

	if(argc != 2) {
		cerr << "Usage: ./<programname> <imagename>.<format>" << endl;
		exit(EXIT_FAILURE);
	}
	
	Mat raw_image = imread(argv[1], IMREAD_GRAYSCALE);
	
	if(raw_image.empty()) {
		cerr << "Image data not valid." << endl;
		exit(EXIT_FAILURE);
	}
	
	Mat dest_image = Mat::zeros(raw_image.rows, raw_image.cols, raw_image.type());

	RobertOperator(raw_image, dest_image);
	
	string name(argv[1]);
	
	imshow("Original Image", raw_image);
	imshow("Robert Edge", dest_image);
	imwrite("RobertMul"+name, dest_image);
	waitKey(0);
	return(0);

}


//Robert Operator is one of the most simple, thus ambigue operator to detect an edge into an image.
//We know that the Gradient is a "substitute" for the derivates of a function in a bidimensional space.
//The Gradient is composed by two element, Partial Derivate of X, and Partial Derivate of Y.
//The derivates calculates how the function is going to change and so do the Gradient: it POINT to the function changing verse.
//That said, the gradient is composed by two element because we're in two dimension: partial derivate of x and y.
//and we also knwo that the edges pixel are ORTOGONAL to the gradient.
//Robert matrix is composed like that:
//	0  1							   1  0
// -1  0 identificable as Gx(x,y) and  0 -1 identificable as Gy(x,y).
//Given those two 2x2 matrix, we need to convolute them on the original image.
//We can also avoid convolution by simply adding and subtracting values between them, avoiding the two additional for needed to convolute
//the matrix on the image.
void RobertOperator(Mat raw_image, Mat& dest_image) {

	//Here, gradient x and y represent the two component of our Gradient. 
	//pixel_intensity represent how the pixel vary on the grayscale (0-255).
	int gradientX, gradientY;
	int pixel_intensity;
	
	//for the size of the entire image
	for(int i=0; i<raw_image.rows; i++) {
		for(int j=0; j<raw_image.cols; j++) {
			//now we need to calculate the two gradient by simulating a convolution on the image pixel:
			//since half of the matrix is composed by 0 (the diagonal) we can just do two sum.
			gradientX = ((1*raw_image.at<uchar>(i,j+1)) + (-1*raw_image.at<uchar>(i+1, j)));
			gradientY = ((1*raw_image.at<uchar>(i,j+1))   + (-1*raw_image.at<uchar>(i+1, j)));
			
			pixel_intensity = abs(gradientX) + (gradientY);
			
			//check for a wrong magnitude
			pixel_intensity = pixel_intensity > 255? 255: pixel_intensity < 0? 0: pixel_intensity;
			
			//build the image at that pixel
			dest_image.at<uchar>(i,j) = pixel_intensity;			
		}
	}

}
