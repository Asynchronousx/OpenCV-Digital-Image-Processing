#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <opencv2/opencv.hpp>
#include <stdlib.h>
#include <stdio.h>

using namespace cv;


int main( int argc, char** argv )
{

	const String window_name = "Sobel Demo - Simple Edge Detector";
    int ksize = 1;
    int scale = 1;
    int delta = 0;
	int ddepth = CV_16S;
	Mat grad;

	Mat src_gray = imread(argv[1], IMREAD_GRAYSCALE);
    //![sobel]
    /// Generate grad_x and grad_y
    Mat grad_x, grad_y;
    Mat abs_grad_x, abs_grad_y;

    /// Gradient X
    Sobel(src_gray, grad_x, ddepth, 1, 0, ksize, scale, delta, BORDER_DEFAULT);

    /// Gradient Y
    Sobel(src_gray, grad_y, ddepth, 0, 1, ksize, scale, delta, BORDER_DEFAULT);
    //![sobel]

    //![convert]
    // converting back to CV_8U
    convertScaleAbs(grad_x, abs_grad_x);
    convertScaleAbs(grad_y, abs_grad_y);
    //![convert]

    //![blend]
    /// Total Gradient (approximate)
    addWeighted(abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad);
    //![blend]

    //![display]
	imshow(window_name, grad);
	waitKey(0);
    imwrite("sobel.png", grad);
 	 
  
 }
