#include <opencv2/opencv.hpp>
#include <iostream>

using namespace std;
using namespace cv;

int main(int argc, char** argv) {

	Mat raw_image, dest_image;
	raw_image = imread(argv[1], IMREAD_GRAYSCALE);
	dest_image = Mat::zeros(raw_image.rows, raw_image.cols, raw_image.type());
	
	
	GaussianBlur( raw_image, dest_image, Size(17, 17), 0, 0 );
		
	imshow("Original Image", raw_image);
	imshow("Filtered Image", dest_image);
	imwrite("ng.png", dest_image);
	waitKey(0);
}
