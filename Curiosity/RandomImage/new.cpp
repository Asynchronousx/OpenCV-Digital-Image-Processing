#include <opencv2/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <math.h>


using namespace cv;
using namespace std;


void create_random_image(Mat&, int, int);

int main() {

	srand(time(0));

	Mat gray_img;
	int height = 1024;
	int width  = 1024;
	

	gray_img=Mat::zeros(width, height, CV_8UC3);
	create_random_image(gray_img, width, height);
	namedWindow("Random IMG", WINDOW_AUTOSIZE);
	imshow("Random IMG", gray_img);
	waitKey(0);
	return(0);

}

void create_random_image(Mat& image, int w, int h) {
	int j;
	for(int i=0; i<w; i++) {
		for(j=0;j<h; j++) {
			image.at<Vec3b>(j,i)[0] = i+j*2;
			image.at<Vec3b>(i,j)[1] = i*20+pow(j,2);
			image.at<Vec3b>(j,j)[2] = i;		
		}
	}
}
