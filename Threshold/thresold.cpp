#include <opencv2/opencv.hpp>
#define M 128

using namespace std;
using namespace cv;

int main(int argc, char** argv) {
	if(argc!=2) {
		cerr << "Must provide an image into the terminal." << endl;
		exit(EXIT_FAILURE);
	}

	Mat image = imread(argv[1], IMREAD_GRAYSCALE);
	
	if(image.empty()) {
		cerr << "Image format is invalid." << endl;
		exit(EXIT_FAILURE);
	}

	//From the theory, a good M should be found analyzing the histogram and taking the two-most occurrence
	//of the grayscale: I.e -> value 23 of the grayscale with value 223 and value 120 of the grayscale
	//with value of 219.
	//Have those two point, let's call them Xa and Xb, we need to get the average out of them.
	//So, Xa-Xb/2. The result is our M.
	
	for(int i=0; i<image.rows; i++) {
		for(int j=0; j<image.cols; j++) {
			image.at<uchar>(i,j) = image.at<uchar>(i,j) < M? 0: 255;
		}
	}
	
	imshow("Display Image", image);
	waitKey(0);
	return(0);
	

}

/* RGB VERSION

Mat image = imread(argv[1], IMREAD_COLOR);
	
	if(image.empty()) {
		cerr << "Image format is invalid." << endl;
		exit(EXIT_FAILURE);
	}

	
	for(int i=0; i<image.rows; i++) {
		for(int j=0; j<image.cols; j++) {
			image.at<Vec3b>(i,j)[0] = image.at<Vec3b>(i,j)[0] < M? 0: 255;
			image.at<Vec3b>(i,j)[1] = image.at<Vec3b>(i,j)[1] < M? 0: 255;
			image.at<Vec3b>(i,j)[2] = image.at<Vec3b>(i,j)[2] < M? 0: 255;	
		}
	}
	
*/
