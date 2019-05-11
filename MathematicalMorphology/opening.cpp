#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

//What is opening? The opening is a special mathematic morphologic operation that have the goal
//of removing small object from an image (such as pixel that doesn't have nothing to do with the coherency
//of an image), or smoothing the angle of the connected components present into the image.
//The opening is composed by an erosion of an image with his structuring elements (4/8 connected), 
//followed by a dilation with the same structuring elements.
//So we can say that: 
//O(A,K) = (A(-)K)(+)K

void Opening(Mat, Mat&, Mat&);

int main(int argc, char** argv) {

	if(argc!=2) {
		cerr << "Usage: ./<programname> <imagename.format>" << endl;
		exit(EXIT_FAILURE);
	}
	
	Mat raw_image = imread(argv[1]);
	
	if(raw_image.empty()) {
		cerr << "Image format not valid." << endl;
		exit(EXIT_FAILURE);
	}

	Mat dest_image_combination = Mat::zeros(raw_image.size(), raw_image.type());
	Mat dest_image_function = dest_image_combination.clone();
	
	Opening(raw_image, dest_image_combination, dest_image_function);
	
	waitKey(0);
	return 0;
	
}

void Opening(Mat raw_image, Mat& dest_image_combination, Mat& dest_image_function) {

	//The opening could be achieved in two ways: with the application of erode and dilate,
	//Or with the function morphologyEx.
	//declaring a tmp dest image to hold the middle result of the opening and a structuring element
	Mat dest_combination_tmp = dest_image_combination.clone();
	Mat s_element = getStructuringElement(MORPH_RECT, Size(3,3), Point(-1, -1));	
	//1) Combination of erode and dilate
	//int iteration = N;
	erode(raw_image, dest_combination_tmp, s_element); //Point(-1,-1), iterations);
	
	//Dilating the above result to obtain the opened image 
	dilate(dest_combination_tmp, dest_image_combination, s_element); //Point(-1,-1), iterations);
	
	//2) Using the function 
	//Using the function morphologyEx
	morphologyEx(raw_image, dest_image_function, MORPH_OPEN, s_element); //, Point(-1,-1), iterations);

	//show the images
	imshow("Original", raw_image);
	imshow("Opened Manual", dest_image_combination);
	imshow("Opened Function", dest_image_function);
	waitKey(0);

}
