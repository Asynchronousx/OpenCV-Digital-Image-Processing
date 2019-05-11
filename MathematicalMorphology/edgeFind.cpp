#include <opencv2/opencv.hpp>

//Extracting the edge given a binary image, is the equivalent of finding the morphological gradient
//of an image.
//We can achieve this operation by simply having an image A and his structuring element, B. 
//Then we should do the erosion of the image A with the element B, and the subtracting the result to the 
//original image A. this would give us the edge of the binary shape. (eroding an image will produce an edge-divoured
//image, then subtracting this image to the original one will give us the little pixels edge that the erosion
//cancelled, giving us the edges shape).
//We can then say: 
//E(A,K) = A - (A(-)K)

using namespace std;
using namespace cv;

void FindEdge(Mat, Mat&, Mat&);

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
	
	FindEdge(raw_image, dest_image_combination, dest_image_function);
	
	waitKey(0);
	return 0;
	
}

void FindEdge(Mat raw_image, Mat& dest_image_combination, Mat& dest_image_function) {

	//We're going to obtain the edge in two ways: with the combination of erosion and subtraction,
	//And from the function morphologyEx.
	//Let's start declaring a TMP dest combination mat object that will be useful when doing the middle
	//step of the closing, then declare the structuring element.
	Mat dest_combination_tmp = dest_image_combination.clone();
	Mat s_element = getStructuringElement(MORPH_RECT, Size(3,3), Point(-1,-1));
	
	//1) Obtaining the edge through the combination of erosion and subtraction
	erode(raw_image, dest_combination_tmp, s_element); //, Point(-1,-1), iterations);
	dest_image_combination = raw_image - dest_combination_tmp;
	
	//2) Obtaining the edge through the function morphologyEx
	morphologyEx(raw_image, dest_image_function, MORPH_GRADIENT, s_element); //, Point(-1,-1), iterations);
	
	//Show the images
	imshow("Original", raw_image);
	imshow("Edge Manual", dest_image_combination);
	imshow("Edge Function", dest_image_function);
	waitKey(0);
	


}
