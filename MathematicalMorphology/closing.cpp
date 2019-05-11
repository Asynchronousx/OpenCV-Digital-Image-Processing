#include <opencv2/opencv.hpp>

//What is Closing? Closing is the process that, as the name suggests, will take care
//of all the holes present into an image filling them.
//The closing operation will take care of all the interconnected objects into an image
//that present a gap or hole in-between, closing those connection edges and resulting in a
//restored image.
//We can obtain a closing by dilating an image A on a filter (4/8 connected) and eroding the result
//of this operation.
//We can say that:
//C(A,K) = (A(+)K)(-)K 

using namespace std;
using namespace cv;

void Closing(Mat, Mat&, Mat&);

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
	
	Closing(raw_image, dest_image_combination, dest_image_function);
	
	waitKey(0);
	return 0;
	
}

void Closing(Mat raw_image, Mat& dest_image_combination, Mat& dest_image_function) {

	//We're going to obtain the closing in two ways: with the combination of erosion and dilation,
	//And from the function morphologyEx.
	//Let's start declaring a TMP dest combination mat object that will be useful when doing the middle
	//step of the closing, then declare the structuring element.
	Mat dest_combination_tmp = dest_image_combination.clone();
	Mat s_element = getStructuringElement(MORPH_RECT, Size(3,3), Point(-1,-1));
	
	//1) Obtain the closing by combination of erode A(-)K and dilation -> (+)K
	//int iterations = 3;
	dilate(raw_image, dest_combination_tmp, s_element); // Point(-1,-1), iterations);
	erode(dest_combination_tmp, dest_image_combination, s_element); // Point(-1,-1), iterations);
	
	//2) Obtain the closing by the function morphologyEx
	morphologyEx(raw_image, dest_image_function, MORPH_CLOSE, s_element); //Point(-1,-1), iterations);
	
	//Show the images
	imshow("Original", raw_image);
	imshow("Closed Manual", dest_image_combination);
	imshow("Closed Function", dest_image_function);
	waitKey(0);


}
