#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

//Dilatation is the process which will "increase" or "dilate" the pixel present into an image.
//What means "dilate"? For explaining the concept, we can refer to the Minkowski sum.
//Having an image A and a structuring element B, we can say that everytime the structuring element 
//B, center his anchor pixel in a pixel of A, it will replicate this pixel into all of his active elements.
//We can say that, the new pixel c is given c = a+b with a€A, b€B.

void Dilate(Mat, Mat&);

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

	Mat dest_image = Mat::zeros(raw_image.size(), raw_image.type());

	Dilate(raw_image, dest_image);
	
	imshow("Original", raw_image); 
	imshow("Dilated", dest_image);
	waitKey(0);
	return 0;
	
}


void Dilate(Mat raw_image, Mat& dest_image) {

	//The dilate function works exatcly like the erode function. For detailed explanation,
	//check out the erode.cpp file.
	//Getting the structured element
	Mat s_element = getStructuringElement(MORPH_RECT, Size(3,3), Point(-1,-1));
	
	//Calling the dilate function
	int iterations = 1;
	dilate(raw_image, dest_image, s_element, Point(-1,-1), iterations); 

}
