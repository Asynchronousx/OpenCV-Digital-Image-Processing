#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

//The Erosion is the procedure which will "eat" all the pixel of a given image if some rules are not met.
//From the Minkowski difference, given an structural element B, doing the traslation of this element on the 
//image A, the erosion will keep the pixel of the original image A which will be present into the structural mask
//B when it pass on the local region area of the image A. if all the pixel of the ACTIVE region belongs to A, 
//then those pixel are saved. If one or more of those active pixel of B does not contain any pixel of A, 
//the relative pixel of A got "erosed", or simply deleted. 

void Erode(Mat, Mat&);

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

	Erode(raw_image, dest_image);
	
	imshow("Original", raw_image); 
	imshow("Eroded", dest_image);
	waitKey(0);
	return 0;
	
}

void Erode(Mat raw_image, Mat& dest_image) {

	//We're using opencv erode: first of all, we need to get a structural element.
	//The function ERODE, if not specified, will initialize a default structural element of size
	//3x3, anchored in the center (Point(-1, -1) interpreted by the function as the center) 8-connected.
	//For the sake of the example, let's initialize our structuring element. 
	
	//We want our element to be 8Connected (rectangular), of size 3x3, anchored in the center.
	Mat s_element = getStructuringElement(MORPH_RECT, Size(3,3), Point(-1,-1));
	
	//Now let's erode the image with this element, then save the result into the dest image.
	//The function erode will take in input the raw image, the dest image, the structuring element, the anchor
	//of the kernel and the iterations to be done on the image.
	//if not specified, they will be set to 1
	int iterations = 1;
	erode(raw_image, dest_image, s_element, Point(-1, -1), iterations);

}
