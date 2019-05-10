#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <iostream>
#include <string>
using namespace cv;
using namespace std;
int main( int argc, char** argv )
{
	if(argc!=2) {
		cerr << "Usage: ./program image.format" << endl;
		exit(EXIT_FAILURE);
	} 
  
	Mat image = imread(argv[1], IMREAD_GRAYSCALE);
	
	if(image.empty()) {
		cerr << "Data not valid. Try with an image." << endl;
		exit(EXIT_FAILURE);
	}	

	namedWindow("Display IMG", WINDOW_AUTOSIZE);
	imshow("Display IMG", image);
	waitKey(0);
	return(0);

}
