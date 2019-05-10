#include <opencv2/opencv.hpp> 
#include <iostream>
#include <algorithm>
#include <math.h>
#define K 3

//Program to remove a noise with a median Filter.
//From the theory, we know that, took a local area of pixel MxN, we would have a series of values
//that represent the pixels. Sorted this local area (i.e: from (0,0) to (0,9) -> 9 pixel elements sorted from their value)
//for example (233, 10, 20, 94, 200, 190, 22, 3 84) -> (3, 10, 20, 22, 84, 94, 190, 200, 233) we now got a sorted array of local pixels.
//What is noise? Noise is something that doesn't is congruent with the image itself. In this case, what would be a non-congruent values? the ones
//at the perifery of the image. Infact, the median represent the medium value of the local area, usually the one that represent the REAL value of the area.
//SO, instead of making an average (and include the noise itself in the calculation if centered) we take the median element that represent the real information
//pixel of the image. 

using namespace std;
using namespace cv;


void padding(Mat, Mat&);
void median_filter(Mat, Mat&);

//Declaring a matrix at runtime: our median filter is composed by all 1 because we need to return the value of the local area itself.
//This kind of declaration assign to a Mat Object a new Mat Object composed by Mat_ (a Mat dataType) of unsigned char of dimension 3x3:
//we proceed to print into this matrix the elements one by one by redirection <<).
Mat median = (Mat_<uchar>(3,3) << 1, 1, 1, 
								  1, 1, 1, 
								  1, 1, 1);
								  
int main(int argc, char** argv) {
	
	//check for a valid image
	if(argc!=2) {
		cerr << "Usage: ./programname imagename.<format>" << endl;
		exit(EXIT_FAILURE);	
	}
	
	//assign to raw_image the content of argv[1] (the image passed in input)
	Mat raw_image = imread(argv[1], IMREAD_GRAYSCALE); //CV_8U
	
	//if is empty exit
	if(raw_image.empty()) {
		cerr << "Image format must be valid." << endl;
		exit(EXIT_FAILURE);
	}
	
	//declare the padded image and the destination image
	Mat padded_image, dest_image;
	
	//now let's call the function to make the zero padding to our padded image
	padding(raw_image, padded_image);
	
	//now let's call the median filter
	median_filter(padded_image, dest_image);
	
	//now let's show the images
	imshow("Original image", raw_image);
	imshow("Filtered image", dest_image);
	waitKey(0);
	return(0);


}


void padding(Mat raw_image, Mat& padded_image) {

	int pad = floor(K/2);
	padded_image = Mat::zeros(raw_image.rows+2*pad, raw_image.cols+2*pad, raw_image.type());
	for(int i=0; i<raw_image.rows; i++) {
		for(int j=0; j<raw_image.cols; j++) {
			padded_image.at<uchar>(i+pad, j+pad) = raw_image.at<uchar>(i,j);
		}
	}
	
}

void median_filter(Mat padded_image, Mat& dest_image) {

	int pad = floor(K/2);
	vector<int> median_vector;
	dest_image = Mat::zeros(padded_image.rows-pad, padded_image.cols-pad, padded_image.type());
	
	for(int i=pad; i<padded_image.rows-pad; i++) {
		for(int j=pad; j<padded_image.cols-pad; j++) {
			for(int u=0; u<median.rows; u++) {
				for(int v=0; v<median.cols; v++) {
					median_vector.push_back((int)padded_image.at<uchar>(i-pad+u, j-pad+v)*median.at<uchar>(u,v));
				}			
			}	
			
			sort(median_vector.begin(), median_vector.end());
			dest_image.at<uchar>(i-pad, j-pad) = median_vector.at(floor(median_vector.size()/2));
			median_vector.clear();	
			
		}
	
	}


} 

/*

	int pad = floor(K/2);
	vector<int> median_vector;
	dest_image = Mat::zeros(padded_image.rows-pad, padded_image.cols-pad, padded_image.type());
	
	for(int i=pad; i<padded_image.rows-pad; i++) {
		for(int j=pad; j<padded_image.cols-pad; j++) {
			for(int u=0; u<median.rows; u++) {
				for(int v=0; v<median.cols; v++) {
					median_vector.push_back((int)padded_image.at<uchar>(i-pad+u, j-pad+v)*median.at<uchar>(u,v));				
				}
			}	

			sort(median_vector.begin(), median_vector.end());
			dest_image.at<uchar>(i-pad, j-pad) = median_vector.at(median_vector.size()/2);
			median_vector.clear();
			
		}
	}

*/
