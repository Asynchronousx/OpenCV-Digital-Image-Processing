#include <opencv2/opencv.hpp>
#include <math.h>
#include <iostream>

using namespace std;
using namespace cv;

//Distance transform:
//The distance transforms is a tool used to many things: above all those things, converting
//a binary image to a grayscale one following particoular rules.
//They grayscale is assigned in a particoular way such as the distance from that very pixel from the background.
//The more the pixel is near the background, the less intensity he has.
//So, a pixel can vary it's color depending on his position and the position of the background.
//The lesser he's near the background, the whiter his intensity is. 
//The more he's near the background, the blacker his intensity is. 
//All the pixel in the middle assume a grayscale value.
//Obviously, in the center of the image, all the pixel will be whiter and brighter than every part in the image.
//The distance transform represent the distance from a pixel from the background.
//Another use the DT could have is to compress image: from a local area, taken the maximum value (alongside the (x,y)
//of the pixel with the max value) and the structuring element, we can compress the image taking a list of those pixel
//and saving the mask. To build the original image again, we should convolve the structuring element with every pixel
//saved into the list into a new image of size MxN (where M,N were the original size of the image). 

void DistanceTransform4(Mat, Mat&);
void DistanceTransform8(Mat, Mat&);

int main(int argc, char** argv) {

	if(argc!=2) {
		cerr << "Usage: ./<programname> <imagename.format>" << endl;
		exit(EXIT_FAILURE);	
	}

	Mat raw_image = imread(argv[1], IMREAD_GRAYSCALE);
	
	if(raw_image.empty()) {
		cerr << "Image format not valid." << endl;
		exit(EXIT_FAILURE);
	}
		
	//Declaring images for both the 4-8 connectivity DT function.
	Mat binary_image4, binary_image8;
	Mat dest_image4, dest_image8;
	
	//Mat gradient_image; only if extracting the shape and filling.
	
	//Before passing the image into the DT function, we should take care
	//of the pre-processing of the image.
	//First pre-processing step: get the threshold image value: we're using
	//OTSU threshold to have an adaptive and realistic threshold.
	//NOTE: We must choose a standard to analyze our image.
	//A binary image can be: BLACK in the background, and WHITE into the object, 
	//OR
	//White into the background and BLACK into the object. 
	//The standard used in this algorithm is BLACK for the background and WHITE for
	//the actual object.
	threshold(raw_image, binary_image4, 127, 255, THRESH_BINARY /*| THRESH_OTSU */ );

	binary_image8 = binary_image4.clone();
	
	imshow("Binary Image", binary_image4);
	waitKey(0);
	
	
	//---OPTIONAL FOR A BETTER RESULT---
	//Now we should take care of some other few pre-processing steps:
	//image segmentation, border extracting and image filling.
	/*
	
	//Before extracting the border, an image segmentation should be ideal: because,
	//extracted the shape of an object, the result of the DT should be even better
	
	//TODO: Implement image segmentation
	
	//Extract border
	morphologyEx(binary_image, gradient_image, MORPH_GRADIENT, Mat());
	imshow("Edge Extraction", gradient_image);
	waitKey(0);
	
	//Filling the image at this point, means that we're filling the shape of the image.
	//With this procedure, all the holes or noise into the image are removed, and the result of the DT
	//will be even better.
	
	TODO: Implement image filling
	
	*/
	
	DistanceTransform4(binary_image4, dest_image4);
	DistanceTransform8(binary_image8, dest_image8);
	
	imshow("DT 4-Conn", dest_image4);
	imshow("DT 8-Conn", dest_image8);
	imwrite("4.png", dest_image4);
	imwrite("8.png", dest_image8);
	waitKey(0);
	return 0;

}


void DistanceTransform4(Mat binary_image, Mat& dest_image) {

	//The logic behind the DT is simple: we need two scanning: 
	//The first, is from up to bottom, left to right: we're analyzing how much a pixel
	//is far from the background. Using a 4-Connected operator, we must analyze, for every pixel,
	//the upper and left pixel into his neighbourhood.
	//In the second scanning, from bottom to up, right to left, we need to check, actually, if the 
	//distance of the scanned pixel is correct: we do that by taking the right and bottom pixel (because
	//of the 4-connected operator and because all the pixels are already calculated.	
	for(int i=1; i<binary_image.rows; i++) {
		for(int j=1; j<binary_image.cols; j++) {
			//If the current pixel isn't a background pixel
			if(binary_image.at<uchar>(i,j) != 0) {
				//Take the left and upper pixel
				int N = binary_image.at<uchar>(i-1, j);
				int W = binary_image.at<uchar>(i, j-1);

				//having the already analyzed neighbourhood, take the min between them and add
				//1 to set the increased distance.
				binary_image.at<uchar>(i,j) = 1 + min(N,W);
			}
		}
	}
	
	//Doing the second scan, bottom to up, right to left. 
	for(int i=binary_image.rows-1; i>=0; i--) {
		for(int j=binary_image.cols-1; j>=0; j--) {
			//as before, considering the current pixel only if not background:
			//Why? because we're analyzing the distance of an OBJECT from the BACKGROUND.
			//so analyzing background pixel will produce a FULL distance image, and not the shape we want.
			if(binary_image.at<uchar>(i,j) != 0) {
				//taking the right and bottom pixel and increase them by one.
				//In the second scanning alle the pixel are defined. So why take only 2 instead of 4?
				//Because we're analyzing the distance from the background from right to left, bottom to up.
				//Taking all the pixel in all the direction would be useless and produce the same exact result.
				int S = binary_image.at<uchar>(i+1, j);
				int E = binary_image.at<uchar>(i, j+1);
			
				//Increment those pixel by one
				S++;
				E++;
		
				//taking the min between those pixel and the current pixel, and putting that
				//into the binary image itself: that's because in the next iteration, we will
				//analyze the corrected pixel at an i-th step.
				binary_image.at<uchar>(i,j) = min((int)binary_image.at<uchar>(i,j), min(S,E));
			}
		}
	}

	//normalize the image
	normalize(binary_image, dest_image, 0, 255, NORM_MINMAX);

}

//TODO: Implement 8Connectivity
void DistanceTransform8(Mat binary_image, Mat& dest_image) {

	//As for the 4-Connectivity, the logic is the same: we must do two scanning,
	//left to right up to bottom, and bottom to up, right to left.
	//In the 8-Connectivity we must consider the top-left pixels when iterating from left to right up to bottom.
	//We start at 1 to avoid padding
	for(int i=1; i<binary_image.rows; i++) {
		for(int j=1; j<binary_image.cols; j++) {
			//If the current pixel isn't a background pixel
			if(binary_image.at<uchar>(i,j) != 0) {
				//In the 8-Conn we need to extract 4 pixels: left, upper left, upper and upper right.
				int W = binary_image.at<uchar>(i, j-1);
				int NW = binary_image.at<uchar>(i-1, j-1);
				int N = binary_image.at<uchar>(i-1, j);
				int NE = binary_image.at<uchar>(i-1, j+1);
				
				//Now we need to take the min between those four: 
				int m = min(min(W,NW), min(N,NE));
				
				//put this pixel into the current index and increment by one, to set and
				//memorize the incremented distance.
				binary_image.at<uchar>(i,j) = m + 1;
			}
		}
	}
	
	//Second scan: bottom to up, right to left
	for(int i=binary_image.rows-1; i>=0; i--) {
		for(int j=binary_image.cols-1; j>=0; j--) {
				//Again, if the current pixel isn't a background pixel
				if(binary_image.at<uchar>(i,j) != 0) {
				//Again, we need to take 4 pixel: bottom left, bottom, bottom right, right
				//for the same reason of the second scanning of the 4-conn.
				//We also must increment them by one.
				int SE = binary_image.at<uchar>(i+1, j+1);
				int S = binary_image.at<uchar>(i+1, j);
				int SW = binary_image.at<uchar>(i+1, j-1);
				int E = binary_image.at<uchar>(i,j+1);
				
				//increment them by one to remember the distance
				SE++;
				S++;
				SW++;
				E++;
				
				//take the min between them 
				int m = min(min(SE,S), min(SW,E));
				
				//confront the min with the current pixel and choose the min
				binary_image.at<uchar>(i,j) = min((int)binary_image.at<uchar>(i,j), m);
			}
		}
	}

	//normalize the matrix
	normalize(binary_image, dest_image, 0, 255, NORM_MINMAX);

}
