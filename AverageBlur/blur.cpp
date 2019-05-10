#include <opencv2/opencv.hpp>
#include <iostream>
#include <math.h>
#define K 3

using namespace std;
using namespace cv;

//Create the kernel: we're creating a convolution/filtering matrix populated with 1 to blur our image.
//For that, we use the functionality of the native dataType Mat_ with the aux of the input stream to print into 
//our matrix the numbers.
Mat average_filter = (Mat_<uchar>(3,3) << 1, 1, 1,
										  1, 1, 1,
										  1, 1, 1);
										  
int calculate_filter_value(Mat);										
void create_padding(Mat, Mat&);
void blur(Mat, Mat&);

int main(int argc, char** argv) {
	
	//check for error
	if(argc!=2) {
		cerr << "Usage: ./<programname> <imagename>.<format>" << endl;
		exit(EXIT_FAILURE);
	}
	
	//declare the images
	Mat raw_image = imread(argv[1], IMREAD_GRAYSCALE);	
	Mat padded_image, dest_image;
	
	//check for an empty data
	if(raw_image.empty()) {
		cerr << "Image format not valid. Try again." << endl;
		exit(EXIT_FAILURE);
	}
	
	//pad the image with 0 padding
	create_padding(raw_image, padded_image);
	
	//blur the image
	blur(padded_image, dest_image);
	
	//show the images
	imshow("Original image", raw_image);
	imshow("Blurred Image", dest_image);
	waitKey(0);
	return(EXIT_SUCCESS);
	
}


//How we can effectively calculate an average filter on an image?
//From the theory, the g(x,y) is given by:
//(Summatory of u=0 to rows, Summatory of v=0 to cols)*(w(u,v)*f(x+u, y+v)) / (Summatory of u=0 to rows, Summatory of v=0 to cols) w(u,v).
//So, we need to calculate the sum of all the multiplication between the filter and the image, and divide this number by the totale value
//of the elements into the filter matrix. To do this, we iterate through the matrix summing up all the numbers and returning the total value.
int calculate_filter_value(Mat filter) {
	
	int value = 0;
	
	for(int i=0; i<filter.rows; i++) {
		for(int j=0; j<filter.cols; j++) {
			value += (int)filter.at<uchar>(i,j);
		}
	}
	
	return value;

}


//in this function we need to create a padded image (done with 0-padding) to pass our filter on. 
//we're doing that because otherwise we would miss some pixel passing the filter on it, because it's the central 
//pixel of the kernel that does the operation on the current pixel on the image.
//So, before we can center the central pixel of an image, we should do a left-top-bottom-right padding to let the convolution
//matrix center every pixel in our f(x,y). Without padding, we would miss many rows as half the dimension of the kernel rows.
void create_padding(Mat raw_image, Mat& padded_image) {

	
	//padding value: equal to Floor(Dim Kernel rows / 2)
	int pad = floor(K/2);
	
	//now let's create a new image with the padded cols/rows filled with 0. 
	//Note that: we need to add the padding in the left and upper part, but also into the rightmost and bottom one.
	//so we need to multiply 2 * the padding value.
	padded_image = Mat::zeros(raw_image.rows+2*pad, raw_image.cols+2*pad, raw_image.type());
	
	
	//now let's copy the raw image into the padded one, and we need to do that starting from the real raw image size:
	//that's because otherwise we would fill the padded image starting from [0][0], nullifing the padding (and lefting 
	//a lot of 0 into the mat element).
	//We could both start from the raw image size (rows and col) or the padded image size MINUS the pad value (padded_image.rows-2*pad) etc.
	for(int i=pad; i<raw_image.rows; i++) {
		for(int j=pad; j<raw_image.cols; j++) {
			//now, assign the current image pixel to the the empty padded one starting from [0+pad][0+pad]
			padded_image.at<uchar>(i+pad, j+pad) = raw_image.at<uchar>(i,j);
		}
	} 

}




//This is the real function that blur our image. What is blurring? BLurring is given by the Average filter.
//We can obtain blur by convolute a filtering Kernel (g(x,y)) on an image (f(x,y)), resulting in a new image with the value of its pixel
//mediated by some mathematical operation. The blurring resemble the Low Pass filter, because we're attenuating the high frequencies.
//We can do this by convolute the Kernel on the Image, and processing for every pixel of the image the mathematical operation in which 
//we assign to a value the sum of all the "mediated" value of the pixel of the image multiplied by the pixel of the kernel, so we need to 
//multiply every kernel pixel (and it's adjacent) to the relative pixel (and it's adjacent) on the image. We're doing this by:
//v += i(u,v)*k(u,v), and from this result we need to find the average by dividing it by the total value of the elements of the 
//kernel: r(i,j) = v/(summatory of all the elements in the kernel). 
void blur(Mat padded_image, Mat& dest_image) {

	//useful variables
	//pad is the padding space of the image 
	int pad = floor(K/2);
	//pixel intensity is used to measure the sum of all the local pixel of the image multiplied by the ones of the window filter.
	int pixel_intensity = 0;
	//average value calculated by the summation of all the value of the filter
	int average_value = calculate_filter_value(average_filter);
	//init the dest image
	dest_image = Mat::zeros(padded_image.rows-2*pad, padded_image.cols-2*pad, padded_image.type());
	
	//for the size of the image (not including the padding, thats why i/j starts from 0 and finish in rows/cols - 2*pad)
	for(int i=pad; i<padded_image.rows-2*pad; i++) {
		for(int j=pad; j<padded_image.cols-2*pad; j++) {
			//for the size of the kernel (row/cols)
			for(int u=0; u<average_filter.rows; u++) {
				for(int v=0; v<average_filter.cols; v++) {
					//sum up into pixel intensity the current product between the current pixel and the pixel of the filter
					pixel_intensity += (int)padded_image.at<uchar>(i-pad+u, j-pad+v)*average_filter.at<uchar>(u,v);
				}				
			}
			
			//calculate the average by dividing the resulting pixel intensity for the total value of the matrix
			dest_image.at<uchar>(i-pad, j-pad) = floor(pixel_intensity/average_value);
			//reset the intensity
			pixel_intensity = 0;
			
		}
	} 



}
