#include <opencv2/opencv.hpp>
#include <iostream>
#include <math.h>

//Gaussian Filters are particoular filter builded with a Sigma value (and a Mu value if not 0) that generate a more/less aggressive
//filter based on how big is the Sigma. Gaussian Filters are useful to remove Gaussian Noise and also blur the image.
//From the theory, we know that a Gaussian Filter is definied like this:
//G(x) = (1/Sigma*sqrt(2*pi))*exp(-(pow(x-Mu, 2))/2*pow(sigma,2)) and if Mu = 0 (the gaussian is centered in his average) then:
//G(x) = (1/Sigam*sqrt(2*pi))*exp(-(pow(x, 2))/2*pow(sigma,2)) with the Mu value missing.
//Now, if we would like to operate into a 2D enviroinment (without the Mu value -> Mu = 0): 
//G(x,y) = (1/(sigmaX*SigmaY)*sqrt(2*pi))*exp(-(pow(x,2)/2*pow(sigmaX,2) + (pow(y,2))/2*pow(sigmaY,2)))
//That said, we can assume that, wanting a CIRCULAR GAUSSIAN FILTER, we're coinsidering that SigmaX = SigmaY.
//With that rule, we can build the equivalent of the Gaussian Filtering into the Spatial Domain:
//G(i,j) = (1/(sigma*sqrt(2*pi))*exp(-(pow(x,2) + pow(y,2))/2*pow(sigma,2)))
//We can write a more convenient formula as: G(i,j) = exp(-(pow(x,2) + pow(y,2))/2*pow(sigma,2)) / (2*M_PI*pow(sigma,2))
//The G(i,j) could be simply written by some for cycles, cumulation of the value of the matrix generated and the normalization through
//the division of each value for this cumulation.
//We can build the Matrix by considering the center, the spot where the pixel will pass to be filtered by the Mean Filter technique.
//We're considering then the Central pixel as 0, and the floor(Size/2)*-1 as the left bound and floor(Size/2)*1 as the right bound,

using namespace std;
using namespace cv;

void printMat(Mat);
void createGaussian(Mat&, double);
void createPadding(Mat, Mat&, int);
void blur(Mat, Mat, Mat&);

int main(int argc, char** argv) {

	if(argc!=2) {
		cerr << "Usage: ./filename <imagename>.<format>" << endl;
		exit(EXIT_FAILURE);
	}

	Mat raw_image = imread(argv[1], IMREAD_GRAYSCALE);
	
	if(raw_image.empty()) {
		cerr << "Image format not valid" << endl;
		exit(EXIT_FAILURE);	
	}
	
	//defining the Mat object
	Mat padded_image, dest_image;
	
	//Generating a kernel 5x5 of double and defining the sigma for creating the gaussian filter
	Mat kernel = (Mat_<double>(5,5));
	//if sigma = 1.73, then sigma^2 = 3
	double sigma = 1.20;
	
	//create gaussian, add padding and blur
	createGaussian(kernel, sigma);
	createPadding(raw_image, padded_image, kernel.rows);
	blur(padded_image, kernel, dest_image);
	
	//show
	imshow("Original Image", raw_image);
	imshow("Dest Image", dest_image);
	imwrite("gauss.png", dest_image);
	waitKey(0);
	return(0);
	
}


void printMat(Mat mat) {
	
	for(int i=0; i<mat.rows; i++) {
		for(int j=0; j<mat.cols; j++) {
			cout << mat.at<double>(i,j) << " ";
		}
		cout << endl;
	}
	
}


void createGaussian(Mat& kernel, double sigma) {

	//Defining the cumulatime sum variable
	double sum = 0.0;
	
	//Defining the boundaries of the interval: boundaries are given by the
	//floor of the division of the size (rows or cols) of the kernel. 
	//The left boundary is -1*left, the right is 1*right. 
	//For example, the bound of a matrix 5x5 are: floor(5/2)*-1 -> -2 and floor(5/2)*1 -> 2.
	int leftBound = -1*floor(kernel.rows/2);
	int rightBound = 1*floor(kernel.rows/2);
	
	//to get integer values of our gaussian filter, we need to divide the smaller number by 1: with that,
	//we can make all the value at least ONE, giving us the possibility to make our gaussian kernel composed by integer.
	//The gaussian kernel is simmetric: this said, we can assume that the smallest values are into the angles of the kernel.
	//Calculate it and divide by one to have our costant.
	//Example: if the value is 0.26, to find the multiplicative costant we need to divide this number by 1: then we have 1/0.26 = 3,8461, 
	//rounded to 3,85.
	double mult_constant = 1/exp(-(pow(leftBound,2) + pow(leftBound,2))/(2*pow(sigma,2)));
	
	//Now we need to generate the gaussian filter. We can do this iterating through the size of the boundaries,
	//assigning values to the current kernel element through the gaussian formula.
	for(int i=leftBound; i<=rightBound; i++) {
		for(int j=leftBound; j<=rightBound; j++) {
			//Now let's applyt the fomula: G(i,j) = (1/(sigma*sqrt(2*pi))*exp(-(pow(x,2) + pow(y,2))/2*pow(sigma,2))).
			//We can write a more convenient formula as: G(i,j) = exp(-(pow(x,2) + pow(y,2))/2*pow(sigma,2)) / (2*M_PI*pow(sigma,2))
			//We can also get rid of / (2*M_PI*pow(sigma,2)) because after the cumulative sum we need to normalize the kernel to 
			//integer values.
			//To make the gaussian kernel composed by integer value, we can directly multiply the resultant number from the 
			//exp(-(pow(x,2) + pow(y,2))/2*pow(sigma,2)), rounding to the next integer.
			//We need to dispatch the kernel to reach the right index: just sum rightBound to i/j.
			kernel.at<double>(i+rightBound, j+rightBound) = cvRound(exp(-(pow(i,2) + pow(j,2)) / (2*pow(sigma,2))) * mult_constant);
			
			//now let's cumulate the sum: here too, we need to dispatch the kernel index summing right bound to it.
			sum += kernel.at<double>(i+rightBound, j+rightBound);
		}
	}


	cout << "Cumulative sum is: " << sum << endl;
	cout << "Sigma is: " << sigma << " - Sigma^2 is: " << pow(sigma,2) << endl << endl;
	cout << "Gaussian kernel without normalization:" << endl;
	printMat(kernel);
	cout << endl;
	
	//let's now normalize the matrix 
	for(int i=0; i<kernel.rows; i++) {
		for(int j=0; j<kernel.cols; j++) {
			kernel.at<double>(i,j) /= sum;
		}
	}

	cout << "Gaussian Kernel with Normalization:" << endl;
	printMat(kernel);

}

void createPadding(Mat raw_image, Mat& padded_image, int kSize) {

	//get the padding and create the 0s matrix
	int pad = floor(kSize/2);
	padded_image = Mat::zeros(raw_image.rows+2*pad, raw_image.cols+2*pad, raw_image.type());
	
	//pad the image
	for(int i=0; i<raw_image.rows; i++) {
		for(int j=0; j<raw_image.cols; j++) {
			padded_image.at<uchar>(i+pad, j+pad) = raw_image.at<uchar>(i,j);		
		}
	}
	
}

void blur(Mat padded_image, Mat kernel, Mat& dest_image) {

	//create the padding, defining a pixel intensity variable to hold the sum of moltiplication,
	//and fill the dest image with 0s
	int pad = floor(kernel.rows/2);
	double pixel_intensity = 0.0;
	dest_image = Mat::zeros(padded_image.rows-2*pad, padded_image.cols-2*pad, padded_image.type());
	
	//do the convolution for the size of the image
	for(int i=pad; i<padded_image.rows-2*pad; i++) {
		for(int j=pad; j<padded_image.cols-2*pad; j++) {
			//for the size of the kernel
			for(int u=0; u<kernel.rows; u++) {
				for(int v=0; v<kernel.cols; v++) {
					//Multiply the current image element dispatched by the index of the kernel minus the padding for
					//the current kernel element, and add it to the pixel intensity value as a sum
					pixel_intensity += padded_image.at<uchar>(i-pad+u, j-pad+v) * kernel.at<double>(u,v);
				}	
			}
			
			//assign the pixel once exited the cycle and reset the intensity
			dest_image.at<uchar>(i-pad, j-pad) = floor(pixel_intensity);
			pixel_intensity = 0;		
		
		}
	}
}
