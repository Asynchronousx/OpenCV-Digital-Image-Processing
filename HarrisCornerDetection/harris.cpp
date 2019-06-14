#include <opencv2/opencv.hpp>
#include <iostream>
#include <math.h>
#include <algorithm>

//Harris Corner Detection
//The Harris Corner Detection Algorithm is useful to take track of all the corner present into an image.
//It take advantage of the changing of the gradient (given by the couple intensity and direction) to know if a pixel region
//contains or not a corner. This is easily understandable by analyzing the gradient in that region: if it contains a STRONG variation
//in ALL the direction, then that region contains a corner.
//If it contains NO VARIATION (or small ones) then it's a flat area, if it contains variation only on the direction of the edge, is an edge:
//(we remember that the edge is ORTHOGONAL to the gradient, to the variation is only perceived around the edge itself.)
//Then, in a corner-region, we can find an HIGH magnitude both on Dx and Dy.
//How we can find all the corner point then?
//By the Harris Formula:
//E(u,v) = (Summation of X,Y)W(x,y)*[I(x+u, y+v) - I(x,y)]^2, where W is a window function that could be either UNITARY (every pixel where a corner is found
//got weight 1) or GAUSSIAN (the corner pixel got the strongest weight, the neighbour less).
//By Taylor expansion, we can say that:
//E(u,v) =~ [u,v]*M*[u].
//					[v]
//Where M is the COVARIANCE Matrix, expressed by:
//M = (summation of X,Y) W(x,y)*[Ix^2, Ixy]
//								[Ixy,  Iy^2]
//Where Ix, and Iy are the DERIVATE of x,y in that specific point in the region, calculated with a differential operator like sobel.
//Now that we have this matrix and the harris formula, how we can CHOOSE if a pixel is edge, flat or a corner? with the HARRIS RESPONSE 
//(or Corner Response), identified by R: it's given by
//R = det(M) + k(trace(M))^2, where
//det(M) is given by Lamba1*Lambda2 (the eigenvalues of M)
//K is a small value (between 0,04 and 0,06)
//and trace(M) is Lambda1+Lambda2.
//Note: we can calculate Lambda 1&2 by the function eigen(M, vector(2)) of the matrix M in each point of the region centered.
//the result R tells us if is a corner, edge or flat: 
//|R| is really small: flat area, Lamba1 & Lambda2 very small BUT still similar L1 =~ L2
//|R| < 0: Lamba 1>>Lambda2 or Lambda2 >> Lambda1: edge area
//|R| is big: Lambda1 & Lambda2 are big and similar, L1 =~ L2.
//This rules are given by the ELLIPTIC APPROSSIMATION of all the gradient vector.
//Thats because L1 & L2 represent the Major and Minor axis of the elliptic circonference.

using namespace std;
using namespace cv;

void printKernel(Mat);
void createGaussian(Mat&, double);
void padding(Mat, Mat&, int);
void blur(Mat, Mat, Mat&);
pair<Mat, Mat> sobel(Mat, Mat&);
void harris(Mat, Mat, Mat&, float);

class EigenValue {
	public:
	int x;
	int y;
	long int lambda2;
	
	EigenValue(int x, int y, long int lambda2) {
		this->x = x;
		this->y = y;
		this->lambda2 = lambda2;	
	}
};

int main(int argc, char** argv) {

	if(argc!=3) {
		cerr << "Usage: ./<programname> <imagename.format> <thresold>" << endl;
		exit(EXIT_FAILURE);
	}
	
	Mat raw_image = imread(argv[1], IMREAD_GRAYSCALE);
	float upper_t = atoi(argv[2]);
	
	if(raw_image.empty()) {
		cerr << "Image format not valid." << endl;
		exit(EXIT_FAILURE);
	}
	
	Mat dest_image = Mat::zeros(raw_image.rows, raw_image.cols, raw_image.type());
	Mat padded_image;
	//Create a non-aggressive Gaussian Kernel with a proper sigma 
	Mat kernel = (Mat_<double>(5,5));
	
	harris(raw_image, kernel, dest_image, upper_t);
	
	
	imshow("Original", raw_image);
	imshow("Harris Image", dest_image);
	waitKey(0);
	

	return 0;

}

void printKernel(Mat mat) {
	
	for(int i=0; i<mat.rows; i++) {
		for(int j=0; j<mat.cols; j++) {
			cout << mat.at<double>(i,j) << " ";
		}
		cout << endl;
	}
	
}

void createGaussian(Mat& kernel, double sigma) {

	double sum = 0.0;
	//Gaussian Kernel is defined by iterating from a left bound to a right bound: 
	//I.e: Kernel 5x5 -> Leftbound = -2, Rightbound = 2, pixel in the middle: (0,0)
	int leftBound = -1*(kernel.rows/2);
	int rightBound = 1*(kernel.rows/2);
	
	//Defining a multiplicative constant to obtaining integer gaussian values: with 1/(exp(-2^2+-2^2)/2*sigma^2)
	//We have found the number that, divided by 1, got the smallest number to mutliply a gaussian values to obtain an integer.
	int mult_const = 1/exp(-(pow(leftBound,2) + pow(leftBound,2)) / (2*pow(sigma,2)));
	
	for(int i=leftBound; i<=rightBound; i++) {
		for(int j=leftBound; j<=rightBound; j++) {
			//G(i,j) = exp(-(pow(x,2) + pow(y,2))/2*pow(sigma,2)) multiplied to mult_const  and rounded to have an integer value
			kernel.at<double>(i+rightBound, j+rightBound) = cvRound(exp(-(pow(i,2) + pow(j,2)) / (2*pow(sigma,2))) * mult_const);
			
			//Cumulate the sum for the normalization
			sum += kernel.at<double>(i+rightBound, j+rightBound);
		}
	}

	//Print info about the kernel
	cout << "Cumulative sum is: " << sum << endl;
	cout << "Sigma is: " << sigma << " - Sigma^2 is: " << pow(sigma,2) << endl << endl;
	cout << "Gaussian kernel without normalization:" << endl;
	printKernel(kernel);
	cout << endl;
	
	//let's now normalize the kernel
	for(int i=0; i<kernel.rows; i++) {
		for(int j=0; j<kernel.cols; j++) {
			kernel.at<double>(i,j) /= sum;
		}
	}

	//print normalized kernel
	cout << "Gaussian Kernel with Normalization:" << endl;
	printKernel(kernel);
	
	
}

//Create a padding to the image equal to kCol or kRow / 2 (floored)
void padding(Mat raw_image, Mat& padded_image, int kSize) {

	int pad = floor(kSize/2);
	padded_image = Mat::zeros(raw_image.rows+2*pad, raw_image.cols+2*pad, raw_image.type());
	
	for(int i=0; i<raw_image.rows; i++) {
		for(int j=0; j<raw_image.cols; j++) {
			padded_image.at<uchar>(i+pad, j+pad) = raw_image.at<uchar>(i,j);
		}
	}

}


//Blur the image using the convolution between the Kernel Matrix and the image itself
void blur(Mat padded_image, Mat kernel, Mat& dest_image) {

	int pad = floor(kernel.rows/2);
	double pixel_value = 0.0;
	
	//convolute G(x,y) & I(x,y)
	for(int i=pad; i<padded_image.rows-pad; i++) {
		for(int j=pad; j<padded_image.cols-pad; j++) {
			for(int u=0; u<kernel.rows; u++) {
				for(int v=0; v<kernel.cols; v++) {
					pixel_value += padded_image.at<uchar>(i-pad+u, j-pad+v) * kernel.at<double>(u,v);		
				}
			}
			
			//assing the cumulate pixel value and reset the counter
			dest_image.at<uchar>(i-pad, j-pad) = floor(pixel_value);
			pixel_value = 0;
			
		}
	}

}

//Differential operator to generate the magnitude intensity and the direction of an image, through
//Their two derivate Dx & Dy
pair<Mat,Mat> sobel(Mat padded_image, Mat& dest_image) {

	Mat sobel_X = (Mat_<char>(3,3) << -1, 0, 1,
									  -2, 0, 2,
									  -1, 0, 1);
									  
	Mat sobel_Y = (Mat_<char>(3,3) << 1,  2,  1,
									  0,  0,  0,
									 -1, -2, -1);
									 
	int pad = floor(sobel_X.rows/2);
	int gradientX = 0;
	int gradientY = 0;
	int pixel_intensity = 0;
	
	//Init DX & DY image: here, dest clone is a Mat filled with zeros with our desired size and type,
	//So just let clone it
	Mat DX = dest_image.clone();
	Mat DY = dest_image.clone();
	
	
	//convolute sobel matrixes on the image
	for(int i=pad; i<padded_image.rows-pad; i++) {
		for(int j=pad; j<padded_image.cols-pad; j++) {
			for(int u=0; u<3; u++) {
				for(int v=0; v<3; v++) {
					//calculate the derivates
					gradientX += padded_image.at<uchar>(i-pad+u, j-pad+v) * sobel_X.at<char>(u,v);
					gradientY += padded_image.at<uchar>(i-pad+u, j-pad+v) * sobel_Y.at<char>(u,v);
				}
			}
				
						
			//calculate the abs value
			
			gradientX = abs(gradientX);
			gradientY = abs(gradientY);
			
			//calculate the pixel intensity and normalize
			//NB: we don't need to normalize the derivate because we'll need the real value into the next step
			//to calculate the matrix M of harris.
			pixel_intensity = abs(gradientX) + abs(gradientY);
			pixel_intensity = pixel_intensity > 255? 255: pixel_intensity < 0? 0: pixel_intensity;
			
			//insert the derivate into the images
			dest_image.at<uchar>(i-pad, j-pad) = pixel_intensity;
			DX.at<uchar>(i-pad, j-pad) = gradientX;
			DY.at<uchar>(i-pad, j-pad) = gradientY;
			
			//reset the counter
			gradientX = 0;
			gradientY = 0;
			pixel_intensity = 0;
			
		}
	}
	

									 
	pair<Mat, Mat> derivates = make_pair(DX, DY);
	return derivates;	
	 
}

//The Harris Algorithm is the one that will help us find the corner of the image
void harris(Mat raw_image, Mat kernel, Mat& dest_image, float upper_t) {

	//Harris is composed by previous step. 
	//1) The first step is to calculate the derivates (X,Y) of an image, and we can do that
	//by using sobel. Also we need a padded image to works on. 
	Mat dest_sobel_image = dest_image.clone();
	Mat padded_sobel_image;
	
	//let's pad the image to be filtered by sobel matrixes (3x3)
	padding(raw_image, padded_sobel_image, 3);
	
	//lets use sobel and retrieve the pair of derivate images
	pair<Mat, Mat> derivates = sobel(padded_sobel_image, dest_sobel_image);

	//retrieve the image
	Mat DX = derivates.first;
	Mat DY = derivates.second;
	
	//show the images
	imshow("DX", DX);
	imshow("DY", DY);
	imshow("Sobel", dest_sobel_image);
	waitKey(0);
	
	//2) The second step is smooth the images, DX and DY to obtain a better image to work on.
	//	 We'll use a gaussian filter with 1.2 sigma to be non-aggressive.
	double sigma = 1.2;
	createGaussian(kernel, sigma);
	
	//let's smooth both image: create their padding and pass them to the blur filter
	Mat DX_P, DY_P;
	padding(DX, DX_P, kernel.rows);
	padding(DY, DY_P, kernel.rows);
	blur(DX_P, kernel, DX);
	blur(DY_P, kernel, DY);
	
	//show the blurred images
	imshow("Blurred DX", DX);
	imshow("Blurred DY", DY);
	waitKey(0);	
	
	//3) Third step is to calculate the covariance matrix M. we can do this by simply iterating 
	//	 through the image with a window function W of size NxN. in each window region, lets calculate
	//   the derivate X^2, Y^2 and XY for the other two diagonal. Then let's find the R (harris response)
	//   and choose if the current point is an iphotetical corner or not. If yes, add it to the list L of point.
	int wSize = 3;
	int pad = floor(wSize/2);
	float det, trace, R;
	float k = 0.04;
	vector<EigenValue> corners;
	Mat padded_image;
	padding(raw_image, padded_image, wSize);
	
	//Build the covariance matrix: for the size of the image
	for(int i=pad; i<padded_image.rows-pad; i++) {
		for(int j=pad; j<padded_image.cols-pad; j++) {
		
			//covariance matrix
			Mat M = Mat::zeros(2, 2, CV_32FC1);
			//eigen value matrix
			Mat eigenVal = Mat::zeros(1, 2, CV_32FC1);
		
			//Windows function w(u,v) to calculate the covariance matrix in a neighbour wSizeXwSize (3x3 in this case)
			for(int u=0; u<wSize; u++) {
				for(int v=0; v<wSize; v++) {
					//First element is derivate X squared 
					M.at<float>(0,0) += pow(DX.at<uchar>(i-pad+u, j-pad+v), 2);
					//Second element is the product of DX x DY
					M.at<float>(0,1) += DX.at<uchar>(i-pad+u, j-pad+v) * DY.at<uchar>(i-pad+u, j-pad+v);
					//Third element is again the product of DX x DY
					M.at<float>(1,1) += DX.at<uchar>(i-pad+u, j-pad+v) * DY.at<uchar>(i-pad+u, j-pad+v);
					//Last element is Derivate Y squared
					M.at<float>(1,0) += pow(DX.at<uchar>(i-pad+u, j-pad+v), 2);
				}
			}
			
			
			//calculate the eigen values
			eigen(M, eigenVal);
			
			//calculate the harris response: 
			//We can calculate the harris response R, then confronting R to the upper thresold
			//Or just confront the minimum of the lambda generated from the eigen calculus to the thresold.
			//Note: In case we should choose to use R for the comparison, the upper thresold should be A LOT
			//bigger: for example, with a chessboard image, the upper thresold in the case of R should be 28Milion,
			//while for the comparison with the minimum lambda, only 2000.
			
			/* Useful only if considering R as the comparison values; we can just pick Lambda2.
			det = eigenVal.at<float>(0,0) * eigenVal.at<float>(0,1);
			trace = eigenVal.at<float>(0,0) + eigenVal.at<float>(0,1);
			R = det + k*(pow(trace, 2)); //also: could be det + k(pow(trace,2))
			*/
			
			//check if R or Lambda2 is big enough to be a corner, or check directly Lambda2 (the minimum of the two eigenvalues)
			float minLambda = min(eigenVal.at<float>(0,0), eigenVal.at<float>(0,1));
			if(minLambda > upper_t) {
				//if so, build the point and pushback into the vector
				EigenValue point(i-pad, j-pad, minLambda);
				corners.push_back(point);
			}
			
			//clear values
			det = 0;
			trace = 0;
			R = 0;
			
		}
	}
	
	
	//4) Non-corner suppression. 
	//Now we need to sort the corner vector and suppress the non-corner into the vector.
	//We can suppress the non-corner iterating a windows of size NxN to the image, and check
	//if analyzing the current pixel we found another corner point in his neighbour, we must 
	//erase that corner from the corner vector.
	
	//sort in descending order from the bigger lambda2
	sort(corners.begin(), corners.end(), [](const EigenValue& l, const EigenValue& r) { return l.lambda2 > r.lambda2; });
	
	//Non-Corner suppression: Here the logic is very simple. 
	//We've sorted the corner vector by lambda2, so now we know that we got from the most important corner (that surely
	//is a corner) to the least important one. From the theory, to suppress the non-maximum corner, is a corner B
	//is in the neighbourhood of a current analyzed corner A, we know for sure that, since we're analyzing A first
	//A is more important than B. so we need to suppress B if is in the range of the neighbourhood of A.
	
	//TO do this, just iterate through the list of corner (0 -> N)
	for(int i=0; i<corners.size(); i++) {
		//Now's the trick: j is i+1: thats because at the i-th iteration, we know that J could not be LESS
		//than than the current element PLUS one (i+1), because, since the corners vector is sorted, we've
		//suppressed the non corner in the last i+1 iteration. So we need only to check if there could be a 
		//non-corner in the remaining j=i+1 elements.
		for(int j=i+1; j<corners.size(); j++) {
			//Now, we say that if the point at j we're analyzing is included in x+pad and x-pad, as well
			//as y+pad and y-pad, then we must erase this point because it is into the neighbourhood. 
			//With the AND condition we're not only analyzing VERTICAL AND HORIZONTAL bound, but also
			//the diagonal, negative and positive
			if(corners.at(j).x <= corners.at(i).x+pad && corners.at(j).x >= corners.at(i).x-pad &&
			   corners.at(j).y <= corners.at(i).y+pad && corners.at(j).y >= corners.at(i).y-pad) {
			   		//if the condition is met, erase the pixel
			   		corners.erase(corners.begin() + j);
					//decrement the j index: that's because we're still into the J for and we need to inform
			   		//the cycle we decremented the index by erasing.
			   		j--;
			   }
		}			
	}
	
	
	//Copy the dest image and draw the circles
	dest_image = raw_image.clone();
	
	//Colorful circles: convert to BGR
	cvtColor(dest_image, dest_image, COLOR_GRAY2BGR);
	
	for(int i=0; i<corners.size(); i++) {
		//Note: we must consider the Y first because Point(x,y) consider x as column and y as row.
		Point center(cvRound(corners.at(i).y), cvRound(corners.at(i).x));
		//Create a circle in the dest image, of coordinate center, radius 5, scalar for identify the color, 
		//1 for the line thickness, and 4 for the line type (last is 0, represent the shift of bit)
		circle(dest_image, center, 5, Scalar(0,255,0), 1, 4, 0);
	}
		
}

