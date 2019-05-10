#include <opencv2/opencv.hpp>
#include <math.h>
#include <iostream>
#include <climits>

using namespace std;
using namespace cv;

void PrintMat(Mat);
void Gaussian(Mat&, double);
void Sobel(Mat, Mat&, Mat&, int&);
void Padding(Mat, Mat&, int);
void Blur(Mat, Mat, Mat&);
Mat Canny(Mat, Mat&, Mat&, Mat&);

int main(int argc, char** argv) {

	if(argc!=2) {
		cerr << "Usage: ./<programname> <imagename>.<format>" << endl;
		exit(EXIT_FAILURE);	
	}
	
	Mat raw_image = imread(argv[1], IMREAD_GRAYSCALE); 
	
	if(raw_image.empty()) {
		cerr << "Image format not valid or Image empty." << endl;
		exit(EXIT_FAILURE);
	}

	Mat padded_image, dest_image;
	Mat kernel = (Mat_<double>(5,5));
	
	Mat canny_image = Canny(raw_image, padded_image, kernel, dest_image);
	
	imshow("Original Image", raw_image);
	imshow("Canny Image", canny_image);
	imwrite("canny.png", canny_image);
	waitKey(0);
	return 0;

}

void PrintMat(Mat mat) {

	for(int i=0; i<mat.rows; i++) {
		for(int j=0; j<mat.cols; j++) {
			cout << mat.at<double>(i,j) <<  " ";
		}
		cout << endl;
	}

}

void Gaussian(Mat& kernel, double sigma) {

	double sum = 0.0;
	double pixel_value = 0.0;
	double leftBound = floor(kernel.rows/2)*-1;
	double rightBound = floor(kernel.rows/2)*1;
	double mult_constant = 1/exp(-(pow(leftBound, 2) + pow(leftBound, 2)) / (2*pow(sigma,2)));
	
	//for the size of the gaussian kernel, build it
	for(int i=leftBound; i<=rightBound; i++) {
		for(int j=leftBound; j<=rightBound; j++) {
			//gaussian formula
			pixel_value = exp(-(pow(i, 2) + pow(j, 2)) / (2*pow(sigma,2)));
			
			//multiplying the value for the costant to make it integer and rounding it
			kernel.at<double>(i+rightBound, j+rightBound) = cvRound(pixel_value * mult_constant);
			
			//cumulate the sum
			sum += kernel.at<double>(i+rightBound, j+rightBound);
		}
	}
	
	cout << "Gaussian Kernel without Normalization:" << endl;
	PrintMat(kernel);
	
	for(int i=0; i<kernel.rows; i++) {
		for(int j=0; j<kernel.cols; j++) {
			kernel.at<double>(i,j) /= sum;
		}
	}
	
	cout << "Gaussian Kernel with Normalization:" << endl;
	PrintMat(kernel);

}


//Maybe its sobel the problem, maybe its the non-max supression.
void Sobel(Mat padded_image, Mat& dest_image, Mat& gradient_orientations, int& max_intensity) {

	//Declaring the Sobel X and Y matrix
	Mat sobel_X = (Mat_<char>(3,3) << -1, 0, 1,
									  -2, 0, 2,
									  -1, 0, 1);
							
	Mat sobel_Y = (Mat_<char>(3,3) << 1,  2,  1,
									  0,  0,  0,
									 -1, -2, -1);
									 
	//Declaring useful variables
	int gradientX = 0;
	int gradientY = 0;
	int pad = floor(3/2);
	int pixel_intensity = 0.0;
	double theta = 0.0;
	double degree = 0.0;
	
	//Filling the empty matrix
	dest_image = Mat::zeros(padded_image.rows-2*pad, padded_image.cols-2*pad, padded_image.type());
	gradient_orientations = Mat::zeros(dest_image.rows, dest_image.cols, dest_image.type());
	
	/* Matrix of the DX & DY images if needed 
	Mat DX = dest_image.clone();
	Mat DY = dest_image.clone();
	*/
	
	//Calculate the Sobel image
	for(int i=pad; i<padded_image.rows-pad; i++) {
		for(int j=pad; j<padded_image.cols-pad; j++) {
			for(int u=0; u<3; u++) {
				for(int v=0; v<3; v++) {
					//build the x derivate
					gradientX += padded_image.at<uchar>(i-pad+u, j-pad+v) * sobel_X.at<char>(u,v);
					//build the y derivate
					gradientY += padded_image.at<uchar>(i-pad+u, j-pad+v) * sobel_Y.at<char>(u,v);
				}
			}
			
			/* Without convolution
			gradientX = abs(((-1*padded_image.at<uchar>(i-1, j-1)) + (-2*padded_image.at<uchar>(i, j-1)) + (-1*padded_image.at<uchar>(i+1, j-1))) +
						(padded_image.at<uchar>(i-1, j+1) + 2*padded_image.at<uchar>(i, j+1) + padded_image.at<uchar>(i+1,j+1)));
			//build the y derivate
			gradientY = abs((padded_image.at<uchar>(i-1, j-1) + 2*padded_image.at<uchar>(i-1, j) + padded_image.at<uchar>(i-1, j+1)) +
						((-1*padded_image.at<uchar>(i+1, j-1)) + (-2*padded_image.at<uchar>(i+1, j)) + (-1*padded_image.at<uchar>(i+1, j+1))));
			*/
			
			//calculate the absolute value
			gradientX = abs(gradientX);
			gradientY = abs(gradientY);
			
			//calculate intensity and orientation by the sobel formula
			pixel_intensity = gradientX + gradientY;
			theta = atan2(gradientY, gradientX);
			degree = theta*180/M_PI;
			
			//checking the intensity 
			if(pixel_intensity > max_intensity) {
				max_intensity = pixel_intensity;
			}
			
			//normalize the pixel
			pixel_intensity = pixel_intensity > 255? 255: pixel_intensity < 0? 0: pixel_intensity;
			
			/* pixel normalizing for the gradient if needed to output the DX & DY images
			gradientX = gradientX > 255? 255: gradientX < 0? 0: gradientX;
			gradientY = gradientY > 255? 255: gradientY < 0? 0: gradientY;
			*/
			
			//normalize the degree
			if(degree > 180) {
				while(degree > 180) degree -= 180;
			} 
			else if (degree < 0) {
				while(degree < 0) degree += 180;
			}			
			
			//discover the angle: here we got 4 cases, from the theory.
			//First case: Horizontal. The degree is included between 0-22.5 OR 157.5-180 -> The degree should be set to 0° (Horizontal = Left or Right)
			if(degree >= 0 && degree < 22.5 || degree >= 157.5 && degree <= 180) {
				degree = 0;
			}
			//Second case: Positive Diagonal. (Upper Left to Lower Right) The degree is included between 22.5-67.5 -> The degree should be set to 45°.
			else if (degree >= 22.5 && degree < 67.5) {
				degree = 45;
			}
			//Third case: Vertical. The degree is included between 67.5-112.5 -> The degree should be set to 90° (Vertical: Upper to Bottom)
			else if (degree >= 67.5 && degree < 112.5) {
				degree = 90;
			}
			//Fourth case: Negative Diagonal. The degree is included between 112.5-157.5. The degree should be set to 135°.
			else if(degree >= 112.5 && degree < 157.5) {
				degree = 135;
			}
			//else: something went wrong, display it. (With the normalization this should never happens).
			else {
				cerr << "Something went wrong." << endl;
			}
			
			//assign the pixel
			dest_image.at<uchar>(i-pad, j-pad) = pixel_intensity;
			gradient_orientations.at<uchar>(i-pad, j-pad) = degree;
			
			
			/* To assign the pixel to the DX and DY images
			DX.at<uchar>(i-pad, j-pad) = gradientX;
			DY.at<uchar>(i-pad, j-pad) = gradientY;
			*/
			
			//reset the values
			pixel_intensity = 0;
			gradientX = 0;
			gradientY = 0;
			theta = 0.0;
			degree = 0.0;
		
		}
	}
	
	/* Show the derivate X / Y images
	imshow("DX", DX),
	imshow("DY", DY);
	*/
	
	waitKey(0);

}

void Padding(Mat raw_image, Mat& padded_image, int Ksize) {

	int pad = floor(Ksize/2);
	padded_image = Mat::zeros(raw_image.rows+2*pad, raw_image.cols+2*pad, raw_image.type());
	
	for(int i=0; i<raw_image.rows; i++) {
		for(int j=0; j<raw_image.cols; j++) {
			padded_image.at<uchar>(i+pad, j+pad) = raw_image.at<uchar>(i,j);
		}
	}
	
}

void Blur(Mat padded_image, Mat kernel, Mat& dest_image) {

	int pad = floor(kernel.rows/2);
	double pixel_value = 0.0;
	dest_image = Mat::zeros(padded_image.rows-2*pad, padded_image.cols-2*pad, padded_image.type());
	
	for(int i=pad; i<padded_image.rows-pad; i++) {
		for(int j=pad; j<padded_image.cols-pad; j++) {
			for(int u=0; u<kernel.rows; u++) {
				for(int v=0; v<kernel.cols; v++) {
					pixel_value += padded_image.at<uchar>(i-pad+u, j-pad+v) * kernel.at<double>(u,v);
				}
			}
			
			dest_image.at<uchar>(i-pad, j-pad) = floor(pixel_value);
			pixel_value = 0.0;
			
		}
	}

}

Mat Canny(Mat raw_image, Mat& padded_image, Mat& kernel, Mat& dest_image) {

	//Canny Step 1: Smooth the image from eventual noise
	//A. Create a Gaussian Kernel (Non-Aggressive Blur: 5x5 with sigma^2 = 1.4 (1.20^2 = 1.44)
	Gaussian(kernel, 1.2);
	
	//B. Pad the image to make the convolution possible
	Padding(raw_image, padded_image, kernel.rows);
	
	//C. Blur the image with the Gaussian Kernel
	Blur(padded_image, kernel, dest_image);
	
	//D. Show the resultant image
	imshow("Gaussian Image", dest_image);
	waitKey(0);
	
	//Canny Step 2: Calculate the gradient and the orientation of the image
	//A. We need to create a new empty Mat object to fill with the information of the 
	//   Gaussian Image, because we need to apply the sobel operator on it. So let's
	//   create a new Padded Gaussian Image Mat object and let's pad it. Also create
	//   a new Mat object for the Sobel result.
	Mat padded_gaussian_image, dest_sobel_image;
	Padding(dest_image, padded_gaussian_image, kernel.rows);
	
	//B. Now let's calculate sobel. With sobel we'll know about the gradient intensity and orientation,
	//   and for that we need a separate Mat object to take track of all the changing happening to the direction.
	int max_intensity = 0;
	Mat gradient_orientations;
	Sobel(padded_image, dest_sobel_image, gradient_orientations, max_intensity);
	
	//C. Show the sobel image
	imshow("Sobel Image", dest_sobel_image);
	waitKey(0);

	//Canny Step 3: Non-Maximum Suppression
	//Now that we got the gradient orientation obtained by expressing the degree of the theta angle to a "comprensible"
	//orientation for the image into the sobel algorithm, we now need to suppress all the non-maximum values by confronting 
	//every pixel centered into the image with his orientation neighbour pixels. We can do that based on the case we're analyzing.
	//We need then to iterate through all the sobel image to suppress every pixel that is not the maximum.
	for(int i=0; i<dest_sobel_image.rows; i++) {
		for(int j=0; j<dest_sobel_image.cols; j++) {
			//find the gradient orientation neighbours values
			int N = dest_sobel_image.at<uchar>(i-1, j);
			int S = dest_sobel_image.at<uchar>(i+1, j);
			int E = dest_sobel_image.at<uchar>(i, j+1);
			int W = dest_sobel_image.at<uchar>(i, j-1);
			int NE = dest_sobel_image.at<uchar>(i-1, j+1);
			int NW = dest_sobel_image.at<uchar>(i-1, j-1);
			int SW = dest_sobel_image.at<uchar>(i+1, j-1);
			int SE = dest_sobel_image.at<uchar>(i+1, j+1);
			int direction_choice = gradient_orientations.at<uchar>(i,j);
			int current_pixel = dest_sobel_image.at<uchar>(i,j);
			
			//switch the choice, and based on which angle we'll suppress the pixel in the right neighbour area.
			switch(direction_choice) {
				//Case if direction is horizontal
				case 0:
					if(current_pixel < E || current_pixel < W)
						dest_sobel_image.at<uchar>(i,j) = 0;
					break;
				
				//Case if direction is into the positive diagonal (upper left to bottom right)
				case 45:
					if(current_pixel < NW || current_pixel < SE) 
						dest_sobel_image.at<uchar>(i,j) = 0;
					break;
				
				//Case if direction is vertical 
				case 90:
					if(current_pixel < N || current_pixel < S) 
						dest_sobel_image.at<uchar>(i,j) = 0;
					break;
				
				//Case if direction is into the negative diagonal (upper right to bottom left)
				case 135:
					if(current_pixel < NE || current_pixel < SW) 
						dest_sobel_image.at<uchar>(i,j) = 0;
					break;
				
				default: 
					break;
			
			}
		}
	}
	
	imshow("Non-Maximum suppress Sobel", dest_sobel_image);
	waitKey(0);
	
	
	//Canny Step 4: Hysteresis
	//The process of Hysteresis is aimed to solve the "hole" created into the edges pixel to recreate a continous edge that
	//represent in the best possible way the edges of the image.
	//We can calculate the upper and the lower thresold multipling a ratio to the highest magnitude gradient pixel into the image.
	//Then we got three cases: The pixel is strong (value > upper thresold) and is certainly an edge, the pixel is weak (value < upper thresold, 
	//but > lower thresold): this could be an edge, and shoul be confronted with his neighbour to find if there is some strong pixel adjacent to him:
	//if yes, then it's a edge pixel too. The last case, the pixel is irrelevant (value < lower thresold) and should be marked to 0.
	//We also define the last matrix to hold the canny image.
	Mat result_image = Mat::zeros(dest_sobel_image.rows, dest_sobel_image.cols, dest_sobel_image.type());
	int lower_t = cvRound(max_intensity*0.05);
	int upper_t = cvRound(max_intensity*0.1);
	
		//for the size of the image
    for(int i=0; i<dest_sobel_image.rows; i++) {
        for(int j=0; j<dest_sobel_image.cols; j++) {
            //First case: check if the pixel is strong (pixel value > upper thresold)
            if(dest_sobel_image.at<uchar>(i,j) > upper_t) {
                result_image.at<uchar>(i,j) = 255;
            }
            //Second case: check if the pixel is irrelevant (pixel value < lower thresold)
            else if(dest_sobel_image.at<uchar>(i,j) < lower_t) {
                result_image.at<uchar>(i,j) = 0;
            }
            //Third case: check if the pixel is weak (pixel vale < upper thresold but > lower thresold)
            else if(dest_sobel_image.at<uchar>(i,j) <= upper_t && dest_sobel_image.at<uchar>(i,j) >= lower_t) {
            	bool has_neighbour = false;
            
            	//check all of his neighbour
                for(int u=-1; u<=1; u++) {
                	for(int v=-1; v<=1; v++) {
                		//if one of his neighbour is a pixel > t1
                		if(dest_sobel_image.at<uchar>(i+u, j+v) > upper_t) {
                			//then the current should be a strong pixel too: has a max neigbour and break
                			has_neighbour = true;
                			break;
                		} 
                	}
                }
                
                //check if one of his neighbour was max: if so, is a strong pixel. irrelevant otherwise
                if(has_neighbour) { 
                	result_image.at<uchar>(i,j) = 255;
                } else {
                	result_image.at<uchar>(i,j) = 0;
                }
                
            }
        }
    }

	//return the final image
	return result_image;
	
}
