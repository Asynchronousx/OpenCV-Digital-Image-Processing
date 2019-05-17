#include <opencv2/opencv.hpp>
#include <iostream>
#include <math.h>
#define RADCOEF CV_PI/180

//Hough Transform:
//The hough transform algorithm is meant to find, given a source image and a certain threshold, 
//all the possible lines into the image. The hough transform is particoular insensible to noise and occlusions.
//Why hough TRANSFORM? Because when operating with Hough, we're analyzing a secondary geometric space called
//PARAMETER SPACE (SP). If a regoular cartesian space is formed by coordinates expressed in (X,Y) (called Image Space), 
//The parameter space operate on a coordinate space of (M,Q), where M is the Slope (Angular Coefficent) and Q the intercept
//of a given line into the image space (SI).
//With that said, given a point into the SI, we can identify all the possible lines that fall into that point (X,Y) into the SP,
//that are potentially infinite. Similalry, given a point into the SP of coordinate (M,Q), we can identify only ONE lines into the SI.
//Then, how do we find a potential line into the SI? 
//First of all, we need to DISCRETIZE our parameter space, because M,Q could be potentially infinite, and we can't calculate every occurrence
//of a line into a given point (they're infinite). So, we transfer the calculus from the geometric space to a POLAR SPACE.
//Coordinate (M,Q) become RO and THETA. So, if M,Q were infinite, RO and THETA are now finite, and began to the interval:
// 0 < RO < R*sqrt(2), -PI <= THETA <= PI
//Also, RO is calculated by: RO = x*Cos(THETA) + y*Sin(THETA), where x and y are the coordinate of a given point, and 
//THETA represent the variation of the angle of a line. 
//THETA then varies between a certain + or - degree. For example 90° should be a nice threshold.
//where R is the max distance that a line could have (called also Support axis, represented by the Maximum distance
//from the origin and the axis perpendicular to the lines where a given point is.
//And theta belong to the interval -PI to PI.
//Having a finite space, we must implement a mechanism that could identify a potential line into our image. This mechanis is called
//VOTING PROCESS, and the logic behind it is simple: Given the entire SP, we must find WHICH point calculated through Hough has the maximum occurrences.
//Thats because, from the theory, we know that a point into the SP is surely a line into the SI if his lines intersection number is maximum into the SP.
//That's works for ONE line. For more than one line, we must set a threshold to let pass all the point with a given number of occurrence > than that threshold.
//Retrieved all those point into the SP that are > than a threshold, we must take their coordinate RO and THETA, convert them to Cartesian and find the line
//with that M,Q into the original image. That done, we found all the lines into our image.
//Note: before operating with Hough we must works on border: let's do a reasonable Gaussian Blurring and a Canny edge detection.

using namespace std;
using namespace cv;

int hough_threshold;

void HoughTransformLine(Mat, Mat&);
void FindPolar(double, double, Point&, Point&);

int main(int argc, char** argv) {

	if(argc!=3) {
		cerr << "Usage: ./<programname> <imagename.format> <threshold>" << endl;
		exit(EXIT_FAILURE);
	}
	
	hough_threshold = atoi(argv[2]);
	Mat raw_image = imread(argv[1], IMREAD_GRAYSCALE);
	
	if(raw_image.empty()) {
		cerr << "Image format not valid." << endl;
	}
	
	Mat blurred_image, edge_image;
	Mat dest_image = raw_image.clone();
	
	//Before calling Hough, we need to apply a Gaussian Smoothing to the image, and find his edges with Canny.
	//At every step, we're displaying the result.
	GaussianBlur(raw_image, blurred_image, Size(5,5), 1.4, 1.4);
	imshow("Gaussian Blur", blurred_image);
	Canny(blurred_image, edge_image, 60, 150, 3);
	imshow("Canny Edge", edge_image);
	waitKey(0);
	
	//Applying Hough transform
	HoughTransformLine(edge_image, dest_image);
	
	//showing the result
	imshow("Hough Lines", dest_image);
	imwrite("Hough.png", dest_image);
	waitKey(0);
	return(0);

}

void HoughTransformLine(Mat edge_image, Mat& dest_image) {
	
	//The hough transform function.
	//First step is the defining how long R can be: we can assure that the longest line
	//can be drawn from the upper left/right angle to the bottom right/left. 
	//So we can safely assume that the max distance a line can obtain is the DIAGONAL of the image.
	//From the pythagora's theorem, we can find a diagonal by using the formula:
	//sqrt(side1^2 + side2^2), that can be simplified into: side*sqrt(2).
	//We'll take the maximum of the two side as preferred side.
	int max_distance = max(edge_image.rows, edge_image.cols) * sqrt(2);

	//Let's then define the votes matrix. This matrix represent our Parameter Space, and should be built in a way
	//that all the possible lines can be drawn in it.
	//Since a point (x,y) is translated into (rho, theta), we can create a matrix of size (MaxRho, MaxTheta).
	//We're using Max Rho and Max Theta because the max values that both of them can be is the size of the matrix itself.
	//This will be useful for accumulate the votes for each point valued from -pi to pi (-90 to 90) and a certain rho.
	//We already have the Max value the lines (represented by rho) can be (max_distance).
	//and also we know that the lines could vary from -90 to 90, so with 180 values max.
	//Also, to take track of the rho in the finding of the polar coordinate, we must create a mechanism that permit us
	//to know EXACTLY which rho was calculated at some point. To do this, at every compute of
	//rho = xcos(theta) + ysin(theta), we'll add max_distance value to each rho to have a costant added tho this value.
	//So, when we'll going to find the rho of each polar coordinate of the votes matrix, we'll just need to subtract the 
	//max_distance to that rho to have the exact value. Let's create a matrix DOUBLE the size of max_distance for the rows
	//that will hold all the values without any segmentation fault.
	//Let's create the matrix with (max_distance, 180).
	Mat votes = (Mat_<uchar>(max_distance*2, 180));
	
	//Populate the votes with 0
	votes = Mat::zeros(votes.size(), CV_8U);
	
	//Let's define rho and theta for every pixel: it's IMPORTANT to define them DOUBLE (or float)
	//because we need to take track of every single precision movement.
	double rho, theta;
	
	//Let's find all the possible lines for a point and start to accumulates votes into the matrix.
	for(int x=0; x<edge_image.rows; x++) {
		for(int y=0; y<edge_image.cols; y++) {
			//if the current pixel is an edge pixel
			if(edge_image.at<uchar>(x,y) > 250) {
				//for all the possible lines for that point (0-180)
				for(theta=0; theta<180; theta++) {			
					//let's calculate RHO: given by x*cos(theta) + y*sin(theta), where theta is the radiant
					//of a given point. Since we're analyzing the degree of a line for a given (x,y), we need to
					//convert theta from degree to radiant. We can do this multiplying theta for PI/180.
					//Note: the x and y are INVERTED into the geometric space. So our y is the x, and our x is the y.
					rho = cvRound((y*cos((theta-90)*RADCOEF) + x*sin((theta-90)*RADCOEF)) +max_distance);					

					//Check: 0 < RHO < R*sqrt(2)? 
					//If RHO is < 0 add max_distance to it. if is > than the max_distance, subtract from it.
				
							
					//vote for that point into the SP votes matrix
					votes.at<uchar>(rho, theta)++;

				}
			}			
		}	
	}
	
	//Once here, all the votes were made. Display the votes
	imshow("Votes", votes);
	waitKey(0);
	
	//declaring a vector of two points that reperesent the line
	vector<pair<Point,Point>> lines;
	
	//define the actual two point P1, P2 which represent the line.
	Point P1, P2;
	
	//To have colorful lines, we must convert the image from the grayscale to the BGR contest.
	cvtColor(dest_image, dest_image, COLOR_GRAY2BGR);
	
	//TODO: Implement non-maximum line suppression?	
	
	//Let's now find the peak point that are major than a defined threshold. If the i-th point meet the condition, 
	//we must translate the polar coordinate to the cartesian one, and create the relative points to display.
	for(int r=0; r<votes.rows; r++) {
		for(int t=0; t<180; t++) {
			if(votes.at<uchar>(r,t) >= hough_threshold) {
				
				//Here, take rho as the current line inder r (that's because for a certain rho, we got a certain theta.
				//Since the theta vary from 0-180, the rho associated with it is the index itself. 
				//Also we need to take the theta as the angle of that specific line at that index: so we'll multiply the 
				//theta for the radiant coefficient.
				//We also need to dispatch the current rho subtracting the max_distance to make it suitable again, and subtract
				//90 from theta to let him be in the interval -90 / 90.
				rho = r-max_distance;
				theta = (t-90)*RADCOEF;
				
				//find the polar coordinate given rho and theta
				FindPolar(rho, theta, P1, P2);
				
				//push the points into the vector
				lines.push_back(make_pair(P1,P2));
			
			}		
		}	
	}
	
	//let's now draw the line, analyzing the line vector
	for(int i=0; i<lines.size(); i++) {
		pair<Point, Point> coordinates = lines.at(i);
		line(dest_image, coordinates.first, coordinates.second, Scalar(0, 0, 255), 1, LINE_AA);
	}
	
	
}

void FindPolar(double rho, double theta, Point& P1, Point& P2) {
	
	//In this function we translate the polar coordinate into the cartesian one.
	//From the math theory, we know that, given a RHO and a THETA, the cartesian coordinate
	//From wikipedia, those values can be found from:
	//x = rho*cos(theta)
	//y = rho*sin(theta);
	//Also, we need to draw those lines for the ENTIRE image. So, given a precise point, we could vary
	//a line lenght using a big scale factor, such as the size of the image.
	//x0 and y0 are our starting point
	int x0 = cvRound(rho * cos(theta));
	int y0 = cvRound(rho * sin(theta));
	
	//let's build the point: we must then add and subtract the big scale value to let the line vary between 
	//the negative and positive scale factor.
	//We're also multiplying for -sin and cos to let them vary between the positive and negative cartesian space.	
	P1.x = cvRound(x0 + 1000 * (-sin(theta)));
	P1.y = cvRound(y0 + 1000 * (cos(theta)));
	P2.x = cvRound(x0 - 1000 * (-sin(theta)));
	P2.y = cvRound(y0 - 1000 * (cos(theta)));
	
}
