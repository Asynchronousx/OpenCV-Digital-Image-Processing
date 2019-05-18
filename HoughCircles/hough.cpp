#include <opencv2/opencv.hpp>
#include <iostream>
#include <math.h>
#define RADCOEF CV_PI/180

//Hough Transform: Circles
//The hough transform algorithm is meant to find, given a source image and a certain threshold, 
//all the possible circles into the image. The hough transform is particoular insensible to noise and occlusions.
//Why hough TRANSFORM? Because when operating with Hough, we're analyzing a secondary geometric space called
//PARAMETER SPACE (SP). 
//In the regoular cartesian space, we operate on coordinates (X,Y), and the space where those coordinates works is called
//IMAGE SPACE (SI). In the case of circles, where a circle is expressed by 3 parameters (Center(x,y) and Radius), we got a 
//similalry space called PARAMETER SPACE (SP), where those circles are expressed as points.
//With that said, given a point into the SI, we can identify all the possible circles that fall into that point (X,Y) into the SP,
//with a THETA variation based on the angle analyzed. That's because we're DISCRETIZING the parameter space to host a limited number of corrispondences.
//Similalry, given a point into the SP of coordinate (A,B), we can identify only ONE circonference into SI.
//Then, how do we find a potential circle into the SI? 
//First of all, we need to DISCRETIZE our parameter space, because the circles found could be potentially infinite.
//So, we transfer the calculus from the geometric space to a POLAR SPACE.
//We can identify a point into the SP by the formula:
//a = x - Rcos(THETA), where R is the radius and x is the x-axis coordinate of the center.
//b = y - Rsin(THETA), as above, y is the coordinate of the y-axis.
//Similalry, we can found a SI coordinate by
//x = b + Rcos(THETA)
//y = b + Rsin(THETA)
//THETA represent the variation of the angle of a circle.
//THETA then varies between a certain + or - degree. For example, a circle could vary between 0-360°.
//Having a finite space, we must implement a mechanism that could identify a potential circle into our image. This mechanis is called
//VOTING PROCESS, and the logic behind it is simple: Given the entire SP, we must find WHICH point into the SP calculated through Hough has the maximum occurrences.
//Thats because, from the theory, we know that a point into the SP is surely a circle into the SI if his circle intersection number is maximum into the SP.
//That's works for ONE circle. For more than one line, we must set a threshold to let pass all the point with a given number of occurrence > than that threshold.
//Also, we need to set a MAXIMUM RADIUS and a MINIMUM RADIUS, to consider only the circle between that radius.
//Note: before operating with Hough we must works on border: let's do a reasonable Gaussian Blurring and a Canny edge detection.


using namespace std;
using namespace cv;

//Class that represent a circle point, composed by:
//Center (of x and y coordinate into the SI) and Radius.
class CirclePoint {
	public:
	int x;
	int y;
	double radius;
	
	CirclePoint(int x, int y, int radius) {
		this->x = x;
		this->y = y;
		this->radius = radius;
	}
};

void HoughTransform(Mat, Mat&);

//Global variables to take track of threshold
int hough_threshold;
int max_radius;
int min_radius;


int main(int argc, char** argv) {

	if(argc!=5) {
		cerr << "Usage: ./<programname> <imagename.format> <threshold> <min radius> <max radius>" << endl;
		exit(EXIT_FAILURE);
	}
	
	Mat raw_image = imread(argv[1], IMREAD_GRAYSCALE);
	
	if(raw_image.empty()) {
		cerr << "Image format not valid." << endl;
		exit(EXIT_FAILURE);
	}
	
	//assign the boundaries
	hough_threshold = atoi(argv[2]);
	min_radius = atoi(argv[3]);
	max_radius = atoi(argv[4]);
	
	Mat blur_image, edge_image;
	Mat dest_image = raw_image.clone();

	//Applying gaussian blur on the source and finding the edges with Canny
	equalizeHist(raw_image, raw_image);
	GaussianBlur(raw_image, blur_image, Size(5,5), 1.5, 1.5);
	Canny(blur_image, edge_image, 80, 150, 3);
	
	//showing the result
	imshow("Blurred Image", blur_image);
	imshow("Canny Image", edge_image);
	waitKey(0);
		
	//Applying hough transform
	HoughTransform(edge_image, dest_image);
	
	//Showing the image
	imshow("Hough Circles", dest_image);
	imwrite("Houhg.png", dest_image);
	waitKey(0);
	return(0);

}

void HoughTransform(Mat edge_image, Mat& dest_image) {

	//Defining the votes matrix, that represent the parameter space SP in which
	//we'll draw any of the possible circles that belongs to the SI.
	//The size of this matrix should be either 2D or 3Dimensional.
	//From the theory, a 3D matrix should be used only when the radius is unknown.
	//Here we got an interval of radius that space from min radius to max radius.
	//To found the right interval and draw the correct circles, we must analyze all of those intervals.
	//We're going to use a 3D matrix to take advantage of his 3-dimensional structure to memorize, at every
	//iteration, if a point into the SP (a,b) should be voted or not. With the 3-dimension we're going to take 
	//track of a certain radius, a certain a coordinate a (x into the polar space) and b coordinate (y into the polar space).
	//That's comes handy when we need to memorize a certain vote into the votes matrix at A(r,a,b).
	//How we're going to init our matrix?
	
	//First, init a int pointer: we can use an array of three element declared inline.
	//At the first value, we're going to set all the POSSIBLE values that float between min and max radius.
	//For example, if min radius is 10 and max radius is 20, the possible analyzable elements are 11. (10...20 inclusive)
	//So, (max radius - min radius) + 1 return exactly 11 possible radius to analyze in this specific case.
	//We're going to pass rows and cols of the edged image to represent the space of the rows and the columns.
	int values3D[] = {(max_radius-min_radius)+1, edge_image.rows+1, edge_image.cols+1};
	
	//Once done that, we're declaring a Mat object of 3 spaces (3D) composed by values3D fields. Each space
	//will contain the right field. (Mat[0] -> Radius of dimension 10, Mat[1] -> Rows of dimension M, Mat[2] -> Cols of dimension N)
	Mat votes3D = Mat(3, values3D, IMREAD_GRAYSCALE, Scalar(0));
	
	//Also, initialize a 2D matrix to take track of the resultant votes image.
	Mat votes = (Mat_<uchar>(edge_image.rows+1, edge_image.cols+1));
	votes = Mat::zeros(votes.size(), votes.type());
		
	//Let's define our polar coordinate for the circle, a and b.
	double a, b;

	//For every pixel in the rows and the columns
	for(int x=0; x<edge_image.rows; x++) {
		for(int y=0; y<edge_image.cols; y++) {
			//if the current pixel is an edge pixel (or a good candidate to be one)
			if(edge_image.at<uchar>(x,y) > 250) {
				//for every possible radius, floating from min radius to max radius passed in input
				for(int r=min_radius; r<=max_radius; r++) { 
					//for every possible theta for that specific pixel I(x,y)
					for(int t=0; t<360; t++) {
						
						//calculate the polar coordinate of that point: this formula is given by
						//a = x - R*cos(theta), b = y - R*sin(theta).
						//Dont forget to INVERT the coordinate to track the right point into the image, because an 
						//image got x and y inverted.
						a = cvRound(y - r*cos(t*RADCOEF));
						b = cvRound(x - r*sin(t*RADCOEF));
						
						//check the boundaries: if a or b exceed, skip this 
						if(a > 0 && a < edge_image.rows && b > 0 && b < edge_image.cols) {
							//update the 2D matrix at (a,b)
							votes.at<uchar>(round(a), round(b))++;
							
							//update the 3d matrix at current radius, a and b.
							votes3D.at<uchar>(r,a,b)++;
						}
					}				
				}			
			}				
		}		
	}
	
	//Show the votes
	imshow("Votes", votes);
	waitKey(0);	
	
	//Convert the image from GRAYSCALE to BGR. This will we useful to draw colorful circles into the image.
	cvtColor(dest_image, dest_image, COLOR_GRAY2BGR);
	
	//Declaring a vector of custom points 
	vector<CirclePoint> circles;
	
	//starting by iterating every possible radius from min radius to max radius
	for(int r=min_radius; r<=max_radius; r++) {
		//for every pixel (x,y) into the image
		for(int x=0; x<edge_image.rows; x++) {
			for(int y=0; y<edge_image.cols; y++) {
				
				//We're now analyzing the votes matrix.
				//If the pixel (x,y) at the current analyzed radius is major than the threshold
				//that means in the voting process was marked multiple times as a possible circumference.
				//Do not convert to polar coordinate (we can draw the circle with that coordinate, since the parameter space
				//is equal to the size of the image space) and push back a new point into the circles vector, composed by
				//the actual x, y and the radius analyzed.
				if(votes3D.at<uchar>(r,x,y) >= hough_threshold) {
					circles.push_back(CirclePoint(x,y,r));
				}
			}
		}
	}
	
	
	cout << "Possible circles found: " << circles.size() << endl;
	
	//Let's draw the circles
	for(int i=0; i<circles.size(); i++) {	
		//fetch the center and the radius
		Point center(circles.at(i).x, circles.at(i).y);
		int radius = circles.at(i).radius;
		
		//draw the center first with a green dot
		circle(dest_image, center, 3, Scalar(0,255,0));		
		
		//draw the real circumference thanks to the radius
		circle(dest_image, center, radius, Scalar(0,0,255), 2);
	}
	
}
