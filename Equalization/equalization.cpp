#include <opencv2/opencv.hpp>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#define L 256

//rk -> k-th gray scale level
//nk -> number of occurrence of the k-th rk
//Pr(rk) -> probability of the rk-th grayscale level pixel to appear 
//			Pr(rk) is equal to: nk/MxN, where M and N size of the image.
//From the theory, the eualization is given by: L-1(Summatory from j=0 to k)Pr(rk).

using namespace std;
using namespace cv;

void calculate_occurrences(Mat&, vector<int>&);
void calculate_probabilities(vector<int>, vector<double>&, int size);
void calculate_cumulative_probabilities(vector<double>, vector<double>&);
void equalize(Mat&, Mat&, vector<double>);

int main(int argc, char** argv) {

	srand(time(NULL));

	if(argc!=2) {
		cerr << "Usage: ./<programname> <imagename>.<format>" << endl;
		exit(EXIT_FAILURE);
		
	}
	
	Mat raw_image = imread(argv[1], CV_8UC1); //IMREAD_GRAYSCALE
	
	if(raw_image.empty()) {
		cerr << "Image must be valid." << endl;
		exit(EXIT_FAILURE);
		
	}	
	
	//Declare two vector for Nk (occurrences) and Pr(rk) (probabilities)
	vector<int> occurrences(L-1);
	vector<double> probabilities(L-1);
	vector<double> cumulative_probabilities(L-1);
	Mat dest_image = Mat::zeros(raw_image.rows, raw_image.cols, raw_image.type());
	Mat openCV_dest_image(dest_image.clone());
	
	//calculate the occurrences of all the rk into the image, and their probability
	calculate_occurrences(raw_image, occurrences);
	calculate_probabilities(occurrences, probabilities, raw_image.rows*raw_image.cols);
	calculate_cumulative_probabilities(probabilities, cumulative_probabilities);
	//Calling custom equalization
	equalize(raw_image, dest_image, cumulative_probabilities);
	
	//Calling openCV equalization
	equalizeHist(raw_image, openCV_dest_image);
	
	imshow("Original Image", raw_image);
	imshow("Custom Equalization", dest_image);
	imshow("OpenCV Equalization", openCV_dest_image);
	
	waitKey(0);
	return(0);
	
}

void calculate_occurrences(Mat& raw_image, vector<int>& occurrences) {

	//Calculating all the occurrences of a given rk.
	for(int i=0; i<raw_image.rows; i++) {
		for(int j=0; j<raw_image.cols; j++) {
			//Just increase the relative scale level into the vector
			occurrences.at((int)raw_image.at<uchar>(i,j))++;
			
		}
	}

}

void calculate_probabilities(vector<int> occurrences, vector<double>& probabilities, int size) {

	//Calculating all the probabilities of a given rk to appear into the image.
	for(int i=0; i<L-1; i++) {
		//Probability is given by: Pr(rk) = nk/MxN
		probabilities.at(i) = (double)occurrences.at(i)/size;
		
	}

}

void calculate_cumulative_probabilities(vector<double> probabilities, vector<double>& cumulative_probabilities) {

	//Cumulate all the probabilities from the formula L-1(Summatory from j=0 to k)Pr(rk)
	double pixel_value = 0;
	for(int i=0; i<probabilities.size(); i++) {
		pixel_value += (L-1)*(probabilities.at(i));
		cumulative_probabilities.at(i) = pixel_value;
		
	}

}

void equalize(Mat& raw_image, Mat& dest_image, vector<double> cumulative_probabilities) {

	//now insert the equalized pixels into the destination image
	for(int i=0; i<raw_image.rows; i++) {
		for(int j=0; j<raw_image.cols; j++) {
			dest_image.at<uchar>(i,j) = floor(cumulative_probabilities.at((int)raw_image.at<uchar>(i,j)));
			
		}		
	}
	
}

