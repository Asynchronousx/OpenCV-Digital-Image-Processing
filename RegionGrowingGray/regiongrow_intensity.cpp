#include <opencv2/opencv.hpp>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <stack>

//Region Grow: Grayscale
//The region grow algorithm is a segmentation algorithm that groups pixels into components called regions.
//The main purpose of this algorithm is assigning to each pixel a "label", that represent the region to which 
//they belong. All the pixel with the same "label", or i-th Region (Ri) share a same charateristic, and then they 
//should be grouped with the same value. The goal of the region grow algorithm is then map individual pixels to set of pixels,
//in which each set share a charateristic of the image. Each set of pixels that share the same charateristic is called Region.
//The approach of the region growing algorithm is simple:
//1) We start at a random pixel (can be the top-left pixel, the middle one or a random choosen pixel), and give that pixel an initial reigon.
//2) Once done that, we analyze all of his neighbour pixels, based on the connectivity we choose (could be 4 or 8), and add them to a stack,
//   if they're not analyzed yet. We can check the state of each pixel by utilizing an additional matrix to store all the result of the labeling,
//   and check if the correspondent pixel into the additional matrix is different than 0. If so, we add that pixel to the stack.
//3) After processing a pixel and checking its state, we should add it to the current region if it does not change much from the analyzed pixel. 
//   If the difference is too much from the current variance of the region, then that pixel should be added in a new region and marked with 
//   an incremented label.
//4) After all the pixel of neighborhood are analyzed, we should move to the next pixel by the stack order.
//5) Then we repeat the same procedure from the step 2 until all the pixel into the image have been analyzed.
//What's a good comparison term to check if the pixel should be added to a region? 
//We can check if the Absolute value between the difference of the current pixel and the variance value of his region, is lower than a threshold T
//Passed in input by the user. If abs(pixel(x,y) - Mean of the neighbor) is lower than a threshold T, then the pixel(x,y) belong to a specific region.
//If not, that pixel forms a new region, and should be labeled with a new incresed number into the additional matrix.
//Once done that, the additional matrix M contains all the labels that represent the actual number of all the regions found.
//We can store into an HashMap the pair (RegionNumber, Variance) for putting into the original image the segmented region with the correspondent value.
//So, iterating through the original image and the region matrix M, we set the pixel at Image(x,y) as the value of the correspondent region into the region matrix:
//If the region matrix(x,y) is 1 for example, we should check into the HashMap which value is represented by the N region: we then set the value of the Image(x,y) to 
//the right approssimation.
//Note: the variance for each pixel is calculated by the mean value of the region - the current pixel intensity.

using namespace std;
using namespace cv;

void RegionGrowing(Mat, Mat&, Mat&);
void grow(vector<Point>, Mat, Mat&, Mat&, int);

//Threshold and possible region of the images
int t, quality;

int main(int argc, char** argv) {

	if(argc != 4) {
		cerr << "Usage: ./<programname> <image.format> <threshold> <quality>" << endl;
		exit(EXIT_FAILURE);
	}
	
	Mat raw_image = imread(argv[1], IMREAD_GRAYSCALE);

	if (raw_image.empty()) {
		cerr << "Image format must be valid." << endl;
		exit(EXIT_FAILURE);
	}
	
	//Declaring our destination image.
	Mat dest_image(raw_image.size(), raw_image.type(), Scalar(0));
	
	//Declaring the additional M matrix that represent the regions: same size of the input image, filled with 0.
	//We declare this float because we need to take track only the region numbers.
	Mat regions = dest_image.clone();
	
	//init the threshold and the quality of the image: since the quality represent the iterations, we're going to mulitply
	//the inserted value with 1000. (i.e: value = 10 -> 100.000 iterations)
	t = atoi(argv[2]);
	quality = atoi(argv[3])*1000;
	
	//Region Growing algorithm: it need the input and output image, and the regions matrix.
	RegionGrowing(raw_image, regions, dest_image);
	
	imshow("Original", raw_image);
	imshow("Region Grow", dest_image);
	waitKey(0);

}


void RegionGrowing(Mat raw_image, Mat& regions, Mat& dest_image) {

	//The Region Growing algorithm is the helper-function that aim to find all the possible regions into an image,
	//given an initial threshold of intensity variance, and the possible regions present into the image.
	
	//Now, since we're working into the Non-Supervisioned, we must choose the accuracy for which we can found a possible region.
	//The higher the seed generated, the higher the quality.
	//That's because we're using a random approach to find the regions. In 1000 iterations of random generated seeds, we can found,
	//for example, 100 possible different seeds not visited. In 10000 iterations we can found 1000 different seed not visited and so on.
	//So, the higher the iterations, the higher the quality of the resolution of the output image.

	//Why this actually happens? Because the initial image is filled with 0, at the beginning of the segmentation the image results total black.
	//With the iterations we then proceed to find a new region everytime we extract a pixel. We do this randomically.
	//So, at the first iteration the image is total black, and we must start to find a region of the image. We extract a random
	//pixel contained between rows and cols of the image, and proceed to expand that region starting from the random seed. 
	//At the end of the first iteration, we found the first segmented region of the image. We still need to iterate N-1 times,
	//(where N is the number of the quality we want into the output image), to found the full segmented image.
	//Note that, since the approach is random, it could happen that the region will be segmented and ready to be displayed in the frist 
	//few iteration, resulting in the waste of computational power in the last N iteration remaining.
	//Or it could happen that it will find the complete segmented regions in the last few iterations. The approach is random.
		
	//Declaring a vector of Points to take track of all the pixel pushed in every iteration
	//of the region grow algorithm.
	vector<Point> pixelVector;
	
	//Init a seed for the x and y position, and a pixel value in that image(x,y). 
	int x_seed, y_seed, pixel_value;
	
	//For every possible region into the image
	for(int i=0; i<quality; i++) {
		//x and y seed are choosen random between 0 and max rows/cols
		x_seed = rand()%raw_image.rows;
		y_seed = rand()%raw_image.cols;
		
		//What we're doing here, is pretty straightforward: if the extracted random pixel was already visited
		//into a previous iteration, we don't need to analyze it because it was already done in the past. This is
		//the necessary condition that allow us to skip useless iterations of the grow algorithm. 
		//If regions(x,y) is not equal to zero, that means that the area starting at that seed is already drawn into the output image.
		if(regions.at<uchar>(x_seed, y_seed) == 0) {
		
			//Since the pixel at image(seed x, seed y) passed the if, that means that there isn't an analyzed pixel in the coordinate
			//(seed x, seed y) of the region matrix. That means that we can expand this pixel and find his entire region.
			//We mark it as visited so it will be not be considered into the growing algorithm.
			regions.at<uchar>(x_seed, y_seed) = 1;
			
			//Assigning the pixel of the image at (x seed, y seed) to a pixel value variable, that will be useful to find all the
			//pixel with similar intensity to this one, that represent the first seed of the current region.
			pixel_value = raw_image.at<uchar>(x_seed, y_seed);
			
			//Since this pixel will represent the entire region, and we've already marked that pixel with "visited", we manually assign
			//this pixel to the dest output image. Otherwise, if not doing this, this pixel will be not considered.
			dest_image.at<uchar>(x_seed, y_seed) = pixel_value;		
			
			//Clear the vector from previous iterations: we're starting at a new seed.
			pixelVector.clear();
			
			//Proceed to push the pixel as a point into the pixel's vector: this will be the first extracted pixel.
			pixelVector.push_back(Point(x_seed, y_seed));
			
			//Calling the grow function to let that pixel grow in that region
			grow(pixelVector, raw_image, regions, dest_image, pixel_value);			
		} 
	}
	
	//Applying median blur to remove the noise: since the initial matrix was filled with 0, is possible that
	//with low quality (and low iterations) the result image could contains salt and pepper noise caused by some 
	//black pixel remaining into the output image. Applying a median blur solve that.
	//Note: with the increment of the median blur, the image should result more clear and segmented.
	medianBlur(dest_image, dest_image, 3);

}

void grow(vector<Point> pixelVector, Mat raw, Mat& regions, Mat& dest, int pixel_value) {		
	//Variables: a pixel representing the point at (x,y) and his coordinate
	Point pixel;
	int x, y;	
	
	//Until the vector is empty (that means, until all the pixel of a region hasn't been analyzed yet)
	while(!pixelVector.empty()) {
	
		//extracting a pixel from the vector: we're using a vector as a stack (using the back of the vector as the top)
		//and fetching his coordinates
		pixel = pixelVector.back();
		x = pixel.x;
		y = pixel.y;
		
		//pop the pixel from the vector: we're done with that pixel because we're gonna analyze it now
		pixelVector.pop_back();

		//Simulating a sliding kernel with the 8-Connectivity
		for(int u=-1; u<=1; u++) {
			for(int v=-1; v<=1; v++) {

				//check the boundaries of the image, assuring that the current analyzed pixel into the sliding kernel does not
				//exceed the rows and cols
				if(x+u >=0 && x+u < raw.rows && y+v >=0 && y+v < raw.cols) {
								
					//Instead of considering the variance of a pixel, gave by |image(x+u, y+v) - average of that region|, we're using a 
					//more soft-approach: we're considering how much a pixel intensity differs from the seed pixel that represent the region,
					//passed in input as pixel_value, the pixel at position image(x seed, y seed).
					//if |seed pixel - current pixel| is lower than a threshold, this means that the intensity difference is not that much:
					//we can assure that the current pixel belong to that region. 
					//If is BIGGER than the threshold, we know for sure that the pixel should not belong to that region because differs too much.
					//We also check that the current pixel has not been analyzed yet (region(current pixel) == 0)
					if(abs(pixel_value - raw.at<uchar>(x+u, y+v)) < t && regions.at<uchar>(x+u, y+v) == 0) {
					
						//We mark the pixel as visited
						regions.at<uchar>(x+u, y+v) = 1;
						
						//Assign to the current pixel the intensity of the seed pixel that represent the region (approssimating the value of the region itself)
						dest.at<uchar>(x+u, y+v) = pixel_value;	
						
						//And pushing the pixel into the vector, this will be useful to find other pixels belonging to the current region.
						pixelVector.push_back(Point(x+u, y+v));
					
					} 
				}
			}
		}
	}	
	
}
