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

class Average3b {
	public:
	double b;
	double g;
	double r;
	int pixel_num;
	
	Average3b(double b, double g, double r, int pixel_num) {
		this->b = b;
		this->g = g;
		this->r = r;	
		this->pixel_num = pixel_num;
	}
	
	Average3b() {;}
	
};

void RegionGrowing(Mat, Mat&, Mat&);
void grow(vector<Point>, map<int, Average3b>, Mat, Mat&, Mat&, Vec3b, int);

//Threshold and possible region of the images
int t, quality;

int main(int argc, char** argv) {

	if(argc != 4) {
		cerr << "Usage: ./<programname> <image.format> <threshold> <quality>" << endl;
		exit(EXIT_FAILURE);
	}
			
	Mat raw_image = imread(argv[1], IMREAD_COLOR);
	
	if (raw_image.empty()) {
		cerr << "Image format must be valid." << endl;
		exit(EXIT_FAILURE);
	}
	
	//OPTIONAL: apply a smoothing filter to omogenize the regions into the image
	GaussianBlur(raw_image, raw_image, Size(3,3), 1, 1);
	
	//Declaring our destination image.
	Mat dest_image(raw_image.size(), raw_image.type(), Scalar(0));
	
	//Declaring the additional M matrix that represent the regions: same size of the input image, filled with 0.
	//We declare this float because we need to take track only the region numbers.
	Mat regions(raw_image.size(), CV_32F, Scalar(0));
	
	//init the threshold and the quality of the image: since the quality represent the iterations, we're going to mulitply
	//the inserted value with 1000. (i.e: value = 10 -> 100.000 iterations)
	t = atoi(argv[2]);
	quality = atoi(argv[3])*1000;
	
	//Region Growing algorithm: it need the input and output image, and the regions matrix.
	RegionGrowing(raw_image, regions, dest_image);
	
	imshow("Original", raw_image);
	imshow("Region Grow", dest_image);
	imwrite("dest.png", dest_image);
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
	
	//HashMap that represents the regions into the image: it contains, at the key the region number, 
	//And our custom class Average3b to mantain the average of the value of the three channels (B,G,R) in the region.
	//We also take track of the number of pixel in that region.
	//The average and the pixel number is useful to calculate the variance of a pixel.
	map<int, Average3b> regionMap;
	
	//Init a seed for the x and y position, and a pixel value in that image(x,y). 
	int x_seed, y_seed;
	
	//Init a pixel value: we're operating into the color space, so we need a Vec3B object.
	//Vec3B is the typedef for Vec<uchar, 3>: this means that is a simple uchar vector with 3 elements: Blue, Green and Red.
	//Because the OpenCV implements colors in the format of BGR instead of RGB.
	Vec3b pixel_value;
	
	//Initially, the current region is only one
	int current_region = 1;
	
	//For every possible region into the image
	for(int i=0; i<quality; i++) {
		//x and y seed are choosen random between 0 and max rows/cols
		x_seed = rand()%raw_image.rows;
		y_seed = rand()%raw_image.cols;
		
		//What we're doing here, is pretty straightforward: if the extracted random pixel was already visited
		//into a previous iteration, we don't need to analyze it because it was already done in the past. This is
		//the necessary condition that allow us to skip useless iterations of the grow algorithm. 
		//If regions(x,y) is not equal to zero, that means that the area starting at that seed is already drawn into the output image.
		if(regions.at<float>(x_seed, y_seed) == 0) {
		
			//If the pixel in the position (seed x, seed y) pass this test, that means that this pixel needs to be analyzed. So we didn't
			//find a proper region linked to that specific pixel. 
			//We mark it as visited with the current region number so it will be not be considered into the growing algorithm.
			regions.at<float>(x_seed, y_seed) = current_region;
			
			//Fetching the pixel value at the seed position
			pixel_value = raw_image.at<Vec3b>(x_seed, y_seed);
			
			//clearing the vector from precedent iterations (bounded to another region)
			pixelVector.clear();
			
			//pushing the pixel at seed coordinate into the vector: this will be the first pixel that will be analyzed
			pixelVector.push_back(Point(x_seed, y_seed));
			
			//insert a new region: current region, mean of the three channel of the regione (currently the BGR pixel itself)
			//and the number of pixel of that region (1).
			regionMap.insert(make_pair(current_region, Average3b(pixel_value[0], pixel_value[1], pixel_value[2], 1)));
			
			//calling the grow function to let that pixel grow in that region
			grow(pixelVector, regionMap, raw_image, regions, dest_image, pixel_value, current_region);
		
			//once done, increase the number of region: that's because once here a region has been analyzed: we must increment to let the iteration find a new
			//region in the next step.
			current_region++;
			
		} 
	}
	
	
	//Once exited the main for loop, we've analyzed all the possible regions found with that quality parameter: 
	//Once here, the region matrix will hold, in each of his pixel, a number representing the region of which that pixel belong.
	//We can easily use that pixel to access the hashmap containing, for that key (the region number), the total sum of his BGR elements
	//and the number of pixel that forms that region: with that, we can simply do average[channel]/pixelnumber to obtain the average of that region.
	
	//We then iterate for the region size (that it's the same as the dest image size)
	for(int i=0; i<regions.rows; i++) {
		for(int j=0; j<regions.cols; j++) {
			//the i-th, j-th pixel of the dest image, is equal to the average of the region map (sum/pixelnumber) obtained with
			//using the i-th, j-th pixel as the key for the hashmap to obtain that value for that specific region number.
			pixel_value[0] = round(regionMap[regions.at<float>(i,j)].b/regionMap[regions.at<float>(i,j)].pixel_num);
			pixel_value[1] = round(regionMap[regions.at<float>(i,j)].g/regionMap[regions.at<float>(i,j)].pixel_num);
			pixel_value[2] = round(regionMap[regions.at<float>(i,j)].r/regionMap[regions.at<float>(i,j)].pixel_num);
			dest_image.at<Vec3b>(i,j) = pixel_value;
		}	
	}
	
	//Applying median blur to remove the noise: since the initial matrix was filled with 0, is possible that
	//with low quality (and low iterations) the result image could contains salt and pepper noise caused by some 
	//black pixel remaining into the output image. Applying a median blur solve that.
	//Note: with the increment of the median blur, the image should result more clear and better segmented.
	medianBlur(dest_image, dest_image, 3);
	
	
	
}

void grow(vector<Point> pixelVector, map<int, Average3b> regionMap, Mat raw, Mat& regions, Mat& dest, Vec3b pixel_value, int region) {		
	//Variables: a pixel representing the point at (x,y) and his coordinate.
	//We also use a Vec3b value to calculate at each iteration the average of the three channels.
	Point pixel;
	int x, y;
	double variance;
	Vec3b average_pixel;
	
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
					
					//Calculate the variance: in the color space, we can calculate the variance with the euclidean distance:
					//The euclidean distance, given the two ternary couple (x,y,z) that represent the current pixel, and (a,b,c) that represent
					//the average of the region, is calculated by:
					//SQRT((x-a)^2 + (y-b)^2 + (z-c)^2). x, y and z represent the Blue, Green and Red value of the current pixel analyzed.
					
					//assigning to the average pixel the average of that region computed until now
					average_pixel[0] = round(regionMap[region].b/regionMap[region].pixel_num);
					average_pixel[1] = round(regionMap[region].g/regionMap[region].pixel_num);
					average_pixel[2] = round(regionMap[region].r/regionMap[region].pixel_num);
					
					//using the opencv norm function to calculate the euclidean distance
					variance = norm(average_pixel, raw.at<Vec3b>(x+u, y+v));
								
					//If the variance of this pixel is contained into the threshold and the analyzed pixel has not been visited yet
					//(region(x+u, y+v) marked as 0
					if(variance < t && regions.at<float>(x+u, y+v) == 0) {
										
						//Mark the pixel as visited with the current region passed in input
						regions.at<float>(x+u, y+v) = region;
						
						//push the new pixel into the stack: we're going to analyze it later
						pixelVector.push_back(Point(x+u, y+v));
						
						//update the average of the region adding to the current value the value of the current pixel (of the sliding windows)
						//for all the three channels: B, G, R.
						regionMap[region].b += raw.at<Vec3b>(x+u, y+v)[0];
						regionMap[region].g += raw.at<Vec3b>(x+u, y+v)[1];
						regionMap[region].r += raw.at<Vec3b>(x+u, y+v)[2];
						
						//update the pixel number in that region
						regionMap[region].pixel_num++;
															
					} 
				}
			}
		}
	}	
	
}
