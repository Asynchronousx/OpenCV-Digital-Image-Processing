#include <opencv2/opencv.hpp>
#include <iostream>
#include <cstdlib>
#include <math.h>

//K-Mean algorithm: a first approach to image clustering
//The K-Mean is an algorithm for the clusterization of an image, belonging to
//the partitive clustering method, in which we do not classify a hierarchy of cluster,
//but we simply divide pixel into dijoint clusters by analyzing the cluster property and mean.
//To begin the algorithm we must first look at the image and wonder: how many cluster i can see
//in that image? The answer is the K of the K-Mean, and that means that we choose how many cluster
//the image should have. Be careful to find a good amount of cluster that isn't neither too big or
//too small. All the pixel will be divided into K-Cluster.
//A simple approach to the K-Means would be: given K input cluster, choose K random pixel and from 
//that pixels, expand the cluster.
//We can somehow categorize the algorithm in this two step:
//1) We confront every pixel of the image to each cluster by analyzing their DISTANCE (careful:
//with distance we do not mean the distance in coordinate from a pixel to another, but the chromatic 
//distance of the colour, calculated easily with the Euclidean distance in the RGB space).
//With that approach we start analyzing this K cluster, initially the single pixel that form the cluster is intended to be 
//the center of the cluster itself (center of the cluster = average of all the pixels presents into the image).
//And confronting each pixel to every cluster to find the best possible cluster to put the analyzed pixel in.
//So we put the pixels into the cluster that dist less from the cluster itself in terms of chromatic distance.

//2)After we put all the pixel into their relative cluster, we recalculate the centroid of the cluster.
//We simply sum all the numeric values of the pixel into the cluster and divide this value by the pixel number. 
//Once done that, we assign the new centroid to the cluster and repeat the algorithm from the step 2, inserting
//the pixels into the new computed cluster that dist less and so on.

//3) The algorithm terminate when the number of max iteration has been met or when the centroid does not change from
//the previous one (we reached the convergence).

using namespace std;
using namespace cv;

//Class useful to take track of the average of the region.
class Average3b {

	public:
	int b;
	int g;
	int r;
	
	Average3b() {;}
	Average3b(int b, int g, int r) {
		this->b = b;
		this->g = g;
		this->r = r;
	}

};

//Class that represent our cluster of pixel. We store sensitive information such as
//the cluster pixel list, the number of pixel and the average of that cluster (the centroid).
//We also insert some method to compute the average.
class Cluster {

	public:
	int pixel_num;
	Average3b region_sum;
	Vec3b centroid;
	vector<Point> pixels;
	
	Cluster() {;}
	Cluster(Vec3b centroid, int x, int y) {
		//Initially, when a cluster is created, the number of his element is zero.
		//It will be increased when the centroid of the cluster will be added
		this->pixel_num = 0;
		
		//Initializing the region sum and the centroid as empty element containing 0
		this->region_sum = Average3b(0,0,0);
		this->centroid = Vec3b(0,0,0);
		
		//adding the pixel to the region and calculating the new center (the pixel itself)
		add_pixel(centroid, x, y);
		calculate_center();

	}
	
	void calculate_center() {
		//The centroid is computed ad the sum of the entire region (of the three channels) 
		//divived by the pixel num.
		this->centroid[0] = region_sum.b/pixel_num;
		this->centroid[1] = region_sum.g/pixel_num;
		this->centroid[2] = region_sum.r/pixel_num;	
	
	}
	
	//Function that add a pixel: for each pixel added, it's BGR values it's added to the region sum
	//of the values.
	void add_pixel(Vec3b pixel, int x, int y) {
		//Adding the value of the pixel to the total region sum
		this->region_sum.b += pixel[0];
		this->region_sum.g += pixel[1];
		this->region_sum.r += pixel[2];
	
		//pushing the pixel coordinate into the vector
		this->pixels.push_back(Point(x,y));
		this->pixel_num++;	
	}
	
	void clear() {
		//resetting the pixels vector: we're going to build it in another iteration
		this->pixels.clear();
		
		//since clear is ALWAYS called after the new compute of the centroid, we assign to the 
		//total region sum the value of the centroid in each of his 3-channel. With that, we're saying
		//that the cluster is again formed by one element.
		this->region_sum.b = centroid[0];
		this->region_sum.g = centroid[1];
		this->region_sum.r = centroid[2];
		
		//since the centroid is set, the pixel count start is 1
		this->pixel_num = 1;
		
 		//the only thing we do not reset is the centroid because it was already modified from another
 		//call into the k-means algorithm.
	}


};


void KMean(Mat, Mat&, int, int, int);

int main(int argc, char** argv) {

	if(argc!=5) {
		cerr << "Usage: ./<programname> <image.format> <K cluster> <Max Iterations> <Threshold>" << endl;
		exit(EXIT_FAILURE); 	
	}
	
	Mat raw_image = imread(argv[1], IMREAD_COLOR);
	
	if(raw_image.empty()) {
		cerr << "Image format not valid" << endl;
		exit(EXIT_FAILURE);	
	}
	
	//Creating dest image filled with 0
	Mat dest_image(raw_image.size(), raw_image.type(), Scalar(0));
	
	//Assigning k cluster and iterations max
	int clusters_num = atoi(argv[2]);
	int iterations = atoi(argv[3]);
	int t = atoi(argv[4]);
	
	//Calling K-Mean
	KMean(raw_image, dest_image, clusters_num, iterations, t);
	
	imshow("Original", raw_image);
	imshow("K-Means", dest_image);
	imwrite("KMeans.png", dest_image);
	waitKey(0);

}

void KMean(Mat raw_image, Mat& dest_image, int clusters_num, int iterations, int t) {

	//Declarations of a vector of cluster to take track of the changing of the image
	vector<Cluster> clusters;
	
	//random coordinate and pixel to extract
	int rand_x, rand_y;
	Vec3b pixel;
	
	//First step: generating K starting centroid as a random point
	for(int i=0; i<clusters_num; i++) {
		rand_x = rand()%raw_image.rows;
		rand_y = rand()%raw_image.cols;
		pixel = raw_image.at<Vec3b>(rand_x, rand_y);
		clusters.push_back(Cluster(pixel, rand_x, rand_y));	
	}
	
	//Beginning of the K-Means: we're putting all the logic of the algorithm inside a for loop: 
	//thats because, we can do a fixed number of iterations: if for some reason the algorithm should
	//converge because it's centroid do not vary anymore, we break the for and proceed to 
	for(int i=0; i<iterations; i++) {
		
		//Let's start analyzing every single pixel of the image to choose to which cluster they belong:
		//iterate through the image
		
		for(int	x=0; x<raw_image.rows; x++) {
			for(int y=0; y<raw_image.cols; y++) {
				//Now, for a generic (x,y) pixel, let's discover to which cluster the pixel belong:
				//We iterate through the clust vector and check if the distance between the pixel and the 
				//centroid is lesser than a threshold: if so, add the pixel to the k-th cluster and break the
				//cycle (because we want to assign one pixel to only one cluster).
				for(int k=0; k<clusters.size(); k++) {
					if(norm(raw_image.at<Vec3b>(x,y), clusters.at(k).centroid) < t) {
						clusters.at(k).add_pixel(raw_image.at<Vec3b>(x,y), x, y);
						break;
					}				
				}
			}
		}
		
		//After all the pixel belong to the cluster, we need, for each cluster, to compute the new centroid:
		//We also check if the new centroid differs by a significative amount from the old one:
		//if differs by a significative amount, that means we need to iterate again with the new centroid
		//value to regroup the pixel into the cluster. If it not differs, and it's true for every centroid
		//of the vector of cluster, break the cycle and go to build the clustered image.
		//We check that with a boolean value changed: if some of the centroid changed by a significative amount
		//we set the changed boolean to true, and keep calculate the new centroid for each cluster in the vector.
		bool changed = false;
		
		for(int k=0; k<clusters.size(); k++) {
			//fetching the current centroid as the old one		
			Vec3b old_centroid = clusters.at(k).centroid;
			
			//compute the new centroid
			clusters.at(k).calculate_center();
			
			//check if the new centroid differs by a significative amount to the old one
			if(norm(old_centroid, clusters.at(k).centroid) > 10) {
				changed = true;
			}	
			
		}
		
		//Now check if some centroid has changed: if so, we need to continue to iterate to find the convergence
		//of the centroid, or untill all the iterations are depleted.
		//We also check if the iterations are finished: if i+1 is equal to the maximum number of iterations,
		//that means in the next iteration we would stop the for cycle: we avoid this by breaking into the if,
		//avoiding clearing the vector of pixels of each cluster.
		//If none of them changed
		if(!changed || i+1 == iterations) {
			//break the cycle and go to build the destination image with the cluster
			break;
			iterations = i;
		}
		
		//If here, at least another iterations should be done to get the convergence: reset the clusters
		//but keep the new centroid calculated: it will be useful in the next iteration to find new pixels
		//to push into the cluster (with the new calculated centroid)
		for(int k=0; k<clusters.size(); k++) {
			clusters.at(k).clear();
		}
	
	}
	
	//Iterations ended: all the possible cluster has been found (either for the convergence or iterations finished).
	//Let's build the final image.
	
	//x and y coordinate to take track of the current analyzed pixel
	int x,y;
	
	//for all the clusters
	for(int k=0; k<clusters.size(); k++) {
		//for the size of the specific K-Cluster (to analyze all of his pixels)
		for(int i=0; i<clusters.at(k).pixels.size(); i++) {
			//fetch the x and y coordinate from the pixels point array
			x = clusters.at(k).pixels.at(i).x;
			y = clusters.at(k).pixels.at(i).y;
			
			//the destination at x,y is equal to the centroid of the cluster: with this assignment, we're going
			//to assign to a specific cluster pixel a fixed value to show into the image the cluster segmentation
			dest_image.at<Vec3b>(x,y) = clusters.at(k).centroid;
		}
	}
	
	//cout the iterations
	cout << "Iterations done: " << iterations << endl;
	
	//apply a median blur to remove the black pixel remained into the image
	medianBlur(dest_image, dest_image, 3);
	
}























