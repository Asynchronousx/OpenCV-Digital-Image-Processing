#include <opencv2/opencv.hpp>
#include <iostream>
#include <cstdlib>

//Note: for K-Means algorithm in-depth theory, refer to the K-Mean algorithm.
//Meng-Hee-Heng algorithm: a better approach to K-Means.
//From the theory, the K-Means is an algorithm that, inputted with K possible cluster into the image
//choosen by the user, will find all this possible clusters and assign each pixel to every K cluster
//found into the image, by finding the convergenge of all the mean of the centroid of the clusters 
//or when the iterations are finished. 
//The big disadvantage of the K-Mean algorithm is that the user should choose the number of clusters. 
//This could be good in the case of very few partition visibile through the image, but when the image
//itself starts to be more complicate, things get hard.
//The Meng-Hee-Heng K-Means algorithm, aims to that: finding all the possible cluster into the image, 
//choosing by his steps who should be a cluster and who not.
//The process of the Meng-Hee-Heng algorithm could be approssimated to the following steps:
//1)We starts with only two cluster, X and Y, formed each by a single pixel (and this pixel is the centroid
//  of the new cluster formed) choosen by a criteria: the pixel that forms X and the pixel that forms Y
//  should be the most distant pixel (chromatically) from each other in the BGR space.
//2)We group all the pixel of the image into the new clusters formed, assigning a pixel to a cluster
//  by analyzing his distance from the centroid. If it's similar, then assign the pixel to that cluster.
//3)We now have all the pixels of the image belonging to some specific cluster. We now began the process
//  of finding new clusters: for each cluster, we take the most distant pixel from the centroid of the cluster
//  itself, let's call it Z, and his distance value: let's call that D, by analyzing all the distance of the pixel
//  from the centroid.
//4) Now, let q be the distance between each pair of cluster calculated by the distance from their centroid.
//   (So q is the average between each pair of distance)
//   If the distance D of the most far pixel into a cluster is > q/2, then the pixel Z forms a new cluster.
//5) If we find a new cluster we start again from step 2. The algorithm terminates when each distance D is
//   lower than q/2: and that have sense, because when D < q/2 we know for sure that the distance from a pixel
//   from his centroid is not so relevant to make a new cluster.

using namespace std;
using namespace cv;

//Class to take track of the total sum of the region: this will be useful to track the centroid value.
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


//Class that represent a cluster into our image. Into this cluster we store sensitive information such as the
//centroid of the cluster, the pixels vector that contains all the pixels belonging to a vector, the total sum
//of the region and the pixel number into that Cluster. We also store useful method to track the mean of the cluster.
class Cluster {
	
	public:
	int pixel_num;
	Average3b region_sum;
	Vec3b centroid;
	vector<Point> pixels;
	
	Cluster() {;}
	Cluster(Vec3b centroid, int x, int y) {
		//initially all is set to zero. It will be job of the function add pixel
		//and compute centroid to init the cluster.
		this->pixel_num = 0;
		this->centroid = Vec3b(0,0,0);
		this->region_sum = Average3b(0,0,0);
		add_pixel(centroid, x, y);
	}

	//Computing the new centroid: it is easily calculated assigning to each centroid channel the
	//division of the appropriate channel by the pixel num in that region.
	void compute_centroid() {
		this->centroid[0] = region_sum.b/pixel_num;
		this->centroid[1] = region_sum.g/pixel_num;
		this->centroid[2] = region_sum.r/pixel_num;
	}
	
	//Adding a pixel to the cluster, means that a new sample is entering the group. 
	//Due to that, we update the region sum, push the pixel	into the cluster's pixel vector and
	//increase the pixel number.
	void add_pixel(Vec3b pixel, int x, int y) {
		//Update the region sum
		this->region_sum.b += pixel[0];
		this->region_sum.g += pixel[1];
		this->region_sum.r += pixel[2];
		
		//push the new pixel into the vector
		this->pixels.push_back(Point(x,y));
		
		//increase the pixel number
		this->pixel_num++;
		
		//compute centroid with the new pixel added
		compute_centroid();
		
	}
	
	//function that takes a pixel in input, and erase it's content from the cluster.
	void erase(Vec3b pixel, int index) {
		//Decrease the total region sum value by subtracting the value of the pixel
		//that's is gonna be erased
		this->region_sum.b -= pixel[0];
		this->region_sum.g -= pixel[1];
		this->region_sum.r -= pixel[2];
		
		//decrease the pixel number
		this->pixel_num -= 1;
		
		//compute the new centroid
		compute_centroid();
		
		//erase the pixel from the pixels vector
		this->pixels.erase(this->pixels.begin() + index);
	
	}
	
	//clearing the cluster 
	void clear() {
		//pixel number is resetted to one: the centroid
		this->pixel_num = 1;
		
		//clearing the pixels vector
		this->pixels.clear();
		
		//assigning to the region sum the new centroid 
		this->region_sum.b = centroid[0];
		this->region_sum.g = centroid[1];
		this->region_sum.r = centroid[2];
		
		 //the only thing we do not reset is the centroid because it was already modified from another
 		//call into the k-means algorithm.
	}
	
};

void MengHeeHeng(Mat, Mat&, int iterations);
vector<Cluster> furthest_point(Mat);

int main(int argc, char** argv) {

	if(argc!=3) {
		cerr << "Usage: ./<programname> <image.format> <iterations>" << endl;
		exit(EXIT_FAILURE);	
	}
	
	Mat raw_image = imread(argv[1], IMREAD_COLOR);
	
	if(raw_image.empty()) {
		cerr << "Image format not valid." << endl;
		exit(EXIT_FAILURE);
	}
	
	//Declaring dest image same as the raw image size and type
	Mat dest_image(raw_image.size(), raw_image.type(), Scalar(0));
	
	//fetching the threshold passed in input
	int iterations = atoi(argv[2]);
	
	//Calling the function
	MengHeeHeng(raw_image, dest_image, iterations);

	//Showing the result
	
	imshow("Original", raw_image);
	imshow("Meng-Hee Heng", dest_image);
	waitKey(0);

}

void MengHeeHeng(Mat raw_image, Mat& dest_image, int iterations) {

	//First step of the algorithm: find the most distant pixel into the image from each other
	//We use the furthest_point function to return the vector of initial cluster composed by X and Y.
	vector<Cluster> clusters = furthest_point(raw_image);
	
	//We're going to enter to the block that will calculate the cluster into the image.
	//We set a limite amount of iterations to make sure the algorithm can converge to a solution at some point.
	for(int it=0; it<iterations; it++) {
	
		cout << "Iteration: " << it << endl;

		//Index is useful to know which cluster should add a pixel in a certain k iteration when dividing the pixels 
		//by cluster.
		//Since we going to create a new cluster with the furthest pixel, we must erase that pixel from the
		//pixels cluster vector. We hold the erase_index to erase it when making the comparison D > q/2
		double distance;
		int index;
	
		//boolean variable to check if we got the convergence: no new cluster formed
		bool changed = false;
		
		//Second step of the algorithm: given our cluster, we're going to insert each pixel of the image into
		//an appropriate cluster, based on how much they dist from the centroid.
		//for each pixel into the image
		for(int i=0; i<raw_image.rows; i++) {
			for(int j=0; j<raw_image.cols; j++) {		
				//Minimum distance is set at each iteration of a new pixel: this will help us choose
				//to which cluster a pixel would belong better, compare at each k-th + 1 iteration
				//the new distance found with the already-found min distance. 
				double min_dist = DBL_MAX;	
				
				//let's compare the distance of every pixel of the image from the 
				//centroid of the found cluster	
				for(int k=0; k<clusters.size(); k++) {
					//check for the distance that should be lower than an already set-up
					//min distance: if so, memorize the index and set the new min distance as the
					//distance just found.
					distance = norm(raw_image.at<Vec3b>(i,j), clusters.at(k).centroid);
					if(distance < min_dist) {
						index = k;
						min_dist = distance;
					}			
				}
				
				//When all the iterations are done, we know for sure to which cluster the pixel belong.
				//so we use the index that memorized the optimal cluster to push the pixel into that cluster
				//pixels vector.
				clusters.at(index).add_pixel(raw_image.at<Vec3b>(i,j), i, j);
				
			}	
		}

		
		//Third step: we got all our cluster set up with all the pixel into the image.
		//We're going to take the most distant pixel Z of all the most distant pixel in
		//all the clusters. Let's call his distance D.
		
		//Max distance D is the distance of the furthest point into a cluster, the new potential cluster
		//centroid, his potential coordinate pot_x and pot_y.
		//We also use a variable distance to take the comparison with the max distance at each iteration
		//and temporary ind_x and ind_y variable to analyze each x,y of the current cluster pixel analyzed.
		//Max Distance is set to zero to let the first iteration discover the distance.
		double max_distance = 0;
		int ind_x, ind_y, pot_x, pot_y, erase_index, k_index;
		Vec3b potential_centroid;
		
		for(int k=0; k<clusters.size(); k++) {
			//for each cluster pixel vector, we find the furthest pixel from the centroid of that cluster:
			//found that, we add the pixel to a "potential cluster" vector, to analyze it later with 
			//the computed average distance q between each cluster
			for(int i=0; i<clusters.at(k).pixels.size(); i++) {
				ind_x = clusters.at(k).pixels.at(i).x;
				ind_y = clusters.at(k).pixels.at(i).y;
				distance = norm(clusters.at(k).centroid, raw_image.at<Vec3b>(ind_x,ind_y));
				
				//check the distance: if it > than the actuall max distance, the distance is
				//the new max distance
				if(distance > max_distance) {
					//take track of the potential centroid that will form a new cluster, his coordinate,
					//the new distance found and the index at which the pixel should be erase, alongside
					//with the k index to specify in which cluster
					potential_centroid = raw_image.at<Vec3b>(ind_x, ind_y);
					max_distance = distance;
					pot_x = ind_x;
					pot_y = ind_y;
					erase_index = i;	
					k_index = k;	
				}
			}					
		}
		
		
		//Fourth step: find the average q distance between each pair of cluster into the vector. 
		//For each cluster, we're going to find his distance between all the other clusters in the vector.
		//For each cluster
		//We initialize the first for to 0 because we want to start at index 0 of our vector.
		//For the second for cycle, we initialize the second scanning on the cluster vector to x+1,
		//because we do not want to find the same distance from each pair of cluster again.
		//Example: We got 3 cluster: 
		// 0 - 1
		//  \2/
		//We want to find the medium distance between each pair: starting at zero, we compare 0 to 0+1 (1) 
		//skipping the comparison with himself. The second for iterate and we compare 0 to 2 and found the second
		//distance.
		//When the first for iterate and x=1, then y=x+1 -> y=1+1 -> 2.
		//We then confront the first cluster with the second, becase we still missed that distance calculation.
		//In the third for, x=2 and y=2+1->3, and we stop because it exceed the clusters size.
		//With that approach, we skipped not only the comparison with a cluster with himself, but also the 
		//already analyzed cluster distance.
		
		//q distance between each couple of cluster into the vector.	
		double q;
		
		for(int x=0; x<clusters.size(); x++) {
			//and for each other cluster in the vector 
			for(int y=x+1; y<clusters.size(); y++) {
				//if the same index, we dont mind calculate the distance from a cluster to himself
				if(x==y) continue;
				
				//adding to q the distance of each cluster to another
				q += norm(clusters.at(x).centroid, clusters.at(y).centroid);
								
			}		
		}
		
		
		//We're dividing q for the clusters size to have the average distance from each pair of cluster
		//found into the image.
		q /= clusters.size();
		
		
		//Fifth step: checking if D > Q/2
		//once here, we know the furthest distance of all the distance into the cluster: 
		//we make the comparison with the average distance q/2, and if D > q/2 we add the potential 
		//cluster to the clusters vector and remove the pixel from the k-th found clusters pixels vector.
		if(max_distance > q/2) {
			
				//new cluster are being formed, changed is true
				changed = true;
				
				//erase the pixel from the cluster: this also update the centroid.
				clusters.at(k_index).erase(potential_centroid, erase_index);
				
				//add the new cluster to the potential cluster vector: 
				//New cluster found, which have the furthest D of all the Distance of the most distant pixel of 
				//each cluster. THis cluster is then, the furthest of the furthest pixel of each cluster.
				clusters.push_back(Cluster(potential_centroid, pot_x, pot_y));
				
		}
		
				
		//check if the iterations are depleted or the convergency has been met
		if(!changed || it+1 == iterations) {
			break;
		}
		
		//if we still need to iterate, let's clear all the cluster from their value lefting only the
		//computed centroid: this will be useful to track the new belonging of the pixel in the next 
		//iteration.
		for(int k=0; k<clusters.size(); k++) {
			clusters.at(k).clear();		
		}
		

	}

	/*
	//Alterntive: random colors
	RNG rng (12345);
	int x, y;
	Vec3b random;
	
	for(int k=0; k<clusters.size(); k++) {
		random = Vec3b(rng.uniform(0,255), rng.uniform(0,255), rng.uniform(0,255));
		clusters.at(k).centroid = random;
	}
	*/
	
	//Build the image with the found cluster
	cout << "Number of cluster found: " << clusters.size() << endl;
	
	int x, y;
	for(int i=0; i<clusters.size(); i++) {
		for(int j=0; j<clusters.at(i).pixels.size(); j++) {
			x = clusters.at(i).pixels.at(j).x;
			y = clusters.at(i).pixels.at(j).y;
			dest_image.at<Vec3b>(x,y) = clusters.at(i).centroid;
		}
	}
	
	//omogenize the image using a medianblur
	medianBlur(dest_image, dest_image, 3);

}

vector<Cluster> furthest_point(Mat raw_image) {
	
	//Helper function that aim to find the most distant pixel between them into the image.
	//Since analyzing all the pixel to every other pixel is really computational-heavy, we use a sthocastic
	//approssimation of the max distance based on a simple concept:
	//We take the mean of all the pixel into the image, and then confronting each pixel with that mean: 
	//the pixel that dist most from that mean, is for sure one of the most distant pixel from the average pixel:
	//from there, we find the most distant pixel from the most distant pixel from the mean, and we got our approssimation.
	
	//Declaring an Average3b object to store the sum of all the pixel into the image
	Average3b image_sum(0,0,0);
	int total_pixel = raw_image.rows * raw_image.cols;
	
	//Declaring distance and max distance to find out the maximum distant pixel into the image:
	//Max distance is initially set to 0 to let the first iteration discover the initial distance
	double max_distance = 0, distance;
	
	//Storing the most distant pixel into two Vec3b element
	Point p1, p2;
	Vec3b mean;
	
	//The vector of cluster to return
	vector<Cluster> initial_cluster;
	
	//Finding the mean of all the pixel into the image
	for(int i=0; i<raw_image.rows; i++) {
		for(int j=0; j<raw_image.cols; j++) {
			image_sum.b += raw_image.at<Vec3b>(i,j)[0];
			image_sum.g += raw_image.at<Vec3b>(i,j)[1];
			image_sum.r += raw_image.at<Vec3b>(i,j)[2];
		}	
	}
	
	//getting the average pixel dividing all the channels by the total pixel number, and assigning that to 
	//the average pixel
	mean[0] = image_sum.b/total_pixel;
	mean[1] = image_sum.g/total_pixel;
	mean[2] = image_sum.r/total_pixel;
	
	//Find the most distant pixel from the average mean
	for(int i=0; i<raw_image.rows; i++) {
		for(int j=0; j<raw_image.cols; j++) {
			//calculating the norm between the current pixel and the mean
			distance = norm(raw_image.at<Vec3b>(i,j), mean);
			
			//if distance found is major than the actual max distance
			if(distance > max_distance) {
				//update the point and the distance
				p1.x = i;
				p1.y = j;
				max_distance = distance;
			}
			
		}	
	}
	
	//Now find the most distant pixel from the last one with the same approach
	max_distance = 0;
	for(int i=0; i<raw_image.rows; i++) {
		for(int j=0; j<raw_image.cols; j++) {
			distance = norm(raw_image.at<Vec3b>(i,j), raw_image.at<Vec3b>(p1.x, p1.y));
			if(distance > max_distance) {
				max_distance = distance;
				p2.x = i;
				p2.y = j;
			} 
		}
	}
	
	//building the vector of initial X,Y cluster with the point found
	initial_cluster.push_back(Cluster(raw_image.at<Vec3b>(p1.x, p1.y), p1.x, p1.y));
	initial_cluster.push_back(Cluster(raw_image.at<Vec3b>(p2.x, p2.y), p2.x, p2.y));

	//return the vector
	return initial_cluster;
	
}


