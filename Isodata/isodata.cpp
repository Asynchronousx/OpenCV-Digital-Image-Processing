#include <opencv2/opencv.hpp>
#include <iostream>

//Isodata clustering: another approach of K-Means
//What is isodata clustering?
//The isodata clustering is a different approach of Clustering based on the logic of the K-Mean.
//In the K-Mean, is the user who choose the number K of cluster, and those K cluster are random choosen
//by selecting random pixel into the image.
//In the Meng-Hee Heng approach, we start from two cluster, formed by the most-furthest pixel into the image.
//We then expand the cluster number from those two cluster to find the optimal number of cluster into the image.
//In the Isodata clustering, the approach is different: We start with an huge number of cluster, and through
//the iteration we converge to the optimal number of cluster. 
//How can we achieve that? 
//The number of initial cluster is limited upperly by the number of pixel into the image (MxN). We can choose
//to pass a number of initial cluster or following a pattern (example: every 50 pixel select a pixel to form a cluster).
//We can define those steps to identify the isodata algorithm:
//Note: the user should pass TWO threshold in input: one for the variance, one for the average
//1) Select an huge amount of initial clusters. We can either do this by passing the number of cluster we want initially, 
//   or following a pattern in which we take a cluster every N pixels. 
//2) We then iterate through the image, and assign to each cluster formed, each possible pixel into the image that have
//   the shortest distance from the centroid.
//3) Once we got our cluster formed, we then proceed to choose if a cluster should be splitted or merged comparing his
//   variance with an user defined threshold. If the variance is higher than the threshold, we proceed to split the cluster
// 	 into two subcluster. Note: we can check the variance by simply calculating the distance of the three channel of all the pixels
//   from the cluster centroid. To assign a pixel to a subcluster, we analyze NOT the variance threshold but an average threshold
//   different from the variance one. 
//4) Once we split all the cluster into N subcluster, we then proceed to compute all the distance between each pair of cluster,
//   and merge who's distance is lesser than an user defined average threshold.
//5) We recompute the means for each cluster, because merging two or more cluster between them altered their centroid.
//6) Repeat step 2-5 until The algorithm should iterate until: the average distance between each cluster falls down the user 
//   define threshold,  the changing of a centroid from his last one is not significative (down the user threshold) or 
//   the max iterations are reached.

using namespace std;
using namespace cv;

Mat public_image;

//Class that represent the total sum of chromatic intensity of a given pixel into a given cluster
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
	
	//funciton that sum up the value of a certain pixel to the current sum
	void add_values(Vec3b pixel) {
		this->b += pixel[0];
		this->g += pixel[1];
		this->r += pixel[2];
	}
	
	void add_values(Average3b p1, Average3b p2) {
		this->b += p1.b + p2.b;
		this->g += p1.g + p2.g;
		this->r += p1.r + p2.r;	
	}
	
};

//Class that represent a cluster into our image. Clusters got useful information such as the number of the 
//pixel contained into it, the total sum of the chromatic intensity of the cluster, the pixels from which the cluster
//is composed, the centroid and the variance of a given cluster to check if a cluster should be split or not.
// We also got some useful methods to take track of the computing of the centroid, adding a pixel etc.
class Cluster {

	public:
	int pixel_num;
	Average3b average_sum;
	double variance;
	Vec3b centroid;
	vector<Point> pixels;
	
	Cluster() {
		//Even if calling the default constructor we need to init the variables
		//because they could be used.
		this->pixel_num = 1;
		this->average_sum = Average3b(0,0,0);
		this->centroid = Vec3b(0,0,0);
		this->variance = 0;
	}
	
	Cluster(Vec3b centroid, int x, int y) {
		this->pixel_num = 0;
		this->average_sum = Average3b(0,0,0);
		this->centroid = Vec3b(0,0,0);
		this->variance = 0;
		add_pixel(centroid, x, y);
	}
	
	//Function that compute the centroid of the cluster
	void compute_centroid() {
		//given by: the toal sum of each channel / the pixel number
		this->centroid[0] = average_sum.b/pixel_num;
		this->centroid[1] = average_sum.g/pixel_num;
		this->centroid[2] = average_sum.r/pixel_num;
	}
	
	//Function that compute the variance of a cluster
	void compute_variance() {
	
		//The variance of a cluster is meant to be the sum of all the distance of the pixel of the cluster
		// from the centroid, then divided for the number of the pixels.
		int x, y;
		for(int i=0; i<pixels.size(); i++) {
			x = this->pixels.at(i).x;
			y = this->pixels.at(i).y;
			
			this->variance += norm(public_image.at<Vec3b>(x,y), this->centroid);

		}
	
		this->variance /= pixel_num;		
			

	}
	
	//Function that takes care a pixel adding into the cluster
	void add_pixel(Vec3b pixel, int x, int y) {
		//Add the pixel value to the region sum
		this->average_sum.add_values(pixel);
				
		//increase the pixel number
		this->pixel_num++;
	
		//add the pixel coordinate to the pixels vector
		this->pixels.push_back(Point(x,y));
		
		//compute the new centroid
		compute_centroid();
	
	}
	
	//clearing the cluster 
	void clear() {
		//pixel number is resetted to one: the centroid
		this->pixel_num = 1;
		
		//clearing the pixels vector
		this->pixels.clear();
		
		//assigning to the region sum the new centroid 
		this->average_sum.b = this->centroid[0];
		this->average_sum.g = this->centroid[1];
		this->average_sum.r = this->centroid[2];
	
		//clear the variance		
		this->variance = 0;
		
		 //the only thing we do not reset is the centroid because it was already modified from another
 		//call into the k-means algorithm.
	}

};


void Isodata(Mat, Mat&, int, int, int);

int main(int argc, char** argv) {

	if(argc!=5) { 
		cerr << "Usage: ./<programname> <image.format> <variance threshold> <average threshold> <Max Iterations>" << endl;
		exit(EXIT_FAILURE); 	
	}
	
	Mat raw_image = imread(argv[1], IMREAD_COLOR);
	
	if(raw_image.empty()) {
		cerr << "Image format not valid" << endl;
		exit(EXIT_FAILURE);	
	}
	
	//Assigning raw image to public image
	public_image = raw_image;
	
	//Creating dest image filled with 0
	Mat dest_image(raw_image.size(), raw_image.type(), Scalar(0));
	
	//Assigning max iteration and threshold
	int vt = atoi(argv[2]);
	int at = atoi(argv[3]);
	int iterations = atoi(argv[4]);

	//Calling isodata
	Isodata(raw_image, dest_image, vt, at, iterations);
	
	imshow("Original", raw_image);
	imshow("Isodata", dest_image);
	imwrite("isodata.png", dest_image);
	waitKey(0);
	
}

void Isodata(Mat raw_image, Mat& dest_image, int vt, int at, int iterations) {

	//Vector that contains all the cluster in the image
	vector<Cluster> clusters;
	vector<Cluster> tmp_clusters;
	
	//Useful variables to take track of the index of cluster at which a point should be inserted,
	//the distance between a pixel and the centroid
	int index, x, y;
	double distance;
	
	//Isodata clustering, step 1: pick an huge amount of cluster.
	//We choose to follow a pattern and choose a cluster each 40 pixel iterated
	for(int i=0; i<raw_image.rows; i+=80) {
		for(int j=0; j<raw_image.cols; j+=80) {
			//Pushing a new cluster formed by the current pattern pixel
			clusters.push_back(Cluster(raw_image.at<Vec3b>(i,j), i, j));		
		}	
	}
		
	cout << "Initial Clusters: " << clusters.size() << endl;
	//Step 1/2/3/4: we're now going to split and merge all the clusters until a convergence is met.
	//The first requirement are the max number of iteration, so let's iterate from 0 to them
	for(int it=0; it<iterations; it++) {
		
		cout << "Iterations: " << it+1 << endl;

		//Step 1.a: we are going to assign every pixel of the image into one of the cluster, finding 
		//the cluster with the centroid that dist less from the analyzed pixel. 
		//At every iteration we calculate again the pixel to assign to each regions
		for(int i=0; i<raw_image.rows; i++) {
			for(int j=0; j<raw_image.cols; j++) {
				//Initializing a min dist each time we analyze a pixel: with that minimum distance,
				//we check every distance from the pixel to the centroid of each cluster, to choose
				//whose should add the pixel into his cluster.
				double min_dist = DBL_MAX;
				
				//for each cluster into the vector
				for(int k=0; k<clusters.size(); k++) {
					//calculating the distance with the norm of the centroid and the current pixel
					distance = norm(raw_image.at<Vec3b>(i,j), clusters.at(k).centroid);
					
					//if the distance is lower than the minimum distance found until now
					if(distance < min_dist) {
						//update the min dist and the index found
						index = k;
						min_dist = distance;			
					}
				}
				
				//once here, we know for sure to which cluster a pixel should belong. 
				//we're using the index found to access the cluster at (index) position,
				//and adding that pixel into the cluster
				clusters.at(index).add_pixel(raw_image.at<Vec3b>(i,j), i, j);	
					
			} 
		}

		//Once we got all the clusters set up, we then compute the variance for each cluster.
		for(int i=0; i<clusters.size(); i++) {
			clusters.at(i).compute_variance();	
		}
	
		//Step 2:		
		//We re going to split the cluster into subcluster if their variance is higher than a user-inputted
		//threshold.
		bool splitted_once = false;	
		
		for(int k=0; k<clusters.size(); k++) {
			//if the variance is higher than the threshold, we're going to split the cluster into two
			//subcluster. We can check this by calculating the norm between the centroid and the variance:
			//if the result is higher than t, that means in that cluster we got pixels that should not belong
			//there.
			if(clusters.at(k).variance > vt) {
				//someone has been splitted
				splitted_once = true;
			
				//initializing two empty cluster that represent subcluster1 and subcluster2
				Cluster split1, split2;
				
				//now, since we know there are pixel that don't belong to the cluster, we make a control:
				//if the norm between the intensity of a pixel of that cluster, and the centroid is higher than t,
				//that means, for example, that the pixel belong to the second subcluster, and if lower to the first.
				//let's iterate all the pixel of the current cluster and verify that
				for(int i=0; i<clusters.at(k).pixels.size(); i++) {
					//fetch the coordinate of that pixel
					x = clusters.at(k).pixels.at(i).x;
					y = clusters.at(k).pixels.at(i).y;
					
					//check the norm
					if(norm(raw_image.at<Vec3b>(x,y), clusters.at(k).centroid) < at) {
						//pixel belong to the first cluster: also, the first pixel inserted will be the centroid.
						split1.add_pixel(raw_image.at<Vec3b>(x,y), x, y);
					} 
					else {
						//pixel belong to the second cluster: also the first pixel inserted will be the centroid.
						split2.add_pixel(raw_image.at<Vec3b>(x,y), x, y);
					}
				}
				
				//erasing the pixel at the k-th index and inserting the two splitted cluster into the tmp vector,
				//and decrease the index
				clusters.erase(clusters.begin() + k);
				tmp_clusters.push_back(split1);
				tmp_clusters.push_back(split2);
				k--;
				
			}
		}
		
		

		//tmp_clusters contains all the splitted cluster from a single one. Let's append it to the cluster
		//vector
		clusters.insert(clusters.end(), tmp_clusters.begin(), tmp_clusters.end());
		
		//clean the vector for the merging phase
		tmp_clusters.clear();
		
		//Step 3: after the splitting of the clusters (and also the new clusters generated), we need to check if some 
		//of the clusters needs to be merged. We compute, for each pair of cluster, their distance (based on the centroid)
		//and if it is lesser than the threshold, we proceed to merge them into one cluster.
		//For each cluster into the vector.
		//We init this bool value to check if at least 1 merge happened
		bool merged_once = false;
		
		for(int k=0; k<clusters.size(); k++) {
			
			//min dist is re-computed at each new iteration
			double min_dist = DBL_MAX;
			
			//For each cluster into the vector starting from i+1 to the cluster size, increasing by 1: this because,
			//at every iteration we're not considering the already analzyed clusters and the cluster with himself,
			//comparing only the clusters that need to be compared. To let the merging be the best possible,
			//we need to compare each cluster to the current one, and then choose who's distance is lesser.
			for(int j=k+1; j<clusters.size(); j++) {
				//getting the distance
				distance = norm(clusters.at(k).centroid, clusters.at(j).centroid);
				
				//if the distance is lesser than the threshold, that means that we got a merging 
				//candidate
				if(distance < at) {
				
					//merged cluster is the result of the merging of two cluster
					Cluster merged_cluster;
				
					//the pixel num of the merge cluster is equal to the sum of both
					merged_cluster.pixel_num = (clusters.at(k).pixel_num + clusters.at(j).pixel_num);
					
					//we're going to append both of the pixel vector of the two cluster to the merged one
					merged_cluster.pixels = clusters.at(k).pixels;			
					merged_cluster.pixels.insert(merged_cluster.pixels.end(), clusters.at(j).pixels.begin(), clusters.at(j).pixels.end());		
					
					//we're going to sum into the average sum of the average both of the average sum of the cluster
					merged_cluster.average_sum.add_values(clusters.at(k).average_sum, clusters.at(j).average_sum);

					//and we're going to do the same for the variance
					merged_cluster.variance = (clusters.at(k).variance + clusters.at(j).variance);
							
					//after the region sum merging, compute the centroid of the new merged region		
					merged_cluster.compute_centroid();
									
					//deleting the two merged cluster, using their relative index
					clusters.erase(clusters.begin() + k);
					clusters.erase(clusters.begin() + j);
					
					//Since we erased 2 element, subtract 2 to the current index
					k -= 2;
					
					//push the new cluster into the tmp vector
					tmp_clusters.push_back(merged_cluster);
					
					//we merged at least once
					merged_once = true;
					
					//break the cycle, we don't want to merge the k-th cluster with another one
					break;
								
				}			
			}
			
		}

		//append the merged clusters to the clusters vector: each time we merged two cluster, we deleted two cluster
		//into the original vector. With that append, we're going to insert into the cluster the new merged one
		clusters.insert(clusters.end(), tmp_clusters.begin(), tmp_clusters.end());
	
		//check for any exit condition
		if(!splitted_once && !merged_once || it+1 == iterations) {
			break;
		}
		
		//once here, we know that a next iteration is gonna happpen: clear each cluster, we're starting from the beginning agains
		for(auto &c: clusters) c.clear();	
	
	}
	
	// Uncomment For random colors
	/*
	RNG rng(12345);
	Vec3b r;
	for(auto &c: clusters) {
		r = Vec3b(rng.uniform(0,255), rng.uniform(0,255), rng.uniform(0,255));
		c.centroid = r;	
	}
	*/
	
	//Once here, the convergency has been met or the maximum number of iteration are depleted. We then proceed, for the 
	//resultant cluster, to process the output image.
	cout << "Cluster Found: " << clusters.size() << endl;
	for(int i=0; i<clusters.size(); i++) {
		for(int j=0; j<clusters.at(i).pixels.size(); j++) {
			x = clusters.at(i).pixels.at(j).x;
			y = clusters.at(i).pixels.at(j).y;
			dest_image.at<Vec3b>(x,y) = clusters.at(i).centroid;
		}
	}	
	
	//applying a median blur
	medianBlur(dest_image, dest_image, 3);
	
}





