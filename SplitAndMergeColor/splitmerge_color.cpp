#include <opencv2/opencv.hpp>
#include <iostream>
#include <math.h>
#include <cstdlib>

//Split and Merge: 
//The split and merge algorithm, like the region growing, is an algorithm that aim to the segmentation of the image.
//Differently from the region growing algorithm, that follow a bottom-up approach, the split and merge algorithm is a 
//top-down methodology.
//The algorithm itself is explanatory: At each pass, if the current analyzed region of the image does not respect a critery
//of omogenity P(Ri), the region is divided in four subregion, and each of those region are processed aswell, getting splitted
//if they do not respect the P(Ri) and so on, until the P(Ri) is respected or the minimum split size has been met.
//Each region splitted is pushed into a vector of regions splitted by this criteria.
//When all the splitting is done, we must analyze this vector of splitted region to merge all the similar region.
//How? we first calculate all the adjacence from a region to all his neighbor (N,S,E,W, excluding the diagonal neighbor)
//And then proceed to merge all the adjacent region with a similar P(Ri), after the adjacence calculation.
//We can calculate the P(Ri) with the variance of the entire region: given the average intensity of a given region, 
//Calculate the variance by the formul pow(average-pixel(x,y),2) or with the deviation standard given by the 
//euclidean distance of two pixel: norm(average_pixel - pixel(x,y)).
//Note: pixel(x,y) represent each pixel of the region that we would like to know the average.

using namespace std;
using namespace cv;

//Public image used at the creation of a new region or in the process of the building of an image to avoid passing
//heavy matrix objects through function.
Mat public_image;

//Class that simply represent the average pixel value utilizied into a specific region to take track of the total sum
//of the element of that region. We do not memorize the pixel number into this class because we already have it into the
//region attributes. We just accumulate the values into those attributes and dividing them for the number of the pixel
//to obtain the average value.
class Average3b {
	public:
	int B;
	int G;
	int R;
	
	Average3b() {;}
	Average3b(int B, int G, int R) {
		this->B = B;
		this->G = G;
		this->R = R;
	}	

};


//Class that represent the attributes of a region, such as his area covered (start col/row -> end col/row), 
//the number of pixel of that region and the average sum/pixel to identify the variance of the region itself.
//We also mantain a list of adjacence created by a set.
class Region {

	public:
	Rect area;
	int pixel_num;
	double variance;
	Average3b average;
	Vec3b average_pixel;
	
	//The adjacancy set contains all the adjacent region: we declare it as a set of pointer to region,
	//because when we'll merge the region by their average, we'll need to access to the exact structure 
	//to make the change permanent.
	set<Region*> adjacency;
	
	Region() {;}
	
	Region(int start_row, int start_col, int end_row, int end_col) {
		//We intend the area as the two corner of the image: Top Left and Bottom right:
		//Note: Rect object is built by inverting rows and cols: x is the column, y is the row.
		//Same as the point structure.
		this->area = Rect(start_col, start_row, end_col, end_row);
		//Example: imagine that the area passed in input is: (sr,er) (sc,ec) -> (210, 219) (220,229).
		//There are 10 rows and 10 cols each. Since the last index is INCLUSE, we got that 219 - 210 is 9, and 229 - 220 is also 9.
		//We then must add 1 to make the sum inclusive aswell, and multiply the result: i.e (9+1) * (9+1) -> 100 element.
		this->pixel_num = ((area.height - area.y) + 1) * ((area.width - area.x) + 1);
		this->average = Average3b(0,0,0);
		compute_average();
		compute_variance();	
	}
	
	//public function to check if a point is inside a specified rectangle.
	//The logic behind it is simple: to check if a point belong to a rectangle, we 
	//check if the coordinate (x,y), intended as (column, row), are contained inside the rectangle bound.
	//If only ONE of the coordinate x or y are lesser than the coordinate x or y of the rectangle, or
	//major than the end_x, end_y that means that the point does not belong to that region. We return
	//false if it is not contained, true if it is.
	bool contains(Point p) {
		if(p.x < area.x || p.x > area.width || p.y < area.y || p.y > area.height) {
			return false;
		} else {
			return true;
		}	
	}
	
	private:
	//Computing the average of the region through his 3-channel pixel
	void compute_average() {
		Vec3b current_pixel;
		for(int i=area.y; i<=area.height; i++) {
			for(int j=area.x; j<=area.width; j++) {
				current_pixel = public_image.at<Vec3b>(i,j);
				this->average.B += current_pixel[0];
				this->average.G += current_pixel[1];
				this->average.R += current_pixel[2];
			}
		}
		
		this->average_pixel[0] = round((this->average.B)/this->pixel_num);
		this->average_pixel[1] = round((this->average.G)/this->pixel_num);
		this->average_pixel[2] = round((this->average.R)/this->pixel_num);
		
	}
	
	//Computing the average of a region by the sum of the variance between the norm of the average pixel
	//and the current pixel at (x,y) analyzed.
	void compute_variance() {
		for(int i=area.y; i<=area.height; i++) {
			for(int j=area.x; j<=area.width; j++) {
				this->variance += norm(average_pixel, public_image.at<Vec3b>(i,j));	
			}
		}

		this->variance /= pixel_num;
	
	}

};


void SplitAndMerge(Mat, Mat&);
void split(int, int, int, int);
void calculate_adjacence();
void merge();
void build_image_rect(Mat&);
void build_image_norect(Mat&);

//Defining a vector of regions and the threshold
int t, pixel_limit;
vector<Region> regions;

int main(int argc, char** argv) {

	if(argc!=4) {
		cerr << "Usage: ./<programname> <imagename.format> <threshold> <minimum pixels limit>" << endl;
		exit(EXIT_FAILURE);
	}
	
	Mat raw_image = imread(argv[1], IMREAD_COLOR);
	t = atoi(argv[2]);
	pixel_limit = atoi(argv[3]);
	
	if(raw_image.empty()) {
		cerr << "Image format not valid." << endl;
		exit(EXIT_FAILURE);
	}
	
	//Resize the raw image: this because we want to operate with a squared image to let the segmentation by
	//splitting be optimal. We choose the max of (rows, cols) to mantain the ratio.
	int opt_size = max(raw_image.rows, raw_image.cols);
	
	//Resizing the image to the desired RowxHeight 
	resize(raw_image, raw_image, Size(opt_size, opt_size));
	
	//create the dest image filled with 0
	Mat dest_image(raw_image.size(), raw_image.type(), Scalar(0));
	
	//Apply gaussian blur to the image
	GaussianBlur(raw_image, raw_image, Size(5,5), 3, 3);
	
	//Display the raw processed image
	imshow("Original", raw_image);
	
	//assign 
	public_image = raw_image;
	
	//call function
	SplitAndMerge(raw_image, dest_image);
	
	//show
	imwrite("out.png", dest_image);
	waitKey(0);
	return 0;

}

void SplitAndMerge(Mat raw_image, Mat& dest_image) {
	
	//Split dest is our intermediate image: it will show the segmentation producted by the splitting of the image itself.
	Mat dest_split = dest_image.clone();
	Mat dest_split_rect = dest_image.clone();
	Mat dest_image_rect = dest_image.clone();
	
	//Fetching the coordinate of the first region: initially, it's the whole image.
	//We're subtracting 1 because we're starting from 0.
	int rows = raw_image.rows - 1;
	int cols = raw_image.cols - 1;
	
	//Calling the split function: we're going to repartition our image in various regions if a certain
	//Omogenity criteria P(Ri) is not respected. We can easily calculate this criteria with the variance of a given region.
	//If the variance does not respect the threshold, we must repartite our image in 4 sub images. 
	split(0, rows, 0, cols);
	
	//Build the splitted image and display it
	build_image_norect(dest_split);
	build_image_rect(dest_split_rect);
	imshow("Splitted Image No Rect", dest_split);
	imshow("Splitted Image Rect", dest_split_rect);
	waitKey(0);
	
	//calulcate the adjacence
	calculate_adjacence();
	
	//merge
	merge();	
	
	//build again to display in the main
	build_image_rect(dest_image_rect);
	build_image_norect(dest_image);
	
	//showing the images
	imshow("Split&Merge No Rect", dest_image);
	imshow("Split&Merge Rect", dest_image_rect);
	
}

void split(int start_row, int end_row, int start_col, int end_col) {

	//Initially, in the split function, we want to know if a certain property is verified: 
	//If the P(Ri) is true. 
	//To know that, we need to know the variance of a certain region, and if the variance is lesser thant the threshold
	//passed in input by the user, the P(Ri) results true.
	//To know the variance of a region, we must do, for every pixel the |average-pixel(x,y)| of the region, where x = 0..M and 
	//Y = 0..N.
	//Since we've incapsulated the computing of the variance of a certain region, given it's coordinate into the class region itself,
	//all we must do is initialize a region with his start row/col and end row/col. 
	//Then we proceed to compare the variance to the threshold.
	//If the variance is lesser than the threshold, this means that the region Ri satisfied the P(Ri), so we need to push that region into 
	//the vector of the regions already computed, that will be useful later for the adjacency and merge check.
	//initialize the region with his starting and ending area
	Region region(start_row, start_col, end_row, end_col);
	
	//Check the area size: we obtain the total size of the currently analyzed area by the difference of start-end, for both row and cols,
	//but, for each value, increasing by one the result: that's because we're starting from zero to (M/N) -1. 
	//I.e: current area is: start_row = 10, end row = 15 and start col = 10, end col = 15. 15-10 = 5 for both rows and col.
	//But we know that from 10 to 15, we got a 6x6 matrix. that's the reason we should add 1 to the result.
	//Creating the region we store this information in pixel_num already.
	
	//First check: is the region big enough to be splitted again? 
	//We try to "split" the actual pixel number of the region by 4: the quadrant of the resulting splitting.
	//If the pixel number is still higher than the pixel limit, that means we can split. Otherwise, splitting
	//could break this bound.
	if(region.pixel_num/4 <= pixel_limit) {
		//we can't do anything with that region: just push the region into the vector of processed region.
		regions.push_back(region);
		return;
	}
	
	//If here, that means that we can split the region in subregions. Before doing that, let's check if the variance of that 
	//regions is smaller than the threshold. if yes, we should put the region into the processed region vector because we don't
	//need to split that anymore.
	if(region.variance < t) {
		
		regions.push_back(region);
		return;
	}
	
	//If here, the region was bigger than the pixel limit threshold but the variance was higher than the threshold.
	//we must split this region in another 4 subregion.
	//We can do that by considering the image as quadrants organized in this way:
	// |Q1|Q2|
	// |Q3|Q4|
	//Let's explain how we can divide the region in that way:
	//Note: r means that the value should be rounded.
	//Q1: R: start row to r((start row + end row)/2) - C: start col to r((start col + end col)/2)
	//Q2: R: start row to r((start row + end row)/2) - C: r((start col + end col)/2)+1 to end col
	//Q3: R: r((start row + end row)/2)+1 to end row - C: start col to r((start col + end col)/2)
	//Q4: R: r((start row + end row)/2)+1 to end row - C: r((start col + end col)/2)+1 to end col
	split(start_row, ((start_row + end_row)/2), start_col, ((start_col + end_col)/2));
	split(start_row, ((start_row + end_row)/2), ((start_col + end_col)/2)+1, end_col);
	split(((start_row + end_row)/2)+1, end_row, start_col, ((start_col + end_col)/2));
	split(((start_row + end_row)/2)+1, end_row, ((start_col + end_col)/2)+1, end_col);

}


void calculate_adjacence() {

	//The computing of the adjacence for each region is the tricky part of the Split and merge algorithm.
	//After the split phase, we got into our regions vector all the splitted region of the image, of various size.
	//A good approach to find all the adjacence without mantaining a list of adjacence for every region, is to declare
	//a set STL container into each region object: with the property of the set, which is an associative container in
	//which each element is unique (cause the key of each element is the element itself) we can mantain a list of adjacence
	//for each point easily. 
	//How? 
	//We can start by sorting the regions vector from the smallest to the largest, or from largest to smallest or
	//simply LEAVING IT as it is. DEPEND on the image. Then, for each region, we need to find 
	//his adjacent neighbor, identified in the Nord, Sud, East and West one (The diagonal are not adjacent).
	//With that said, how can we achieve that? Sorted the vector to the smallest to the largest, we extract the first region
	//from the vector and calculate the midpoint of each of his side orientation: Nord, sud, east, west, in form of a Point(x,y).
	//We then add -1 or +1 to the rows or cols, based on which side we're facing, to retrieve a point of another region.
	//NOTE: -1 should result, sometimes, still in the midpoint. For better results, add -2 or 2.
	//Before checking if a point belong to a certain region, we must check another thing: the current analyzed region, should be
	//smaller or equal to the region being extracted from the list. Thats because if a region is smaller the the current one, that
	//means that the region extracted from the list is splitted into various region aswell, so on a certain side we could have N neighbor
	//instead of just one. We left the calculation of the adjacency of the smaller region to the smaller region itself, because we can't 
	//decide how many regions are on a single side. We must compare with JUST ONE, that is equal or bigger than the current region.
	//Then If the resultant point BELONG to one of the neighbor region (N,S,E,W), then the region is adjacent to the current region analyzed.
	//We then proceed to push both the region extracted from the list into the current analyzed region adjacency list, and also we're going to 
	//push the current analyzed region into the adjacency list of the extracted region to respect the adjacency criteria
	
	//Optional: sort
	//sort(regions.begin(), regions.end(), [](const Region& r1, const Region& r2) { return r1.pixel_num < r2.pixel_num; });
	
	//for the size of the regions list
	for(int i=0; i<regions.size(); i++) {
		
		//current region
		Region& analyzed = regions.at(i);
		
		//calculating the center of the region: since we're operating on squared images, we can extract the midpoint in this way:
		//calculate the center of the square by subtracting, for each row/col: (endIndex - startIndex)/2, and the distance from the center
		//to one of his border: given the simmetry of the square, we can use this distance to move and shift alongside the border.
		//Note: remember that Point, Rect, Circle and all the drawning object represent on the x coordinate the columns and on the y the 
		int cen_row = round((analyzed.area.y + analyzed.area.height)/2);
		int cen_col = round((analyzed.area.x + analyzed.area.width)/2);
		Point c(cen_col, cen_row);
		
		///Calculate the distance from a side (doesn't matter who: is symmetric) to the center. We choose the rows to do this.
		//We then do: (center.x/y - start_row / start_cols)
		//For example, if the region is a square from row/col 90 to row/col 180, the center is located in (135,135). 
		//We calculate the distance for the row as: (Center.y - region.area.y).
		//Where area is the Rect object of the region that holds the coordinate inverted as explained before. 
		//Then the calculus is: (135 - 90) -> 45. The distance from the center to any of his side is 45.
		int distance = (c.y - analyzed.area.y);
		
		//Building the midpoint: 
		//Nord: the midpoint on the nord side relay on the same column of the center, but on the start row of the region.
		//We subtract -2 from the row to shift from the current region to the nord one. 
		//So we use the same columns as the center,  and the center row (135) minus the distance (45) -> 90: the start row.
		//if exists. We still consider x as columns and y as row.
		Point N(c.x, (c.y - distance)-2);
		
		//East: the midpoint on the east side relay on the same row of the center, but on the start column of the region.
		//We subtract -2 to shift the column from the current region to the east one. 
		//So we use the same row of the center and the center column (135) minus the distance (45) -> 90: the start col
		Point E((c.x + distance)-2, c.y);
		
		//Sud: the midpoint on the south side relay on the same column of the center, but on the end row of the region.
		//We add +2 to the row to shift from the current region to the south one.
		//So we use the same column as the center, and the center row(135) plus the distance(45) -> 180.
		Point S(c.x, (c.y + distance)+2);
		
		//West: the midpoint on the west side relay on the same row of the center, but on the end col of the region.
		//We add +2 to the col to shift from the current region to the west one.
		//So we use the same row as the center, and the center column(135) plus the distance(45) -> 180.
		Point W((c.x - distance)-2, c.y); 
		
		//analyze all the possible adjacent regions
		for(int j=0; j<regions.size(); j++) {
			
			//skip the iteration if i==j
			if(i==j) continue;
			
			//extract a region
			Region& extracted = regions.at(j);
			
			//if the pixel are of the current region is <= than the neighbor (smaller area)
			if(analyzed.pixel_num <= extracted.pixel_num) {
				//if the extracted area contains one of those point
				if(extracted.contains(N) || extracted.contains(S) || extracted.contains(E) || extracted.contains(W)) {
					//The extracted region contains one of the N,S,E,W point: it's adjacent to the analyzed region.
					//Push back the region for both the analyzed and extracted region
					regions.at(i).adjacency.insert(&extracted);
					regions.at(j).adjacency.insert(&analyzed);				
				}		
			}	
		}
	}
	
	

}


void merge() {

	//Calculated the adjacence, merging is an easy process:
	//We must check, for each adjacent region to a current analyzed one, if the difference of their
	//average is lesser than a certain value.
	//If so, we proceed to assign the analyzed neighbour the average of the current analyzed region.

	//for the size of the regions
	for(int i=0; i<regions.size(); i++) {
		//current region is
		Region& region = regions.at(i);
		//analyze his neighbor
		for(auto& n: region.adjacency) {
			//if the distance is lesser than a threshold
			if(norm(n->average_pixel, region.average_pixel) < t) {
				//merge the regions assigning the average to the neigh
				n->average_pixel = region.average_pixel;
			}		
		}
	}

}



void build_image_rect(Mat& output) {

	//For each region in the regions vector
	for(int i=0; i<regions.size(); i++) {
		//Starting at the beginning row of the region (y because rect store x and y inverted into his structure, and height
    	//that represent the end row of the region)    	
		for(int x=regions.at(i).area.y; x<regions.at(i).area.height; x++) {
			//Starting at the beginning col of the region (x because rect store x and y inverted into his structure, and width
        	//that represent the end col of the region)
			for(int y=regions.at(i).area.x; y<regions.at(i).area.width; y++) {
				//Assign to the output image the average of that region.
				output.at<Vec3b>(x,y) = regions.at(i).average_pixel;		
			}		
		}
	}
		
}

void build_image_norect(Mat& output) {

	//For each region in the regions vector
	for(int i=0; i<regions.size(); i++) {
		//Starting at the beginning row of the region (y because rect store x and y inverted into his structure, and height
    	//that represent the end row of the region)    	
		for(int x=regions.at(i).area.y; x<=regions.at(i).area.height; x++) {
			//Starting at the beginning col of the region (x because rect store x and y inverted into his structure, and width
        	//that represent the end col of the region)
			for(int y=regions.at(i).area.x; y<=regions.at(i).area.width; y++) {
				//Assign to the output image the average of that region.
				output.at<Vec3b>(x,y) = regions.at(i).average_pixel;		
			}		
		}
	}
		
}








































