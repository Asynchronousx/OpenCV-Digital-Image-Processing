#include <opencv2/opencv.hpp>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>  

//Split and Merge: 
//The split and merge algorithm, like the region growing, is an algorithm that aim to the segmentation of the image.
//Differently from the region growing algorithm, that follow a bottom-up approach, the split and merge algorithm is a 
//top-down methodology.
//The algorithm itself is explanatory: At each pass, if the current analyzed region of the image does not respect a critery
//of omogenity P(Ri), the region is divided in four subregion, and each of those region are processed aswell, getting splitted
//if they do not respect the P(Ri) and so on, until the P(Ri) is respected or the minimum split size has been met.
//Each region splitted is pushed into a vector of regions splitted by this criteria.
//When all the splitting is done, we must analyze this vector of splitted region to merge all the similar region.
//How? we first calculate all the adjacence from a region to all his neighbor (N,S,E,O, excluding the diagonal neighbor)
//And then proceed to merge all the adjacent region with a similar P(Ri), after the adjacence calculation.
//We can calculate the P(Ri) with the variance of the entire region: given the average intensity of a given region, 
//Calculate the variance by the formul pow(average-pixel(x,y),2) or with the deviation standar abs(average-pixel(x,y)).
//The pixel(x,y) represent each pixel of the region that we would like to know the average.

using namespace std;
using namespace cv;

//Declaring public matrix that represent the raw image. So each object of the class region can acceed this matrix without
//having it as a member or accessing through reference pass.
Mat public_image;

//Class that represent the attributes of a region, such as his area covered (start col/row -> end col/row), 
//the average pixel value in that region, the number of pixel of that region and the label that represent the
//identifier of that region (useful for the adjacency property).
class Region {

	public:
	//Area represent the rectangle of the region
	Rect area;
	//Pixel num represent how many pixels are into the image
	int pixel_num;
	//Average represent the average value of the pixel of the region
	float average;
	//Variance represent the average variage value of the region
	float variance;
	//The adjacancy set contains all the adjacent region: we declare it as a set of pointer to region,
	//because when we'll merge the region by their average, we'll need to access to the exact structure 
	//to make the change permanent.
	set<Region*> adjacency;
	
	Region() {;}
	Region(int start_row, int start_col, int end_row, int end_col, float average=0, float variance=0) {
		//We intend the area as the two corner of the image: Top Left and Bottom right:
		//Note: Rect object is built by inverting rows and cols: x is the column, y is the row.
		//Same as the point structure.
		this->area = Rect(start_col, start_row, end_col, end_row);
		this->average = average;
		this->variance = variance;
		//Example: imagine that the area passed in input is: (sr,er) (sc,ec) -> (210, 219) (220,229).
		//There are 10 rows and 10 cols each. Since the last index is INCLUSE, we got that 219 - 210 is 9, and 229 - 220 is also 9.
		//We then must add 1 to make the sum inclusive aswell, and multiply the result: i.e (9+1) * (9+1) -> 100 element.
		this->pixel_num = ((area.height - area.y) + 1) * ((area.width - area.x) + 1);
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
	//Method to calculate the variance of the region, given a pixel in input.
	void compute_variance() {
		for(int i=area.y; i<=area.height; i++) {
			for(int j=area.x; j<=area.width; j++) {
				this->variance += pow(public_image.at<uchar>(i,j) - this->average, 2);
				//we could also calculate the variance as: this->variance += abs(this->average - public_image.at<uchar>(i,j));
			}		
		}
		
		this->variance /= pixel_num;
		
	}
	
	//Method to obtain the average of the region
	void compute_average() {
		for(int i=area.y; i<=area.height; i++) {
			for(int j=area.x; j<=area.width; j++) {
				this->average += public_image.at<uchar>(i,j);					
			}		
		}	
		
		this->average /= pixel_num;
		
	}
		
};

void SplitAndMerge(Mat, Mat&);
void split(int, int, int, int);
void calculate_adjacency();
void merge();
void build_image_rect(Mat&);
void build_image_norect(Mat&);

//Defining a vector of regions and the threshold
vector<Region> regions;
int t;
int pixels_limit;

int main(int argc, char** argv) {
	
	if(argc!=4) {
		cerr << "Usage: ./<programname> <imagename.format> <threshold> <pixel split limit>" << endl;
		exit(EXIT_FAILURE);	
	}
	
	Mat raw_image = imread(argv[1], IMREAD_GRAYSCALE);
	t = atoi(argv[2]);
	pixels_limit = atoi(argv[3]);
	
	if(raw_image.empty()) {
		cerr << "Image format not recognized." << endl;
		exit(EXIT_FAILURE);	
	}
	
		
	//Resize the raw image: this because we want to operate with a squared image to let the segmentation by
	//splitting be optimal. We choose the max of (rows, cols) to mantain the ratio.
	int opt_size = max(raw_image.rows, raw_image.cols);
	
	//Do any post processing here if necessary
	
	//Resizing the image to the desired RowxHeight 
	resize(raw_image, raw_image, Size(opt_size, opt_size));
	
	//Declaring our dest image with the same size and type of the raw image resized, filled with 0
	Mat dest_image(raw_image.size(), raw_image.type(), Scalar(0));
	
	//Applying a gaussian blur to smooth the image into omogeneous regions
	GaussianBlur(raw_image, raw_image, Size(5,5), 3, 3);
	
	//Assigning it to the public image: this image will be used from the region class to calculate variance and average
	//everytime a new region has been created.
	public_image = raw_image;
	imshow("original", raw_image);

	//Calling the function Split and merge: the only thing we need to pass to that function is the raw image itself,
	//and the destination. All the values will be calculated and fetched inside this function.
	SplitAndMerge(raw_image, dest_image);	
	
	//Display the dest image
	imshow("Split & Merge", dest_image);
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
	int rows = raw_image.rows-1;
	int cols = raw_image.cols-1;
	
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
	
	//Set the adjacency for all the region
	calculate_adjacency();
	
	//Merge
	merge();
	
	//After the merging, the new average value are present in each regions. 
	//Let's re-build our image.
	//build again to display in the main
	build_image_rect(dest_image_rect);
	build_image_norect(dest_image);
	
	//showing the images
	imshow("Split&Merge No Rect", dest_image);
	imshow("Split&Merge Rect", dest_image_rect);
	
	imwrite("out1.png", dest_image);
	imwrite("out2.png", dest_image_rect);
	
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
	if(region.pixel_num/4 <= pixels_limit) {
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

void calculate_adjacency() {

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
	//push the current analyzed region into the adjacency list of the extracted region to respect the adjacency criteria.

	//Sort the regions vector inline
	//sort(regions.begin(), regions.end(), [](const Region& r1, const Region& r2) { return r1.pixel_num > r2.pixel_num; });

	//for the entire regions vector, extracting the current analyzed region
	for(int i=0; i<regions.size(); i++) {
	
		//analyzed is the region being analyzed to all the other region	
		Region& analyzed = regions.at(i);
		
		//calculating the midpoint of the region: since we're operating on squared images, we can extract the midpoint in this way:
		//calculate the center of the square by subtracting, for each row/col: (endIndex - startIndex)/2, and the distance from the center
		//to one of his border: given the simmetry of the square, we can use this distance to move and shift alongside the border.
		//Note: remember that Point, Rect, Circle and all the drawning object represent on the x coordinate the columns and on the y the 
		//rows.
		int row_cen = round((analyzed.area.height + analyzed.area.y)/2);
		int col_cen = round((analyzed.area.width + analyzed.area.x)/2);
		Point c(col_cen, row_cen);
		
		
		//Calculate the distance from a side (doesn't matter who: is symmetric) to the center. We choose the rows to do this.
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
		Point E((c.x - distance)-2, c.y);
		
		//Sud: the midpoint on the south side relay on the same column of the center, but on the end row of the region.
		//We add +2 to the row to shift from the current region to the south one.
		//So we use the same column as the center, and the center row(135) plus the distance(45) -> 180.
		Point S(c.x, (c.y + distance)+2);
		
		//West: the midpoint on the west side relay on the same row of the center, but on the end col of the region.
		//We add +2 to the col to shift from the current region to the west one.
		//So we use the same row as the center, and the center column(135) plus the distance(45) -> 180.
		Point W((c.x + distance)+2, c.y); 
		
		//for the entire regions vector, extracting a region to compare with the current analyzed region from the vector
		for(int j=0; j<regions.size(); j++) {
			//if j == i that means we're analyzing the same region being analyzed. skip this iteration
			if(j==i) continue;
			
			//extracted is the i-th region extracted to compare to the analyzed region
			Region& extracted = regions.at(j);
				
			//TODO: check if neighborhood are correct		
			
			//check if the analyzed region is smaller or equal than the extracted region
			if(analyzed.pixel_num <= extracted.pixel_num) {
				//check if the extracted region contains any of the N,S,E,W point: if yes, the region are adjacent.
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
	
	//For the size of the regions vector
	for(int i=0; i<regions.size(); i++) {
		//Region is the current analyzed region
		Region& region = regions.at(i);
		//range-for on the adjacency set of the current region
		for(auto& neighbor: region.adjacency) {

		
		
			//If the difference between current and neighbor average is lesser than a value (30 in this case)
			if(abs(neighbor->average - region.average) < 30) {
				//Assign to the neighor the average of the current region 
				neighbor->average = region.average;
			}			
		}
	}

}


void build_image_rect(Mat& output) {

	
    //For each region in the regions vector
    for(int k=0; k<regions.size(); k++) {
    	//Starting at the beginning row of the region (y because rect store x and y inverted into his structure, and height
    	//that represent the end row of the region)    	
        for(int i=regions.at(k).area.y; i<regions.at(k).area.height; i++) {
        	//Starting at the beginning col of the region (x because rect store x and y inverted into his structure, and width
        	//that represent the end col of the region)
            for(int j=regions.at(k).area.x; j<regions.at(k).area.width; j++) {
            	//Assign to the output image the average of that region.
                output.at<uchar>(i,j) = regions.at(k).average;
            }
        }
    }
          
}

void build_image_norect(Mat& output) {

	
    //For each region in the regions vector
    for(int k=0; k<regions.size(); k++) {
    	//Starting at the beginning row of the region (y because rect store x and y inverted into his structure, and height
    	//that represent the end row of the region)    	
        for(int i=regions.at(k).area.y; i<=regions.at(k).area.height; i++) {
        	//Starting at the beginning col of the region (x because rect store x and y inverted into his structure, and width
        	//that represent the end col of the region)
            for(int j=regions.at(k).area.x; j<=regions.at(k).area.width; j++) {
            	//Assign to the output image the average of that region.
                output.at<uchar>(i,j) = regions.at(k).average;
            }
        }
    }
    
    public_image = output;
       
}
