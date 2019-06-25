#include <opencv2/opencv.hpp>
#include <iostream>
#include <math.h>

//Olhander Algorithm: A clustering algorithm for color image
//The Olhander algorithm bases his logic on a recursive method of the image analysis.
//The algorithm logic is based on the Histogram analysis of a given masked image at every i-th step.
//We build a stack containing BINARY MASKS from each histogram suddivision (Example: given an image, we
//split the image into his 3 channel and for each channel we calculate the histogram: we find the two peak of
//the histogram (if present), and choose the midpoint of the histogram: we generate then 2 binary mask for each
//histogram (binary mask lower than midpoint, and binary mask major than midpoint) and push all those binary masks
//into the stack. At each iteration, we pop a mask from the stack, and apply it on the image, calculating again his
//histogram etc, until no more histogram can be calculated because there are no more significative 2-peak to find.
//That means the histogram with no 2-peak should be considered as one of the histogram that segment our image into 
//regions.

//The approach is the following:
//We start with our original image, and a mask with all ones that we push into the stack. The applying of this mask
//on the image will produce the image itseflf. 
//We then proceed to extract the mask pushed into the stack, applying it to the image, and splitting the image into
//his 3-channel: R,G,B. For each resultant splitted image, we calculate the histogram of the relative channel, and proceed
//to threshold the histogram (of R,G,B) by finding his MIDPOINT between the TWO PEAK. We produce the two binary masks
//for each channel, and push them into the mask stack if the mask HAS the midpoint (more than one cluster).
//If the mask has no midpoint (no 2-peak) we put that into a resultant mask vector: that mask will isolate all the 
//pixel of a specific region.
//So, from each channel we should produce (if there are two peaks) two masks: those masks will be added to the mask
//stack an popped into the next iteration (or recursion), applied on the image, and from that masked image we will 
//compute the three histogram again and following the procedure written above.
//With those mask then, we active, respectively, all the pixel that are higher than the threshold, and all the pixel 
//that are lower than the threshold, isolating a specific region into our image.
//The algorithm will REACH CONVERGENCE EVERYTIME, because there will be a time when an histogram could not be splitted
//into two masks anymore. Those hisogram are our region-segmentation masks of a specific region.

using namespace std;
using namespace cv;

//Mask stack and result stack
vector<Mat> result_masks;
vector<Mat> mask_stack;
vector<Mat> result_vector;

//Variables for the calcHist function
int hist_size = 256;
float ranges[] = {0, 256};
const float* hist_ranges = {ranges};

//threshold for the minimum histogram splitting range
int t;

//BGR channels to apply masks and analyze histograms at every recursion, and the
//histogram matrix: since we do not need to calculate at every step the splitted 
//channels histogram, we can keep only the initial matrix and keep applying resultant
//mask on that: we then proceed to overwrite at each recursion the channels_hist to 
//calculate the new masks.
Mat channels[3];
Mat channels_hist[3];

//Public image
Mat public_image;

void Ohlander();
void generate_mask(Mat, int, int);
void calculate_midpoint(Mat hist, int&, int&);

int main(int argc, char** argv) {
	
	if(argc!=3) {
		cerr << "Usage: <./programname> <image.format> <histogram threshold>" << endl;
		exit(EXIT_FAILURE);
	}
	
	Mat raw_image = imread(argv[1], IMREAD_COLOR);
	t = atoi(argv[2]);
	
	if(raw_image.empty()) {
		cerr << "Image format must be valid." << endl;
		exit(EXIT_FAILURE);	
	}

	GaussianBlur(raw_image, raw_image, Size(5,5), 8);
	//assigning the image to the public image
	public_image = raw_image;
	
	//splitting the image into the bgr vector
	split(raw_image, channels);
	
	//creating a dest image
	Mat dest_image(raw_image.size(), raw_image.type(), Scalar(0));
		
	//Generating the first mask to push into the pixel stack: at the beginning we want to consider the whole image,
	//so let's generate a all-ones mask: that means the mask got all his pixel active.
	Mat initial_mask(raw_image.size(), CV_8UC1, Scalar(1));
	
	//pushing the mask into the vector: initially, the vector will contain only one mask, with all ones.
	mask_stack.push_back(initial_mask);
	
	//calling ohlander
	Ohlander();
	
	//Once here, Ohlander finished his procedure. We now start, for each mask into the result vector, to
	//assign a random color to that region.
	RNG rng(12345);
	for(int i=0; i<result_vector.size(); i++) {
		Vec3b color = Vec3b(rng.uniform(0,255), rng.uniform(0,255), rng.uniform(0,255));
		for(int x=0; x<raw_image.rows; x++) {
			for(int y=0; y<raw_image.cols; y++) {
				if(result_vector.at(i).at<uchar>(x,y) != 0) {
					dest_image.at<Vec3b>(x,y) = color;
				}
			}		
		}
	}
	
	
	//show olhander
	imshow("Original", raw_image);
	imshow("Ohlander", dest_image);
	waitKey(0);

}




void Ohlander() {

	//If the mask stack is empty, we do not have nothing more to do
	if(mask_stack.empty()) {
		return;	
	}
	
	//maybe clone local copy of Mat channels
	//What if itsn't empty? We re going to pop a mask from the vector and remove it.
	Mat mask = mask_stack.back();
	mask_stack.pop_back();

	//We then apply calcHist to R,G,B channel and calculating the histogram. We can apply this
	//mask directly on the histogram to have, in result, the three masked histogram that we'll send
	//to the threshold find function.	
	//We already got all the variables needed to achieve that, so we're going to calculate
	//the three histogram applying the mask popped from the stack on them. We can achieve that
	//by passing to calcHist the current mask. Calc hist takes in input those argument:
	int valley, peak, index, max_valley, max_peak = 0;
	bool found = false;
	for(int i=0; i<3; i++) {		
		calcHist(&channels[i], 1, 0, mask, channels_hist[i], 1, &hist_size, &hist_ranges);
		calculate_midpoint(channels_hist[i], valley, peak);
		
		//if the midpoint couldn't be found (valley = -1), that means we should consider the current analyzed 
		//masked channel as one of the results mask of the image. let's push it into the resulting mask
		if(valley <= 0) {
			//push into the result mask
			result_masks.push_back(mask);
		} 
		//else, we need to fetch only the highest-peak histogram to send it to the generation of two masks
		else if(peak > max_peak) {
			//fetch infos and found = true
			found = true;
			index = i;
			max_valley = valley;
			max_peak = peak;	
		}		
	}
	
	//once here we got two cases: we found something to push, or we didnt. 
	//if we found something, lets generate his submasks
	if(found) {
		generate_mask(channels_hist[index], max_valley, index);		
	} 
	//if we found nothing, we need to do two things: its possible that the stack mask still contain something.
	//lets assign to the result vector the mask stack first, because those mask are the "MOST GREZZE" and so need
	//to be pushed before the "accurate mask" found into the result mask vector.
	//we then append the result mask vector
	else {
		result_vector = mask_stack;
		result_vector.insert(result_vector.end(), result_masks.begin(), result_masks.end());
		return;
	} 
	
	
	//calling again olhander to check if there are still mask to analyze into the image
	Ohlander();

}

void calculate_midpoint(Mat hist, int& valley, int& peak) {
   
    //Applying a strong gaussian filter on the histogram that is going to be analyzed: that's because, an histogram
    //could present a lot of peaks into it. With that gaussian blur, we're going to transform him into a possible
    //bimodal histogram, composed only by two peak: from that, is easy to find the valley (and the two peak)
    GaussianBlur(hist, hist, Size(55,55), 0, 0);
   
    //Declaring a map formed by, at the key the histogram occurrences, and as value the index at which the occurrences
    //are found.
    map<int, int> hist_map;
   
    //inserting into the map all the occurrence of the histogram with their relative index
    for(int i=0; i<hist.rows; i++) {
        hist_map.insert(make_pair(hist.at<float>(i,0), i));
    }
   
    //Let's find the two peaks based on the peak approximation of the valley, and recompute the exact valley from
    //the index of the two peak found.
    int fmax = 0, smax = 0;
   
    //first scanning for the firt max, from 0 to the hist rows: since the histogram is a Mat(256,1) we know that each column will
    //be 0.
    for(int i=0; i<hist.rows; i++) {
        //if the current element is <= than the max, that means we reached the max value of a possible peak,
        //so we take the index found until now and
        if(hist.at<float>(i, 0) >= fmax) {
            fmax = hist.at<float>(i, 0);
        } else break;
    }
   
    //After finding the first max, we know that we should iterate at most to the first maximum index found: this should
    //not happens if there are two peak into the histogram.
    for(int i=hist.rows-1; i>hist_map[fmax]; i--) {
        //same approach as above
        if(hist.at<float>(i,0) >= smax) {
            smax = hist.at<float>(i, 0);
        } else break;
    }
   
    //the index of the peak are found at the fmax/smax key value: calculate the new exact valley: but first,
    //check if the distance between a peak and another is enough significative to split the histogram:
    if(abs(hist_map[fmax] - hist_map[smax]) < t) {
    	valley = -1;
    	peak = -1;
    } else {
		//if here, the difference from the peak is significative.
		//we just need to return the index of the found valley and the peak
		valley = (hist_map[fmax] + hist_map[smax]) / 2;
		peak = max(fmax, smax);    
    }
    
}

void generate_mask(Mat hist, int valley, int index) {

	//otherwise, valley represent a real index at which the histogram should be split.
	//Let's split the histogram into two subhistogram, from which we'll generate a segmentation mask.
	//To generate the left and right histogram, we clone the original one into both, and deactivate
	//the pixel above or lesser than the valley threshold.
	Mat hist_left = hist.clone();
	Mat hist_right = hist.clone();
	
	
	//divide the hist passed in input in two regions: left and right, using the valley index as delimiter
	//Left region
	for(int i=valley; i<hist.rows; i++) {
		hist_left.at<float>(i,0) = 0;
	}
	
	//Right region
	for(int i=valley; i>=0; i--) {
		hist_right.at<float>(i,0) = 0;
	}
	
	//Generating the two mask, big as the image size, initially with all the pixels unactive(0)
	Mat left_mask(channels[index].size(), CV_8UC1, Scalar(0));
	Mat right_mask = left_mask.clone();
	
	//Now, we're going to analyze all the image to check whenever a pixel in the left or right mask should be
	//activated (1) or deactived(0). 
	//How can we achieve that?
	//Remember that, an histogram is relative to a specific channels, B, G or R. Analyzing the B histogram, we should
	//analyze the B channel of the image, Analyzing the R histogram the R channel etc. We can check which channel we're 
	//analyzing by using the index passed in input.
	//So we iterate through the spllited channel image, and check for the intensity level of every given pixel. 
	//Since the left and right histogram will present a region with all zeros and a region with all values > 0, 
	//we can use the values of the intensity of the channels to check if a specific value is > 0 or not. If
	//with that, we'll build or mask, activaging and deactivating the right pixels.
	for(int i=0; i<channels[index].rows; i++) {
		for(int j=0; j<channels[index].cols; j++) {
			//fetching the intensity of the pixel at (i,j)
			int intensity = channels[index].at<uchar>(i,j);
			//using that intensity to check if the specific value into the specific histogram is > 0. 
			//If so, should be activated into the relative mask, otherwise deactivated.
			left_mask.at<uchar>(i,j) = hist_left.at<float>(intensity, 0) > 0? 1: 0; 
		}	
	}
	
	//and doing the same process for the right mask
	for(int i=0; i<channels[index].rows; i++) {
		for(int j=0; j<channels[index].cols; j++) {
			//fetching the intensity
			int intensity = channels[index].at<uchar>(i,j);
			//and using the intensity as index to access the specific histogram, to check
			//if a value is deactivated or not.
			right_mask.at<uchar>(i,j) = hist_right.at<float>(intensity, 0) > 0? 1: 0;		
		}
	}

	//push the built masks into the mask stack
	mask_stack.push_back(right_mask);
	mask_stack.push_back(left_mask);

	
	
}









