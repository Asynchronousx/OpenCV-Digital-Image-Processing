#include <opencv2/opencv.hpp>
#define L 256

using namespace cv;
using namespace std;

//The hit or miss transorm is useful to find pattern of object into a given binary image.
//This morphological operation is useful also to reduce "binary noise", expressed as discontinued
//dots into the image. The logic of this algorithm is to convolve two structuring element, B1 and B2
//on an Image A at the same time, to find all the MATCHING elment of the first structuring element B1
//and all the UNMATCHING element of the filter B2. Then we INTERSECT (bitwise and of the images) the two result
//to obtain the unconnected pixels into the image.
//The formula is:
//A(x)(B1,B2) = (A(-)B1) AND (Ac(-)B2), where Ac is the complementary of A (the negative).
//Also, we can express the two structuring elements B1 and B2 as ONE filter, which will have three kind of pixels
//into it: Foreground pixels (1), Background pixels (-1), and irrelevant pixels. (0)
//The erosion, from the theory, will delete all the pixel that fall into the anchor point who's neighbour
//does not fall properly into the structuring element active pixels. In this case, the active pixels are
//the BG and FG one.

int main(int argc, char** argv){

	//Creating a custom input image in pixels, to understand better about the Hit or Miss transform.
    Mat raw_image = (Mat_<uchar>(8, 9) <<
       	0, 0, 0, 0, 0, 0, 255, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 255, 0, 255, 255, 255, 0,
        0, 0, 0, 255, 0, 255, 255, 255, 0,
        0, 0, 255, 255, 255, 255, 255, 255, 0,
        0, 0, 255, 255, 255, 255, 255, 255, 0,
        0, 255, 255, 255, 255, 255, 255, 255, 0,
        0, 255, 255, 255, 255, 255, 255, 255, 0);

    //Defining a single kernel B (union of B1, B2) to use into the morphologyEx operation.
    Mat B = (Mat_<int>(3,3) <<
    -1, -1, 0,
    -1,  1, 1,
     0,  1, 0);  //structuring element for searching the pattern of a top left angle

	//Let's indicate with F the foreground pixel, B the background and I the irrilevant ones.
	//This kernel B is the union of two kernel B1 and B2:
	/*
		|I F I|			|I I I|	   |I F I|
	 B1 |F I F|  and B2 |I B I| -> |F B F|
		|I F I|			|I I I|    |I F I|
	
	This special pattern, will find all the background pixel contained in a Foreground region (B1).
	
	For example, a pattern to find the upper left angles of an object would be
	
		|B B I|			|I I I|    |B B I|
	 B1 |B I I|  and B2 |I F F| -> |B F F|
		|I I I|			|I F I|	   |I F I|
		

	Example: if we're into a region area of an image A composed like this:
	
		        |B I B|										 |B B B|
	region of A |B F F| and a structuring element like this: |B F F|
			    |I F I|										 |I F I|
			    
	The central element of the image will be erased because it don't match the structuring element pattern
	into the element A(0,1) (is I, while the structuring element (0,1) got B. Then the pixel A(centerX, centerY)
	will be erased by the erode operation.
	
	Since the Kernel is anchored by default at the center, the result of the hit or miss will return EXACTLY
	the point centered into the kernel, ergo the point who match our pattern searching thanks to the erosion
	property.
	
	*/    	
    
    Mat dest_image;
    
    morphologyEx(raw_image, dest_image, MORPH_HITMISS, B);
    
    //Showing the images
    const int scaleRate = 10;
    
    //Set the kernel color to display it properly
    B  = (B  + 1) * 127;
    
	//Scaling the kernels
	resize(B, B, Size(), scaleRate*5, scaleRate*5, INTER_NEAREST);
    
    //Scaling the input and dest images
    resize(raw_image, raw_image, Size(), scaleRate*5, scaleRate*5, INTER_NEAREST);
    resize(dest_image, dest_image, Size(), scaleRate*5, scaleRate*5, INTER_NEAREST);

    //Showing
    //Before printing B, we need to convert B to a image-readable kernel.
    B.convertTo(B, CV_8U);
    imshow("Kernel B", B);

    imshow("Raw Image", raw_image);
    imshow("Hit or Miss", dest_image);
    waitKey(0);
    return 0;

}

