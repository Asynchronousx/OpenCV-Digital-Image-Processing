#include <opencv2/opencv.hpp>
#include <iostream>

using namespace cv;
using namespace std;

int main(int argc, char** argv) {



	double alpha = 0.5; double beta; double input;

 Mat src1, src2, dst;

 /// Ask the user enter alpha
 std::cout<<" Simple Linear Blender "<<std::endl;
 std::cout<<"-----------------------"<<std::endl;
 std::cout<<"* Enter alpha [0-1]: ";
 std::cin>>input;

 /// We use the alpha provided by the user if it is between 0 and 1
 if( input >= 0.0 && input <= 1.0 )
   { alpha = input; }

 /// Read image ( same size, same type )
 src1 = imread(argv[1]);
 src2 = imread(argv[2]);
 
 Size size(src1.cols, src1.rows);
 resize(src1, src2, size);

 if( !src1.data ) { printf("Error loading src1 \n"); return -1; }
 if( !src2.data ) { printf("Error loading src2 \n"); return -1; }

 /// Create Windows
 namedWindow("Linear Blend", 1);

 beta = ( 1.0 - alpha );
 addWeighted( src1, alpha, src2, beta, 0.0, dst);

 imshow( "Linear Blend", dst );

 waitKey(0);
 return 0;

	/*
	Mat raw1, raw2, dest;
	
	raw1 = imread(argv[1]);
	raw2 = imread(argv[2]);
	
	Size size(raw1.cols, raw1.rows);
	
	resize(raw1, raw2, size);
	
	double alpha = 0.5;
	double beta = (1.0 - alpha);
	
	addWeighted(raw1, alpha, raw2, beta, 0.0, dest);
	
	imshow("Blend", raw2);
	waitKey(0);
	return 0;
	*/
	
}
