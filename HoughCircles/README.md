Important note: <b>Pre-Processing</b> with the Hough algorithms are important quite a lot. <br>
For example, without a good threshold on Canny or without an histogram equalization we could have sometimes<br>
worst result. Be careful on the Blur, Canny and other pre-processing threshold.<br>

Some samples of circle recognition with my Hough implementation:

OpenCV image with Threshold of 160, min radius of 20 and max radius of 60

<br>
Coins image with Threshold of 170, min radius of 20 and max radius of 50

<br>
Concentric circles withThreshold of 170, min radius fo 20 and max radius of 180 (Note: for the sake of the execution time,<br>
multiple Hough iterations with smaller radius would be computationally faster.)
