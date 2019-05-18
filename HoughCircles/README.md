Important note: <b>Pre-Processing</b> with the Hough algorithms are important quite a lot. <br>
For example, without a good threshold on Canny or without an histogram equalization we could have sometimes<br>
worst result. Be careful on the Blur, Canny and other pre-processing threshold.<br>

Some samples of circle recognition with my Hough implementation:

OpenCV image with Threshold of 160, min radius of 20 and max radius of 60<br>
<img src="https://i.ibb.co/pvXVXdR/ocv.png" width="300">
<img src="https://i.ibb.co/J7qnNWG/h1.png" width="300">
<br>
Coins image with Threshold of 170, min radius of 20 and max radius of 50<br>
<img src="https://i.ibb.co/vdrkZP7/coins1.jpg" width="300">
<img src="https://i.ibb.co/L5Wv93V/h2.png" width="300">
<br>
Concentric circles withThreshold of 170, min radius fo 20 and max radius of 180 (Note: for the sake of the execution time,<br>
multiple Hough iterations with smaller radius would be computationally faster.)<br>
<img src="https://i.ibb.co/yycDmg1/conc.jpg" width="300">
<img src="https://i.ibb.co/drkmm61/h3.png" width="300">
<br>
