A simple approach to the fundaments of the image segmentation, the region growing algorithm.<br>
This algorithm is implemented both with the analysis of the variance, that offers a more complex approssimation<br>
of the average for each region, and the analysis of the intensity, which will choose if a pixel belong or not to a region<br>
based on the difference of the pixel intensity.<br>
Note: the quality input paramter means the number of the iterations wanted * 1000. Also,<br>
Also, one pass of median filter is applied. More passing of median filter should result in a better image.<br>

For the sake of the example, both are showed with their input, output and threshold values.

#### Image 1 with variance analysis
Coins image with a Threshold of 80 and quality of 10.<br>
<img src="https://i.ibb.co/Kw2NH0H/circle.jpg" width="300">
<img src="https://i.ibb.co/LNDgmjH/cirvar.png" width="300">
<br>

#### Image 1 with intensity analysis
Coins image with a Threshold of 80 and quality of 10.<br>
<img src="https://i.ibb.co/Kw2NH0H/circle.jpg" width="300">
<img src="https://i.ibb.co/c2462ZV/cirint.png" width="300">
<br>

#### Image 2 with variance analysis
Geometric shape with a Threshold of 30 and quality of 10.<br>
<img src="https://i.ibb.co/WW339Pz/rg.png" width="300">
<img src="https://i.ibb.co/85XKx7M/destrgint.png" width="300">
<br>

#### Image 2 with intensity analysis
Geometric shape with a Threshold of 30 and quality of 10.<br>
<img src="https://i.ibb.co/WW339Pz/rg.png" width="300">
<img src="https://i.ibb.co/2nMLHx3/destrg.png" width="300">
<br>
