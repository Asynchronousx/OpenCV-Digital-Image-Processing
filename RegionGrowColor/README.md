##### Region Growing Color

A simple approach to the fundaments of the image segmentation, the region growing algorithm for the color images.<br>
This algorithm is implemented both with the analysis of the variance, that offers a more complex approssimation<br>
of the average for each region, and the analysis of the intensity, which will choose if a pixel belong or not to a region<br>
based on the difference of the pixel intensity.<br><br>
The image showed here are the output of a <b>random seed</b> approach.<br>
<br>Note: the quality input paramter means the number of the iterations wanted * 1000. Also,<br>
Also, one pass of median filter is applied. More passing of median filter should result in a better image.<br>

For the sake of the example, both are showed with their input, output and threshold values.

#### Image 1 with variance analysis
NY city image with a Threshold of 50 and quality of 10.<br>
<img src="https://i.ibb.co/SfPcn4t/ny.jpg" width="300">
<img src="https://i.ibb.co/0yMPXDz/nycva.png" width="300">
<br>

#### Image 1 with intensity analysis
NY city image with a Threshold of 80 and quality of 10.<br>
<img src="https://i.ibb.co/SfPcn4t/ny.jpg" width="300">
<img src="https://i.ibb.co/4PXrgXj/nycin.png" width="300">
<br>

#### Image 2 with variance analysis
Mixed fruits image with a Threshold of 50 and quality of 10.<br>
<img src="https://i.ibb.co/89JNrP0/mixed-fruits.jpg" width="300">
<img src="https://i.ibb.co/sQtfnL1/mixedvar.png" width="300">
<br>

#### Image 2 with intensity analysis
Mixed fruits image with a Threshold of 80 and quality of 25.<br>
<img src="https://i.ibb.co/89JNrP0/mixed-fruits.jpg" width="300">
<img src="https://i.ibb.co/kMhtyps/mixedint.png" width="300">
<br>
