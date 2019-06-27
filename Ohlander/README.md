Simple example of Ohlander Algorithm for image clustering.<br>
Known Bug: The algorithm sometimes does not converge (infinite loop) with very coloured and detailed images,<br>
and images with large sizes. If the algorithm does not show *istantly* the solution, try running with an higher threshold.<br>

# Results
### Image 1
Apples on grass, 2 clusters found with 10 histogram threshold<br>
<img src="https://i.ibb.co/v1S6KmH/1.jpg" width="300">
<img src="https://i.ibb.co/n8wKGrb/output.png" width="300">

### Image 2
Apples on yellow background, 8 clusters found with 40 histogram threshold<br>
<img src="https://i.ibb.co/PmQ4ywr/mele.jpg" width="300">
<img src="https://i.ibb.co/NFgHkqx/output.png" width="300">

### Image 3
Building with grass, 12 clusters found with 20 histogram threshold <br>
<img src="https://i.ibb.co/k3JkCfj/km.png" width="300">
<img src="https://i.ibb.co/r5B2V0J/output.png" width="300">
