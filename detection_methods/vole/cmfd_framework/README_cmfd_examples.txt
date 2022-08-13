Examples:
---------

1. Block-based example:

bin/vole cmfd -I tree.png -O . --chan GRAY --featvec ZERNIKE --kdsort --fastsats --minSameShift 1000 --markRegions --useOrig

Load image tree.png, apply grayscale conversion, extract features using Zernike moments, match the nearest neighbors, and
verify the matches with fastsats. Finally, mark the regions in the original image and save it in the current folder.


2. Keypoint-based example:

bin/vole cmfd -I tree.png -O . --chan GRAY --detectorType SURF --descriptorType SURF 
			--cluster --numRows 10 --th 0.6 --crit EUCLDIAN_RATIO --useOrig --blockSize 1 --nnBrute --normalizeCluster
			--estimateMethod ownransac --ranReprojTh 0.05 --ranIterations 2000 --computeCorrelMap --binaryTh 0.4
			--pp.rmSmallAreas --pp.fillHolesByContour --pp.areaThreshold 1000 --writePost --markContour

This converts the image horses_copy.png to grayscale, computes keypoints and features of SURF, gets the nearest neighbor by brute force, 
and checks the 10 nearest neighbors if the euclidian-ratio is below 0.6 according to the g2nn method of Amerini et al., 
and clusters these results (cutThreshold is adjusted internally). These clusters are normalized and the transformation
between the clusters is estimated by our implementation of the RANSAC method using 2000 iterations and a reprojection-threshold
of 0.05. Now a correlation map is computed (see Pan et al.) with a binarization threshold of 0.4, small areas are removed,
and holes in the contour are filled. Areas containing less than 1000 pixels are removed. From these areas, only the contour
is marked and the output file is written out to the current directory.


3. Save correspondence file using config-file + load correspondence-files
   (thus splitting up the computation process, good solution to defer the execution of the first (expensive) part to a compute cluster)
a)	bin/vole cmfd zernike.conf -V 0 -I tree.png --numThreads 4 --writeCorresFile
Computes the zernike features on the image using the zernike.conf config-file and saves the matches in a correspondence file called
'tree_info.ser'.

b)	bin/vole cmfd -V 0 --loadCorresFile tree_info.ser --fastsats --minSameShift 1000 --writeMatrix 
Load the matches by giving the ser-file as parameter, now compute the verification step and write out the computed cmfd-map
which can be compared with the ground-truth
