/--------\
| OPTIONS |
\--------/

GENERAL
------
-V / --verbosity (int) : 1
	Verbosity level, from 0 (silent - no output) to 4 (heavy output, commenting everything)

-I / --input (string) : ""
	The input-file which is processed, every common 2D-image file should be
	readable. The original image won't be overwritten.

-O / --output (string) : "/tmp/"
	The output-directory. Note: the output file will have the same name but
	additionally a suffix indicating the output behavior specified by the
	various --write* commands.

--log-times (boolean) : false	
	Writes the times of different sections to a log-file (specified by --log-file)

--log-file (string) : "" 
	File where the measured times are stored to.

--graphical	(bool) : false
	Shows some graphical output, not really integrated everywhere. It may not
	show what you'd like to see. 

--suffix (string) : ""
	Appends a suffix to the output filename which is then 'basename' +
	'internal_suffix' (according to the write command --write*) + '.png'

--groundTruthFile (string): ""
	Used only in combination with --writePairs to mark pairs which are correct
	in different colors.  

--numThreads (int) : 0
	Number of threads used for feature extraction (only for block-based features). 
	Note: 0 = use all possible threads.


PREPROCESSING
-------------
--checkEntropy (bool) : false
	Should the entropy be checked -> blocks with too low entropy are removed.
	For Bravo Solorio features (--featvec=BRAVO) please use just an entropyTh >
	0.0 and don't use this option.

--loadCorresFile (string) : ""
	Load a saved correspondence-file (*.ser) and omits preprocessing, feature
	extraction and matching. The saved matches (stored as 16bit pngs are loaded
	according to the .ser-File - Don't rename any of them!)

--writeCorresFile (string)	: ""
	Write a correspondence-file (*.ser + *.png), the most options are saved in
	the ser-file and the matches are stored as 16bit png images. Can later be loaded.

--chan (enum: ALL, GRAY, RED, GREEN, BLUE) : GRAY
	Chose a specific channel or convert to a gray-scaled version. Standard is
	GRAY, i.e. a the image is converted to grayscale.  Note: Some features like
	BRAVO, LUO need all color channels.

--pyramidLvl (int) : 0
	You can set a pyramid-level, a gaussian pyramid will be build up to this
	level and the highest (= level specified) will be chosed as new image.

--waveletLvl (int)	: 0 
	You can set a wavelet-level, a wavelet decomposition (Haar = Daubechies D2)
	will be build up to this level and the highest (= level specified) will be
	chosen as new image.

-B --blockSize (int) : 16
	One dimension of a square which will be used as the block / window which is
	then slided over the window with a certain stepsize (--step). 

--step (int) : 1
	The stepsize in which the windows are shifted over the image. 

--roi (vector<int> ) : empty
	Extract blocks only from this region of interest (doesn't work for
	keypoint-based features)	
	Form: x1 y1 x2 y2, the pair x1 and y1 have to be one corner of the window
	the other pair has to be the opposite corner (diagonal) you can
	specify more than one window, every 4 values are one region of
	interest.


FEATURE-EXTRACTION
------------------
--featvec (enum: LUO, BLUR, DCT, PCT, SVD, HU, DWTFEAT, FMT, LIN, BRAVO, CIRCLE, ZERNIKE, NO, NOTHING) : NO
	The feature vectors for the block-based methods which can be used. See the
	related publications for informations.  Note: 'NO' specifies the raw image
	block. If 'NOTHING' is chosen really nothing happens.
		 
--detectorType (string) : ""
	All possible detector-types from the OpenCV-library.

--descriptorType (string) : ""
	All possible descriptor-types from the OpenCV-library.

--normalize (bool) : false
	Normalizes the feature vector between 0 and 1.

--maxFeatureSize (int) : 512
	Truncate the feature's size to the given value. Only for block-based methods.

--minFeatureSize (int) : 0
	Size a feature should have. Only for (K)PCA.

--modified (bool) : false
	Some features have a second version implemented which can be specified by this value.
	Currently a modified version of the BLUR-features and SVD-features are implemented.

--circle (bool) : false [deprecated / needs revision]
	If true, only a circle around the midpoint of the block will be used for feature calculation.
	Note: only for block-based features.

--qf (double) : 0.5
	Quality-factor for DCT if --featvec=DCT

--numHU (int): 4
	Number of Hu-Moments used if --featvec=HU

--decomposeLvl (int) : 0
	Decomposition-level of DWT-features if --featvec=DWT

--degZernike (int) : 5
	Degree of Zernike-moments if --featvec=ZERNIKE.

--dstwidth (int) : 6
	Size of the destination width of the log-polar transformation if --featvec=FMT

--usePCA (bool) : false
	Apply a PCA after feature extraction.

--dim (int) : 0
	New dimension for PCA (either feature or PCA after feature extraction) or
	SVD. 
	If set to '0', than this will computed dinamically according to --eps.
	If set to '-1' (only SVD) then the dimension will be the average dimension of all blocks.
	If set to '-2' (only SVD) then the dimension is computed every time.
		Features with size smaller than the maximal dimension will be padded with zeros.

--eps (double) : 0.01
	Fraction of ignored variance along the principal axis - only used for PCA/SVD if --dim=0.

--sigma (double) : 25.0
	Kernel-variance for KPCA features (--featvec=KPCA).

--numSamples (int) : 144
	Number of samples used for training the KPCA. 

--useCps (bool) : false
	Uses cross power spectrum to search for duplicated regions.	
	Note: This method doesn't work and has to be revised.
	Note: you have to choos NO as featvec!

--partsInX (int) : 2
	For cps: parts in x direction

--partsInY (int) : 2
	For cps: parts in y direction

--correlTh1 (double) : 0.025 
	For cps: correlation threshold (in feature).

--correlTh2 (double) : 0.2
	For cps: correlation threshold (in matching for crit=CORREL).

--pixelDiff (int) : 2
	For cps: maximal difference between similar pixel")


QUANTIZATION
------------
General note: this is not a histogram-based quantization.

--qStart (int) : 0
	Start of quantization of the featurevector (zero based).

--qEnd (int) : 0
	End of quantization of the featurevector (zero-based and not including the end).
	If is set to 0 , no quantization is performed.
	If set to -1, the whole feature will be quantized - qEnd is set to featureSize.

--qNum (int) : 1
	Number of quantization bins in which of the featurevector is quantisized - so if set to 1, nothing changed.

--qRound (bool) : false
	If you quantize you can set this option to round before quantization - will be rounded down (floor)

--qRound2 (bool) : false
	If you quantize you can set this option to round after quantization - will be rounded down (floor)

--qAbs (bool) : false
	If you quantize you can set this option so take the absolute value of the features.



MATCHING
--------
--numRows (int) : 1
	Attention: This is not only for row-wise sorting.
	It denotes the number of nearest neighbors which will be tested further.
	So, if --kdsort=false && --nnBrute=false it denotes the number of consecutive rows of sorted blocks which will be tested further.
	Otherwise, it will check the next nearest neighbor according to the other criteria further.

--minDistEuclidian (int) : 50
	Denotes the minimum Euclidian distance two blocks need to be away from each other in image space.

--minDistCheby (int) : 0 
	Denotes the minimum Chebychev distance two blocks need to be away from each other in image space.
	Can be used instead of the minimum Euclidian distance or additionally to it.

--crit (vector<enum>: NONE, LUOCRIT, EUCLIDIAN, CPS, CORREL, BRAVOCRIT, EUCLIDIAN_RATIO)
	The similarity criterion to decide which block (respectivly the featurevector of a block) is similar to another\n"
	Notes: 	LUOCRIT (should only be used in combination with featvec LUO) -> then you have to set 9 thresholds (--th).
			EUCLIDIAN computes eculdian distance in featurespace and compares it with the threshold (--th).
			CPS computes cross power spectrum between 2 blocks and compares the max with phaseCorrelTh.
			CORREL computes correlation coefficient of the fourier magnitudes of the log polar trafo of a block.
			BRAVOCRIT compares the first 3 avg colors with thresholds and then performs CORREL criterion.
			EUCLIDIAN_RATIO compares the ratio of the distance between the i and the i+1 feature vector. 

--th (vector<double>)
	Threshold(s) for the similarity criterion.

--euclidianTh (bool) : true
	If set to true: the first Threshold (--th) will be seen as an euclidian distance in
	feature space, otherwise it is seen as a similarity measurement in percent.

--kdSort (bool): false
	Sort the features according to approximated nearest neighbor search (FLANN) or with lexicographic sorting (if false).

--nnBrute (bool) : false
	Compute the nearest neighbor by brute force - should only be used with keypoints or small images.

--entropyTh (double) : 0.0
	Entropy-threshold for the method of Bravo-Solorio et al. (crit=BRAVOCRIT)

--phaseCorrelTh (double) : 0.0
	Phase-correlation threshold if crit=CPS.


VERIFICATION
------------
--cluster (bool) : false
	 Use clustering to group near features which lie close to each other. You
	 should only use it with keypoint-based methods

--cutTh (double) : 0.0
	Cut threshold for clustering. 0.0 means that this threshold is adjusted
	automatic depending on the number of matches and the method used

--normalizeCluster (bool) : false
	normalize the cluster points such that the points are translated to the
	cluster-center and scaled that avg-distance = sqrt(2) --estimateMethod",

--estimateMethod (string: 'ownransac', 'ransac', 'lmeds', 'pure', 'goldstandard') : ""
	Method to find a valid transformation between copied and pasted regions.
	OR: Method to finetune already found transformations.
	That means there exist 3 possibilities: 
	1. use the method after clustering (as Amerini et al. did it)
	2. use the method after fastsats to finetune parameters
	3. use method on the correpondence maps (as Pan et al. did it)
	'ownransac' is an own implementation of RANSAC in contrast to 'ransac' which uses
	the OpenCV-implementation. 'lmeds' is searches for a homography using LMEDS
	(least median robust method). 'pure' uses a regular method using all the points. 
	'goldstandard' computes the homography matrix according to Hartley. 

--ranReprojTh (double) : 3.0
	RANSAC reprojection threshold if --estimateMethod=ransac.

--ranIterations (int) : 100
	Number of iterations used for RANSAC.

--fastsats (bool) : false
	Fast same affine transformation selection on shift vectors.

--withTree (bool) : false
	Fast SATS with kd-tree, works better for keypoint based methods. 

--maxDist (int) : 5
	Maximum distance another block should be away from the analyzed one (used in fastsats).

--recomputationTh (int) : 2
	Recompute transformation matrix when recomputationTh points are added (used in fastsats).

--checkShift (bool) : false
	Check shift vectors according to their similarity, i.e. build a histogram of all shift vectors and take 
	either all which are higher than --minSameShift or take the --numMax ones. 

--minSameShift (int) : 50
	How many blockpairs have to have the same (or similar if fastsats is used)
	shiftvector. Can be used for fastsats or with checkshift.

numMax (int) : 0
	Number of maximum shift vectors/

shiftVariance (int) : 0
	Variance between shift vectors (as stated by Bravo-Solorio et al.).

--computeCorrelMap (bool) : false
	Compute a correlation map with the help of earlier computed transformation
	matrices. See Pan et al .


MARKING (writing output files)
-------
--writeChrom (bool) : false
	Writes an image where the block positions is color coded as chromaticity-value.

--writePairs (bool) : false
	Marks block pairs, visualizes pairs in GT differently (needs --groundTruthFile).

--writeOutput (bool) : false
	Writes an ouput image with the non-duplicated regions are made darker. 
	Currently, this option doesn't seem to work.						  

--writeMatrix (bool) : false
	Writes down the matrix, a value != 0 in the matrix indicates a match.

--writeCluster (bool) : false
	Writes an image, in which every cluster-point has a different color.

--writeCorrelation (bool) : false
	Writes the not further prozed binarized correlation map. 
	Note: works of course only if we have a transformation matrix. 

--writePost (bool) : false
	Write the post-processed output matrix. 

--scaleUp (bool) : false
	If img is downscaled, e.g. due to gaussian pyramid, this scales it up again.

--useOrig (bool) : true
	Use original image instead of modified one e.g. due to gray-scale conversion
	and marks then the blocks in the original image (must be
	given via -I). 

--minVal (int) : 0
	only pixels greater than this value will be marked, makes only sense if you
	can have more than one match at a location. 
	Note: currently not used
	
--markRegions (bool) : false
	Marks forged regions each in different colors in the image.

--MarkShift (bool) : false
	Marks shift vector in images. 

--markContour (bool) : false
	Draw the countour of the areas in the original image. Looks pretty when
	having post-processed areas.

--binaryTh (double) : 0.3
	Threshold for binarization of the correlation maps.
	
POST-PROCESSING
---------------
Options which stem from the post-processing module.
These options can also be used independent from the cmfd-command (by removing
the pp-prefix), i.e. instead of 'bin/vole cmfd -I tree.png pp.rmSmallAreas...'
you can apply the post-processing to the output of the cmfd-pipeline:
	'bin/vole pp -I tree_matrix.png --rmSmallAreas...'

Note: If you use the correlation-maps you should apply the postprocessing to it, 
	but can additionally provide the match-matrix for further verification that each region
	must contain a matched point, i.e.:
	'bin/vole pp -I tree_correl.png --matrix tree_matrix.png --rmSmallAreas ...'
	This is done automatically if you use correlation-maps during cmfd.

--pp.rmSmallAreas (bool) : false
	Removes small areas using the areaThreshold and/or if a matrix is passed
	(which is the case for cmfd) it removes areas in which no match occurs. 

--pp.areaThreshold (int) : 0
	Duplicated areas must have at least the size: if areaThreshold (0,1] -> size =
	areaThreshold * image.width * image.height, else size = areaThreshold

--pp.fillHolesByContour (bool) : false
	Filles the holes regarding only the most extreme outer contours.

--pp.morphOp (vector<string>) 
	Morphological operations which can be applied (ERODE, DILATE, OPEN, CLOSE,
	ADIEND, TOPHAT, BLACKHAT).  See OpenCV-documentation what they do.

--pp.morphIter (vector<int>)
	How many iterations each morphological operation transformation should be 
	applied. Default: 1 iteration

--pp.matrix (image)
	See note above. Normally not neccessary when using 'vole' with the 'cmfd' command.
