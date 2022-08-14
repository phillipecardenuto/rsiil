/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "copymoveframework_config.h"
//#include "post_process/post_process_config.h"

using namespace boost::program_options;

ENUM_MAGIC(featureType)
ENUM_MAGIC(similarityCriterion)

CmfdConfig::CmfdConfig(const std::string &prefix)
	: Config(prefix) , post("pp") // init post processinc config

{
#ifdef WITH_BOOST
	initBoostOptions();
#endif // WITH_BOOST
}

CmfdConfig::~CmfdConfig() {
}

std::string CmfdConfig::getString() const {
	std::stringstream s;
	if (prefix_enabled){
		s << "[" << prefix << "]" << std::endl;
	}
	s << "input=" << inputfile << " # Image to process" << std::endl
	  << "output=" << outputdir << " # Working directory" << std::endl
	  << "suffix=" << suffix
	  << " suffix - outputfile = img-basename + internal-suffix + suffix + ending" << std::endl
	  << "groundTruthFile=" << groundTruthFile << " # Groundtruth file" << std::endl

		 // GENERAL
	  << "numThreads=" << numThreads
	  << " # Number of threads which should be used, -1 == num of cpu cores"<< std::endl
	  << "checkEntropy=" << checkEntropy
	  << "should the entropy be checked -> blocks with too low entropy are removed"<< std::endl
	  << "loadCorresFile=" << loadCorresFile
	  << " load a serialized correspondence map from file"<< std::endl
	  << "writeCorresFile=" << writeCorresFile
	  << " write a serialized correspondence map to file"<< std::endl
		 // PREPROCESSING
	  << "chan=" << p.chan
	  << "# Possible are: ALL, GRAY, RED, GREEN, BLUE\n" << std::endl
	  << "# e.g. if you choose red only the red channel will be used for the feature computations"<< std::endl
	  << "pyramidLvl=" << p.pyramidLvl
	  << "You can set a pyramid-level, a gaussian pyramid will be build up to this level\n" << std::endl
	  << " and the highest (= level specified) will be chosed as new image"<< std::endl
	  << "waveletLvl=" << p.waveletLvl
	  << "You can set a wavelet-level, a wavelet decomposition (Haar = Daubechies D2)\n" << std::endl
	  << " will be build up to this level\n" << std::endl
	  << " and the highest (= level specified) will be chosen as new image"<< std::endl
	  << "Applies Gray-coding to the image"<< std::endl
		 // BLOCKS
	  << "blockSize=" << b.blockSize
	  << "The size a block has, has to be quadratic"<< std::endl
	  << "step=" << b.step
	  << "How many pixels should lie between two blocks" << std::endl
	  << "in x and y direction"<< std::endl
	  << "roi=";
	for (unsigned int i = 0; i < b.roi.size(); i++)
		s << b.roi[i] << " ";
	s << std::endl
	  << "You can specify region(s) of interest (ROI) in which it'll be searched" << std::endl
	  << " for similar blocks - if you don't specify a ROI, the whole image will be used\n" << std::endl
	  << "Form: x1 y1 x2 y2\n" << std::endl
	  << "the pair x1 and y1 have to be one corner of the window\n" << std::endl
	  << "the other pair has to be the opposite corner (diagonal)\n" << std::endl
	  << "you can specify more than one window, every 4 values are one region of interest"<< std::endl
		 // FEATURE
	  << "detectorType" << f.detectorType
	  << "method used for keypoint detections - for keypoint based methods only"<< std::endl
	  << "descriptorType" << f.descriptorType
	  << "method for keypoint based feature extraction"<< std::endl
	  << "normalize" << f.normalize
	  << "normalize the feature matrix"<< std::endl		  
	  << "featvec=" << f.featvec
	  << "Kind of Feature Vector\n" << std::endl
	  << "\tPossible are:\n\t  LUO, BLUR, DCT, PCT, SVD, HU, DWTFEAT, FMT, LIN, BRAVO, CIRCLE, ZERNIKE\n" << std::endl
	  << " and the keypoint based features SIFT, SURF\n" << std::endl
	  << " and NO for being the block itself the data"<< std::endl
	  << "modified=" << f.modified
	  << "If this option is set, the modified versions (if available of the feature algorithm will be used)"<< std::endl
	  << "circle=" << f.circle
	  << "if true, only a circle around the midpoint of the block will be used  for feature calculation"<< std::endl

	  << "qf=" << f.qf
	  << "Qualityfactor for DCT"<< std::endl
	  << "numHu=" << f.numHu
	  << "how many hu moments shall be used"<< std::endl
	  << "decomposeLvl=" << f.decomposeLvl
	  << "If DWTFEAT is chosen, this specifies the decomposition level"<< std::endl
	  << "degZernike=" << f.degZernike
	  << "Degree of zernike moments" << std::endl
	  << "usePCA=" << f.usePCA
	  << "If set to true: a pca of the features is done" << std::endl
	  << "dim=" << f.dim
	  << "New dimension for PCA/SVD - if set to 0, than this will computed dinamically with the eps\n"  << std::endl
	  << "if set to -1 (only SVD) then the dimension will be the average dimension of all blocks\n"  << std::endl
	  << "if set to -2 (only SVD) then the dimension is computed every time. Features with size smaller"  << std::endl
	  << " than the maximal dimension will be padded with zeros" << std::endl
	  << "eps=" << f.eps
	  << "fraction of ignored variance along the principal axis - only used if no dimension is set" << std::endl
	  << "dstwidth=" << f.dstwidth
	  << "used for FMT , size of the destination width of the log polar trafo" << std::endl
		 // QUANTIZATION
	  << "qStart=" << f.qStart
	  << "Start of quantization of the featurevector (zero based)" << std::endl
	  << "qEnd=" << f.qEnd
	  << "End of quantization of the featurevector (zero-based and not including the end)\n"  << std::endl
	  << "If is set to 0 , no quantisation is performed\n"  << std::endl
	  << "If set to -1, the whole feature will be quantized - qEnd is set to featureSize" << std::endl
	  << "qNum=" << f.qNum
	  << "Number of quantization bins in which of the featurevector is quantisized - so if set to 1, nothing changed)" << std::endl
	  << "qRound=" << f.qRound
	  << " if you quantize you can set this option to round before quantization - will be rounded down (floor)" << std::endl
	  << "qRound2=" << f.qRound2
	  << " if you quantize you can set this option to round after quantization - will be rounded down (floor)" << std::endl
	  << "qAbs=" << f.qAbs
	  << " if you quantize you can set this option so take the absolute value of the features" << std::endl

		 // CPS
	  << "useCps=" << f.useCps
	  << "uses cross power spectrum to search for duplicated regions\n"  << std::endl
	  << "you have to choos NO as featvec!" << std::endl
	  << "partsInX=" << f.partsInX
	  << "for cps: parts in x direction" << std::endl
	  << "partsInY=" << f.partsInY
	  << "for cps: parts in y direction" << std::endl
	  << "correlTh1=" << f.correlTh
	  << "for cps: correlation threshold (in feature)" << std::endl
	  << "correlTh2=" << m.correlTh
	  << "for cps: correlation threshold (in mat)" << std::endl
	  << "pixelDiff=" << f.pixelDiff
	  << "for cps: maximal difference between similar pixel" << std::endl	  

		 // MATCHING
	  << "numRows=" << m.numRows
	  << "Number of consecutive rows of the sorted blocks which will be tested further\n" << std::endl
	  << "minDistEuclidian=" << m.minDistEuclidian
	  << "minimal Euclidian distance between to block positions\n"  << std::endl
	  << "should at least be as big as the blockSize" << std::endl
	  << "minDistCheby=" << m.minDistCheby
	  << "minimal Chebychev distance between to block positions\n"  << std::endl
	  << "should at least be as big as the blockSize" << std::endl
	  << "crit=";
	for (unsigned int i = 0; i < m.crit.size(); i++)
		s << m.crit[i] << " ";
	s << std::endl;
	s << "The similarity criterion to decide which block (respectivly the featurevector of a block) is similar to another\n"  << std::endl
	  << "Possible are:\n"  << std::endl
	  << "LUOCRIT (should only be used in combination with featvec LUO),\n"  << std::endl
	  << "\t->then you have to set 9 thresholds\n"  << std::endl
	  << "EUCLIDIAN computes eculdian distance in featurespace and compares it with the threshold\n"  << std::endl
	  << "CPS computes cross power spectrum between 2 blocks and compares the max with correlTh\n"  << std::endl
	  << "CORREL computes correlation coefficient of the fourier magnitudes of the log polar trafo of a block\n"  << std::endl
	  << "BRAVOCRIT compares the first 3 avg colors with thresholds and then performs CORREL criterion\n" << std::endl
	  << "th=";
	for (unsigned int i = 0; i < m.th.size(); i++)
		s << m.th[i] << " ";
	s << std::endl
	  << "Threshold(s) to decide if two blocks are similar. "  << std::endl
	  << "This is a vector, so if a method needs more thresholds just add more." << std::endl
	  << "euclidianTh=" << m.euclidianTh
	  << "If set: the first ths will be seen as an euclidian distance\n"  << std::endl
	  << " # radius of sphere around featureposition\n"  << std::endl
	  << "otherwise it is seen as a similarity measurement in percent" << std::endl
	  << "kdsort=" << m.kdsort
	  << "if set: kd sort instead lexicographic sort will be used\n"  << std::endl
	  << "this is helpful with EUCLIDIAN as simCrit" << std::endl
	  << "nnBrute=" << m.nnBrute
	  << "nearest neighbor by brute-force" << std::endl
	  << "entropyTh=" << m.entropyTh
	  << "entropy threshold (e.g. for bravo-solorio et al.)" << std::endl
	  << "phaseCorrelTh=" << m.phaseCorrelTh
	  << "for cps: correlation threshold" << std::endl

		 //VERIFICATION=
	  << "cluster" << v.cluster
	  << "use clustering for finding and marking clusters in the x,y positions of the features" << std::endl
	  << "cutTh" << v.cutTh
	  << "cluster cut threshold" << std::endl
	  << "estimateMethod=" << v.estimateMethod
	  << " # method for estimation of the transformation matrix in the verification" << std::endl
	  << "ranReprojTh" << v.ranReprojTh
	  << "reprojection threshold for RANSAC" << std::endl
	  << "ranIterations" << v.ranIterations
	  << "number of iterations for RANSAC" << std::endl
	  << "fastsats=" << v.fastsats
	  << "use fast-SATS algorithm for verification" << std::endl
	  << "withTree=" << v.withTree
	  << "fast-SATS with kd-tree" << std::endl
	  << "checkShift=" << v.checkShift
	  << "If this is set to false it wont be checked for the same shift vector" << std::endl
	  << "minSameShift=" << v.minSameShift
	  << "How many blockpairs have to have the same shiftvector" << std::endl	
	  << "numMax=" << v.numMax
	  << "number of maximum shift vectors" << std::endl
	  << "shiftVariance" << v.shiftVariance
	  << "variance between shift vectors" << std::endl	
	  << "maxDist=" << v.maxDist
	  << "maxDistance maximum distance of pixels from the analyzed one" << std::endl	  
	  << "recomputationTh=" << v.recomputationTh
	  << " recompute transformation matrix when recomputationTh points are added" << std::endl
		 // MARKING
	  << "scaleUp=" << mk.scaleUp
	  << "if img is downscaled, e.g. due to gaussian pyramid, this scales it up again" << std::endl
	  << "useOrig=" << mk.useOrig
	  << "use original image instead of modified e.g. due to grey conversion" << std::endl	
	  << "minVal=" << mk.minVal
	  << "only pixels greater than this value will be marked" << std::endl
	  << "markShift=" << mk.markShift
	  << " # marks shift vector in image"
	  << "writeChrom=" << mk.writeChrom
	  << " # writes an image with chromaticity = block positions" << std::endl
	  << "writePairs=" << mk.writePairs
	  << " # marks block pairs, visualizes pairs in GT differently" << std::endl
	  << "writeMatrix=" << mk.writeMatrix
	  << " # writes down the matrix, a value != 0 in the matrix indicates a match" << std::endl
	  << "writeOutput=" << mk.writeOutput
	  << " # writes an output image with the duplicated regions marked" << std::endl
	  << "writeCluster=" << mk.writeCluster
	  << " # writes a matrix, every cluster has a different color" << std::endl
	  << post.getString();
	return s.str();
}

#ifdef VOLE_GUI
QWidget *CmfdConfig::getConfigWidget() {
	// create a qt widget here, fill the widget of the parent class
	this->initConfigWidget();
	QVBoxLayout *cmfd_config = new QVBoxLayout();
	// ...
	layout->addLayout(cmfd_config);
	configWidget->setLayout(layout);
	return configWidget;
}

void CmfdConfig::updateValuesFromWidget() {
	// pull values out of the widget defined in getConfigWidget()
}
#endif //VOLE_GUI

#ifdef WITH_BOOST
void CmfdConfig::initBoostOptions() {
	options.add(post.options);
	options.add_options()
			("input,I", value(&inputfile)->default_value(""),
			 "Image to process")
			("output,O", value(&outputdir)->default_value("/tmp/"),
			 "Working directory")
			("log-times", bool_switch(&log_times)->default_value(false),
			 "create a log file where times are protocolled")
			("log-file", value(&log_file)->default_value(""),
			 "log to a specific file")
			("graphical", bool_switch(&graphical)->default_value(false),
			 "Show any graphical output during runtime")
			("suffix", value(&suffix)->default_value(""),
			 " suffix - outputfile = img-basename + internal-suffix + suffix + ending")
			("groundTruthFile",value(&groundTruthFile)->default_value(""),
			 "path to the ground truth image")
			// GENERAL
			("numThreads", value(&numThreads)->default_value(-1),
			 "Number of threads which should be used, -1 == num of cpu cores")
			("checkEntropy", bool_switch(&checkEntropy)->default_value(false),
			 "should the entropy be checked -> blocks with too low entropy are removed\n"
			 "(for bravo solorio plz use just an entropyTh > 0.0 and don't use this optioon)")
			("loadCorresFile", value(&loadCorresFile)->default_value(""),
			 "load a serialized file of the correspondence map")
			("writeCorresFile", bool_switch(&writeCorresFile)->default_value(false),
			 "write a serialized file of the correspondence map - may be reused")
			// PREPROCESSING
			("chan", value(&p.chan)->default_value("GRAY"),
			 "Possible are: ALL, GRAY, RED, GREEN, BLUE\n"
			 "e.g. if you choose red only the red channel will be used for the feature computations")
			("pyramidLvl", value(&p.pyramidLvl)->default_value(0),
			 "You can set a pyramid-level, a gaussian pyramid will be build up to this level\n"
			 " and the highest (= level specified) will be chosed as new image")
			("waveletLvl", value(&p.waveletLvl)->default_value(0),
			 "You can set a wavelet-level, a wavelet decomposition (Haar = Daubechies D2)\n"
			 " will be build up to this level\n"
			 " and the highest (= level specified) will be chosen as new image")
			// BLOCKS
			("blockSize,B", value(&b.blockSize)->default_value(16),
			 "The size a block has, has to be quadratic")
			("step", value(&b.step)->default_value(1),
			 "How many pixels should lie between two blocks"
			 "in x and y direction")
			("roi", value(&b.roi),
			 "You can specify region(s) of interest (ROI) in which it'll be searched"
			 " for similar blocks - if you don't specify a ROI, the whole image will be used\n"
			 "Form: x1 y1 x2 y2\n"
			 "the pair x1 and y1 have to be one corner of the window\n"
			 "the other pair has to be the opposite corner (diagonal)\n"
			 "you can specify more than one window, every 4 values are one region of interest")
			// FEATURE			
			("featvec", value(&f.featvec)->default_value(NO),
			 "Kind of Feature Vector\n"
			 "\tPossible are:\n\t  LUO, BLUR, DCT, PCT, SVD, HU, DWTFEAT, FMT, LIN,\n"
			 " BRAVO, CIRCLE, ZERNIKE, NO(=block itself)\n"			
			 " and NO for being the block itself the data")
			("detectorType", value(&f.detectorType)->default_value(""),
			 "type of detector - keypoint extraction - only for keypoint based methods")
			("descriptorType", value(&f.descriptorType)->default_value(""),
			 "type of descriptor - feature - only for keypoint based methods")
			("normalize", bool_switch(&f.normalize)->default_value(false),
			 "normalizes the feature-vectors between 0 and 1")
			("maxFeatureSize", value(&f.maxFeatureSize)->default_value(512),
			 "maximal feature size if feature is smaller or if maxFeatureSize==0\n"
			 "it will be ignored internally and the storage won't change,\n"
			 "else the featur-matrix will be reduced to col = maxFeatureSize\n"
			 "Note: only for block-based features")
			("minFeatureSize", value(&f.minFeatureSize)->default_value(0),
			 "can be used for (K)PCA")
			("modified", bool_switch(&f.modified)->default_value(false),
			 "If this option is set, the modified versions (if available of the feature algorithm will be used")
			("circle", bool_switch(&f.circle)->default_value(false),
			 "if true, only a circle around the midpoint of the block will be used \n"
			 "for feature calculation")
			("qf", value(&f.qf)->default_value(0.5),
			 "Qualityfactor for DCT")
			("numHu", value(&f.numHu)->default_value(4),
			 "how many hu moments shall be used")
			("decomposeLvl", value(&f.decomposeLvl)->default_value(0),
			 "If DWTFEAT is chosen, this specifies the decomposition level")
			("degZernike", value(&f.degZernike)->default_value(5),
			 "degree of Zernike Moments")
			("dstwidth", value(&f.dstwidth)->default_value(6),
			 "used for FMT , size of the destination width of the log polar trafo")
			// (K)PCA
			("usePCA", bool_switch(&f.usePCA)->default_value(false),
			 "If set to true: a pca of the features is done")
			("dim", value(&f.dim)->default_value(0),
			 "New dimension for PCA/SVD - if set to 0, than this will computed dinamically with the eps\n"
			 "if set to -1 (only SVD) then the dimension will be the average dimension of all blocks\n"
			 "if set to -2 (only SVD) then the dimension is computed every time. Features with size smaller"
			 " than the maximal dimension will be padded with zeros")
			("eps", value(&f.eps)->default_value(0.01),
			 "fraction of ignored variance along the principal axis - only used if no dimension is set")
			("sigma", value(&f.sigma)->default_value(25.0),
			 "sigma for kernel-varianz of KPCA")
			("numSamples", value(&f.numSamples)->default_value(144),
			 "number of training samples of KPCA")
			// QUANTIZATION
			("qStart", value(&f.qStart)->default_value(0),
			 "Start of quantization of the featurevector (zero based)")
			("qEnd", value(&f.qEnd)->default_value(0),
			 "End of quantization of the featurevector (zero-based and not including the end)\n"
			 "If is set to 0 , no quantisation is performed\n"
			 "If set to -1, the whole feature will be quantized - qEnd is set to featureSize")
			("qNum", value(&f.qNum)->default_value(1),
			 "Number of quantization bins in which of the featurevector is quantisized - so if set to 1, nothing changed")
			("qRound", bool_switch(&f.qRound)->default_value(false),
			 " if you quantize you can set this option to round before quantization - will be rounded down (floor)")
			("qRound2", bool_switch(&f.qRound2)->default_value(false),
			 " if you quantize you can set this option to round after quantization - will be rounded down (floor)")
			("qAbs", bool_switch(&f.qAbs)->default_value(false),
			 " if you quantize you can set this option so take the absolute value of the features")

			// CPS
			("useCps", bool_switch(&f.useCps)->default_value(false),
			 "uses cross power spectrum to search for duplicated regions\n"
			 "you have to choos NO as featvec!")
			("partsInX", value(&f.partsInX)->default_value(2),
			 "for cps: parts in x direction")
			("partsInY", value(&f.partsInY)->default_value(2),
			 "for cps: parts in y direction")
			("correlTh1", value(&f.correlTh)->default_value(0.025, "0.025"),
			 "for cps: correlation threshold (in feature)")
			("correlTh2", value(&m.correlTh)->default_value(0.9, "0.9"),
			 "for cps: correlation threshold (in mat)")
			("pixelDiff", value(&f.pixelDiff)->default_value(2),
			 "for cps: maximal difference between similar pixel")			

			// MATCHING
			("numRows", value(&m.numRows)->default_value(1),
			 "Number of consecutive rows of the sorted blocks which will be tested further\n")
			("minDistEuclidian", value(&m.minDistEuclidian)->default_value(50),
			 "minimal Euclidian distance between to block positions\n"
			 "should at least be as big as the blockSize")
			("minDistCheby", value(&m.minDistCheby)->default_value(0),
			 "minimal Chebychev distance between to block positions\n"
			 "should at least be as big as the blockSize")
			("crit", value(&m.crit),
			 "The similarity criterion to decide which block (respectivly the featurevector of a block) is similar to another\n"
			 "Possible are:\n"
			 "LUOCRIT (should only be used in combination with featvec LUO),\n"
			 "\t->then you have to set 9 thresholds\n"
			 "EUCLIDIAN computes eculdian distance in featurespace and compares it with the threshold\n"
			 "CPS computes cross power spectrum between 2 blocks and compares the max with correlTh\n"
			 "CORREL computes correlation coefficient of the fourier magnitudes of the log polar trafo of a block\n"
			 "BRAVOCRIT compares the first 3 avg colors with thresholds and then performs CORREL criterion\n"
			 "EUCLIDIAN_RATIO compares the ratio of the distance between the i and the i+1 feature vector"
			 )
			("th", value(&m.th),
			 "Threshold(s) to decide if two blocks are similar. "
			 "This is a vector, so if a method needs more thresholds just add more.")
			("euclidianTh", value(&m.euclidianTh)->default_value(true),
			 "If set: the first ths will be seen as an euclidian distance\n"
			 " = radius of sphere around featureposition\n"
			 "otherwise it is seen as a similarity measurement in percent")
			("kdsort", bool_switch(&m.kdsort)->default_value(false),
			 "if set: kd sort instead lexicographic sort will be used\n"
			 "this is helpful with EUCLIDIAN as simCrit")
			("nnBrute", bool_switch(&m.nnBrute)->default_value(false),
			 "compute nearest neighbor by brute-force")
			("entropyTh", value(&m.entropyTh)->default_value(0.0),
			 "entropy threshold for the method of bravo-solorio et al.")
			("phaseCorrelTh", value(&m.phaseCorrelTh)->default_value(0.9, "0.9"),
			 "for cps: correlation threshold")

			//VERIFICATION
			("cluster", bool_switch(&v.cluster)->default_value(false),
			 "use clustering to find and mark clusters of x,y coordinates of features")
			("cutTh", value(&v.cutTh)->default_value(0.0),
			 "cluster cut threshold")
			("normalizeCluster", bool_switch(&v.normalizeCluster)->default_value(false),
			 "normalize the cluster points such that the points are translated "
			 " to the cluster-center and scaled that avg-distance = sqrt(2)")
			("estimateMethod", value(&v.estimateMethod)->default_value(""),
			 "method to find a valid transformation between copied and pasted regions")
			("ranReprojTh", value(&v.ranReprojTh)->default_value(3.0),
			 "reprojection threshold for ransac method")
			("ranIterations", value(&v.ranIterations)->default_value(100),
			 "number of iterations used for RANSAC")
			("fastsats", bool_switch(&v.fastsats)->default_value(false),
			 "fast same affine transformation selection on shift vectors")
			("withTree", bool_switch(&v.withTree)->default_value(false),
			 "fast SATS with kd-tree")
			("checkShift", bool_switch(&v.checkShift)->default_value(false),
			 "If this is set to false it wont be checked for the same shift vector")
			("minSameShift", value(&v.minSameShift)->default_value(50),
			 "How many blockpairs have to have the same shiftvector")			
			("numMax", value(&v.numMax)->default_value(0),
			 "number of maximum shift vectors")
			("shiftVariance", value(&v.shiftVariance)->default_value(0),
			 "variance between shift vectors")

			("maxDist", value(&v.maxDist)->default_value(5),
			 "maxDistance maximum distance of pixels from the analyzed one")		
			("recomputationTh", value(&v.recomputationTh)->default_value(2),
			 "recompute transformation matrix when recomputationTh points are added")
			("computeCorrelMap", bool_switch(&mk.computeCorrelMap)->default_value(false),
			 "compute a correlation map with the help of earlier computed transformation matrices")
			// MARKING
			("writeChrom", bool_switch(&mk.writeChrom)->default_value(false),
			 "writes an image with chromaticity = block positions")
			("writePairs", bool_switch(&mk.writePairs)->default_value(false),
			 "marks block pairs, visualizes pairs in GT differently")
			("writeOutput", bool_switch(&mk.writeOutput)->default_value(false),
			 "writes an ouput image with the duplicated regions marked")
			("writeMatrix", bool_switch(&mk.writeMatrix)->default_value(false),
			 "writes down the matrix, a value != 0 in the matrix indicates a match")
			("writeCluster", bool_switch(&mk.writeCluster)->default_value(false),
			 "writes a matrix, every cluster has a different color")
			("writeCorrelation", bool_switch(&mk.writeCorrelation)->default_value(false),
			 "write a binarized correlation map, transformation estimation needed")
			("writePost", bool_switch(&mk.writePost)->default_value(false),
			 "write the post-processed output matrix")
			("scaleUp", bool_switch(&mk.scaleUp)->default_value(false),
			 "if img is downscaled, e.g. due to gaussian pyramid, this scales it up again")
			("useOrig", bool_switch(&mk.useOrig)->default_value(false),
			 "use original image instead of modified e.g. due to grey conversion")			
			("minVal", value(&mk.minVal)->default_value(0),
			 "only pixels greater than this value will be marked")
			("markRegions", bool_switch(&mk.markRegions)->default_value(false),
			 "marks forged regions each in different colors in image")
			("markShift", bool_switch(&mk.markShift)->default_value(false),
			 "marks shift vector in image")
			("markContour", bool_switch(&mk.markContour)->default_value(false),
			 "draw the countour of the areas in the original image")
			("binaryTh", value(&mk.binaryTh)->default_value(0.3, "0.3"),
			 "threshold for binarization of the correlation maps")
			;
}
#endif // WITH_BOOST

