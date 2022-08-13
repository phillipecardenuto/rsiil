/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef COPYMOVEFRAMEWORK_CONFIG_H
#define COPYMOVEFRAMEWORK_CONFIG_H

#include <streambuf>

#include "vole_config.h"
#include "post_process/post_process_config.h"

#ifdef WITH_BOOST
using namespace boost::program_options;
#endif // WITH_BOOST

/// Type of Feature
enum featureType {
	LUO,
	BLUR,
	DCT,
	PCT,
	SVD,
	HU,
	DWTFEAT,
	FMT,
	LIN,
	BRAVO,
	CIRCLE,
	ZERNIKE,
	KPCA,
	NO,
	NOTHING
};
/// the Criterion which is used for the similarity map
enum similarityCriterion {
	NONE = 0,
    LUOCRIT = 1,
    EUCLIDIAN = 2,
	CPS = 4,
	CORREL = 8,
	BRAVOCRIT = 16,
	EUCLIDIAN_RATIO = 32
};

// to use this enum in configuration, following has to be defined
// see also ENUM_MAGIC in ?.cxx file
#define featureTypeString {"LUO", "BLUR", "DCT", "PCT", "SVD",\
	"HU", "DWTFEAT", "FMT", "LIN", "BRAVO", "CIRCLE",\
	"ZERNIKE", "KPCA", "NO", "NOTHING"}
#define similarityCriterionString {"NONE", "LUOCRIT", "EUCLIDIAN", "CPS", \
	"CORREL", "BRAVOCRIT", "EUCLIDIAN_RATIO"}

// Preprocessing
struct pr_cfg {
	/// which channel should be used
	std::string chan;
	/// pyramid level
	int pyramidLvl;
	/// wavelet level
	int waveletLvl;
};

// BLOCKS
struct bl_cfg {
	/// Size of a block
	unsigned int blockSize;
	/// The stepsize between two blocks (x and y direction)
	int step;
	/** coordinates for region of interest:
	 * two coordinate-pairs (x and y) has to be given for one window
	 * -> so 4 values for one window
	 * every coordinate pair has to be a corner of the same diagonal of the window
	 */
	std::vector<unsigned int> roi;
};

// FEATURE
struct ft_cfg {	
	/// Type of feature vector
	featureType featvec;
	/// Detector type of keypoint based features 
	std::string detectorType;
	/// Descriptor type of keypoint based features
	std::string descriptorType;

	/// normalize feature such that L2-Norm == 1
	bool normalize;
	/// maximal feature size if feature is smaller or if maxFeatureSize==0
	/// it will be ignored internally the storage won't change,
	/// just the feature-head will be changed (-> not continuous any more)
	/// Note: only for block-based features atm
	int maxFeatureSize;
	/// can be used for (K)PCA
	int minFeatureSize;
	/// if set, the modified feature calculation (if available) will be used
	/// some features are implemented in a modified version which can be
	/// enabled this way
	bool modified;
	/// if true, only a circle around the midpoint of the block will be used
	/// for feature calculation
	bool circle;
	/// quality factor for DCT
	double qf;
	/// how many hu-moments shall be computed
	int numHu;
	/// decomposeLvl for dwt features
	int decomposeLvl;
	/// degree of zernike moments
	int degZernike;

	/// uses the pca after a feature calculation
	bool usePCA;
	/// new dimension for PCA/SVD
	int dim;
	/// eps - if no dimension set for PCA/SVD
	double eps;
	/// for FMT
	int dstwidth;
	/// sigma for kernel-varianz of KPCA
	double sigma;
	/// number of training samples of KPCA
	int numSamples;

	// QUANTIZATION of the feature vector
	/// start of quantization (zero-based)
	int qStart; 
	/// end  of quantization (zero-based)
	int qEnd;
	/// Number of quantization bins
	double qNum;
	/// round down before quantization
	bool qRound;
	/// round down after the quantization
	bool qRound2;
	/// take the absolute value
	bool qAbs;
	
	/// use cross power spectrum for phase correlation
	bool useCps;
	/// parts for Cps
	int partsInX;
	int partsInY;
	/// correlation threshold for Correlation coefficient
	double correlTh;
	/// pixel Difference threshold for Cps
	int pixelDiff;
};

// MATCHING
struct mat_cfg{
	/// number of consecutive rows or if 'kdsort' is true these are the
	/// number of (approxiamated) nearest neighbors in the kd-tree
	int numRows;
	/// the minimal euclidian distance between two block positions
	int minDistEuclidian;
	/// the minimal Chebychv distance between two block positions
	int minDistCheby;
	/// similarity criterion to decide which blocks are similar to each other
	std::vector<similarityCriterion> crit;
	/// vector of thresholds for checking similarity
	std::vector<double> th;
	/// interpretation of the th for euclidian distance -> true: its given as a radius of the sphere
	bool euclidianTh;
	/// kdtree for sorting instead lexicographic sorting
	bool kdsort;
	/// nearest neighbor w brute-force
	bool nnBrute;
	/// correlation threshold for Correlation coefficient
	/// (this is a copy of correlTh of struct ft_cfg)
	double phaseCorrelTh;
	/// correlation threshold for the similarity criterion CORREL
	double correlTh;
	/// entropy threshold (e.g.for BRAVO)
	double entropyTh;
};

// VERIFICATION
struct very_cfg {
	/// builds cluster of the feature coordinates
	bool cluster;
	/// cluster cut threshold
	double cutTh;
	/** normalize the cluster points such that the points are translated
		to the cluster-center and scaled that avg-distance = sqrt(2)
	 */
	bool normalizeCluster;
	/// method for estimation the transformation matrix:
	/// ransac, ownransac, lmeds, pure
	std::string estimateMethod;
	/// reprojection threshold for ransac
	double ranReprojTh;
	/// number of iterations for RANSAC
	double ranIterations;
	/// using fast same affine transformation selection
	bool fastsats;
	/// fast sats with the tree-variant
	bool withTree;
	/// if true then the shift vectors will be checked
	bool checkShift;
	/// how many block have to have the same shiftvector
	unsigned int minSameShift;
	/// number of maximum shift vectors;
	int numMax;
	/// variance between shift vectors
	int shiftVariance;

	/// maxDistance maximum distance of pixels from the analyzed one
	int maxDist;	
	/// recompute transformation matrix when recomputationTh points are added
	int recomputationTh;
};

// MARKING
struct mark_cfg {
	/// computes a correlation map
	/// Note: this is kinda in-between of verification and marking
	bool computeCorrelMap;
	/// if img is downscaled, e.g. due to gaussian pyramid, this
	/// scales it up again
	bool scaleUp;
	/// use original image instead of modified e.g. due to gray conversion
	bool useOrig;	
	/// only pixels greater than this value will be marked
	int minVal;
	/// threshold for binarization of the correlation maps
	double binaryTh;
	/// mark the regions w. different colors
	bool markRegions;
	/// mark the shift vectors
	bool markShift;
	/// draws the contours of the detected regions in the original image
	bool markContour;
	/// writes an image with chromaticity = block positions
	bool writeChrom;
	/// writes block pairs - pairs in ground truth marked differently
	bool writePairs;
	/// writes an output image with the duplicated regions marked
	bool writeOutput;
	/// makes a dump of the map
	bool writeMatrix;
	/// makes a colorful dump of the clusters
	bool writeCluster;
	/// write a binarized correlation map, transformation estimation needed
	bool writeCorrelation;
	/// write post-processed matrix
	bool writePost;
};

class CmfdConfig : public vole::Config {
public:

	CmfdConfig(const std::string & prefix = std::string());
	~CmfdConfig();

	/// graphical output
	bool graphical;
	/// input image (should be 8bit 3 or 1 channel)
	std::string inputfile;
	/// output directory
	std::string outputdir;
	/// outputfile = img-basename + internal-suffix (write-command dependend) +
	/// suffix + file-ending (write-command dependend)
	/// Note: atm only used for the end-files _matrix & output not for the serialized files
	std::string suffix;
	/// path to the ground truth image
	std::string groundTruthFile;
	/// shell the times of each step be logged?
	bool log_times;
	/// specific log-file
	std::string log_file;

	/// number of threads
	int numThreads;
	/// should the entropy be checked?
	bool checkEntropy;

	/// load correspondence map(s) from the specified file
	std::string loadCorresFile;
	/// file to write the correspondence map
	bool writeCorresFile;
	
	/// cfg for preprocessing
	struct pr_cfg p;
	/// cfg for the blockhandling
	struct bl_cfg b;
	/// cfg for features
	struct ft_cfg f;
	/// cfg for matching
	struct mat_cfg m;
	/// cfg for verification
	struct very_cfg v;
	/// cfg for marking
	struct mark_cfg mk;

	/// postprocessing config
	PPConfig post;

	virtual std::string getString() const;

    #ifdef VOLE_GUI
        virtual QWidget *getConfigWidget();
        virtual void updateValuesFromWidget();
    #endif// VOLE_GUI

protected:
    #ifdef WITH_BOOST
        virtual void initBoostOptions();
    #endif // WITH_BOOST

    #ifdef VOLE_GUI
//        QLineEdit *edit_sigma, *edit_k_threshold, *edit_min_size;
 //       QCheckBox *chk_chroma_img;
    #endif // VOLE_GUI

}; 

#endif // COPYMOVEFRAMEWORK_CONFIG_H
