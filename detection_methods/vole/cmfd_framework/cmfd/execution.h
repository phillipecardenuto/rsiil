/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef EXECUTION_H
#define EXECUTION_H

#include "copymoveframework_config.h"

#include <opencv2/nonfree/features2d.hpp>
#include <opencv2/core/core.hpp>

#include <iosfwd>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

// Forward declarations:
class BlockHandling;
class Feature;
class Matching;
class Verification;
class Mark;
class Cps;

// make things easier for other classes
#define vout(x) (Execution::verbosity_out(x))

/*! this class encapsulates functions for execution and modularize the
 *	general execute() function
 */
class Execution{
public:
	Execution(CmfdConfig &_config);
	~Execution();

	/// returns a reference to the actual image with all its manipulations
	inline const cv::Mat& getCurImage(void) const { return img; }

	/// returns a reference to the original loaded image, could be distorted due to marking step
	inline const cv::Mat& getImage(void) const { return imgorig; }

	/// get the actual config
	inline CmfdConfig getConfig(void) const { return config; }

	/// sets the config
	/*! Note: should only be called in between of execute functions
   */
	// removed due to problems in the boost copy constructor of the "options_description options" field in Config
	/*
  inline void setConfig(CmfdConfig _config) { config = _config; }
*/

	/*! loads the image, specified in global_config
   *	has to be called before you can work with the image
   */
	void loadImage(void);

	/*! performs a preprocessing step
   */
	void preprocess(void);

	/*! this method generates the overlapping block
   *	has to be done at first before you can calculate features,
   *	after that you can do different manipulations,
   *	lik calculating features and match the blocks with the
   *	similarity functions
   */
	void generateBlocks(void);

	/// computes the features for every block
	void computeFeatures(void);

	/*! computes similare blockpairs with respect to our similarity
   *	criterion which was chosen
   */
	void computeMatching(void);

	/*! verifies the similar block pairs -> filters them
   *	can only be used after computeSimilarity
   *	\param corres vector of correspondence matrices for every channel
   *	       of the image you get a correspondence matrix
   */
	void verifySimilarity(void);

	/// compute correlation map as proposed in Pans paper
	cv::Mat_<float> computeCorrelMap(const cv::Mat_<uchar> &plane, const cv::Mat &trans);

	/// marks similar block pairs to visualize it better
	void markIt(void);

	/// static verbosity level 0-3 s. core/common/config.h just for simpler usage
	static int verbosity;

	/// static outputdir for simpler usage
	static std::string outputdir;
	/// static image basename for simpler usage
	static std::string image_basename;
	/// static suffix for simpler usage
	static std::string suffix;

	/// converts the correspondence matrix to unsigned short matrix
	/// which can then be used to be printed as an image
	std::vector<cv::Mat> corresToImg(void) const;
	/// the other way around of corrresToImg (converts
	/// 16bit image to correspondence matrix)
	void imgToCorres(const std::vector<cv::Mat> & corres_in_img);

	/// ease the output
	static std::ostream & verbosity_out(int x){
		class null_out_stream : public std::ostream
		{
		public:
			null_out_stream() : std::ios(0), std::ostream(0){}
		};

		static null_out_stream cnul;
		return (verbosity >= x) ? std::cout : cnul;
	}
	/// return the number of blocks/keypoints processed
	int getBlockCount() const;
	int getMatchCount() const;
private:
	friend class boost::serialization::access;
	// for serialization
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & dimx;
		ar & dimy;
		ar & config.b.blockSize;
		ar & config.m.numRows;
		ar & config.p.chan;
		ar & image_basename;
		// TODO: sometime, not now in the middle of our experiments! -> compatibility!		
		// TODO: also the suffix!
		// ar & config.suffix
	}
	/// updates the time and returns relative time
	inline double relativeTime(){
		double t = static_cast<double>(cv::getTickCount());
		double elapsed = (t - old_time) / cv::getTickFrequency();
		old_time = t;
		return elapsed;
	}
	/// local copy of config
	CmfdConfig config;
	/// working copy of img, may be altered due to preprocessing etc
	cv::Mat img;
	/// dimension in x direction (= width) of img
	int dimx;
	/// dimension in y direction (= height) of img
	int dimy;
	/// copy of the original image, will only be modified in the marking process
	cv::Mat_<cv::Vec3b> imgorig;
	/// reference to the blocks, contain perhaps features (depends at the actual step of our chain)
	BlockHandling *imgBlocks;
	/// feature matrix of keypoints
	cv::Mat feature_matrix;
	/// the actual feature class (i.e. DCT, LUO, BLUR,..)
	Feature * feature;
	cv::Ptr<cv::FeatureDetector> feature_detector;
	cv::Ptr<cv::DescriptorExtractor> descriptor_extractor;
	std::vector<cv::KeyPoint> keypoints;
	/// instead of features use cross power spectrum
	Cps *cps;
	/// correspondence matrix
	std::vector<cv::Mat_<cv::Point> > corres_matrix;
	/// Matching instance
	Matching *match;
	/// Verification instance;
	Verification *veri;
	/// number of threads used
	unsigned int threadnum;
	/// starting time stamp
	double time;
	/// relative time (between two time stamps)
	double old_time;

	/// the matrix to be marked in the marking step
	cv::Mat_<uchar> markMatrix;
	/// individual clusters
	std::vector<std::vector<std::pair<cv::Point,cv::Point> > > point_groups;
	/// transformation matrices
	std::vector<cv::Mat> transformations;
};

#endif
