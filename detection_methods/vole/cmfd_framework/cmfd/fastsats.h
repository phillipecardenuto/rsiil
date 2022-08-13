/*
   Copyright(c) 2010 Vincent Christlein <vincent.christlein@web.de>
 and Christian Riesse <christian.riess@cs.fau.de>

   This file may be licensed under the terms of of the GNU General Public
   License, version 3, as published by the Free Software Foundation. You can
   find it here: http://www.gnu.org/licenses/gpl.html

 If you use this code in your research, please cite:
 V.Christlein and C.Riess and E.Angelopoulou."On Rotation Invariance in Copy-Move Forgery Detection",
 Workshop on Information Forensics and Security, Seattle, December 2010

 Suggestions for improvements or other feedback are very welcome!
*/

#ifndef FASTSATS_H
#define FASTSATS_H

#include <vector>
#include <deque>
#include <stack>
#include <queue>

// TODO: change this to opencv 2.2 header convention if required, i.e. opencv2/imgproc/imgproc.hpp
#include <cv.h>

class Fastsats {
public:
	/*! Constructor
	 * \param corres point correspondences between coordinates of blocks.
	 *   The format is as follows: Assume that e.g. for every image pixel,
	 *   one feature vector is extracted from the neighborhood of this
	 *   pixel (i.e. "one block"). This feature vector is matched to its
	 *   closest neighbor in feature space. The image coordinates of the
	 *   closest neighbor of the feature vector from point (x,y) are stored
	 *   in corres[y][x], otherwise (-1,-1).
	 * \param blockSize size of a block
	 * \param verbosity verbose mode for debugging (ranging from 0 to 3, 0 =
	 *  no output at all, 1 = some output, 2 = additional output
	 * \param maxDist Upper limit how far (in pixels) the matches of
	 *   neighbored blocks may be away from each other. Implicitly, this
	 *   imposes an upper bound on the detection of scaled CMFD
	 * \param step Step size over the image pixels; useful to omit
	 * 	 blocks, e.g. with a step size of 2, the number of considered block
	 * 	 pairs is reduced by a factor of 4 (4 = 2*2, in x- and
	 * 	 y-direction).
	 * \param color [TODO: Not in use atm] If true we have for every individual
	 *   color-channel 1-n
	 *   correspondence-matrices, in the moment every channel is handled
	 *   individually, too. It could make sense to merge them at this stage
	 * \param with_tree Indicates whether a kd-tree should be used to get
	 *   the nearest neighbors around the current point, otherwise this will
	 *   be determined local - is VERY SLOW! if many points have correspondences
	 *   but may be fast(er?) with keypoint based features
	 */
	Fastsats(std::vector<cv::Mat_<cv::Point> > &corres,
			 int block_size = 16,
			 int verbosity = 1,
			 int max_dist = 7,
			 int step = 1,
			 bool color = false,
			 bool with_tree = false);

	/// Destructor
	~Fastsats(void);

	/*! Perform clustering of matches that adhere to the same affine
	 * transformation. Discard small groups of matches. Create a binary map
	 * of copy-moved regions.
	 *
	 * \param minSameTransformPairsTh Minimum number of block pairs with
	 *   the same linear transformation; this is equivalent to the minimum
	 *   number of shift vectors if only translation is considered, see the
	 *   paper for additional information. A set of block pairs that is
	 *   larger than this threshold is considered as copy-moved.
	 * \param recomputationTh If the number of block pairs for one
	 *   transform hypothesis exceed this threshold, the transformation
	 *   matrix is recomputed (for better stability). This implementation
	 *   differs from the paper, the threshold grows exponentially with
	 *   increasingly accurate estimations (Thus, the more pairs we have,
	 *   the more rarely the transformation is recomputed).
	 * \return a binary map of copy-moved points. A pixel is either 0 or
	 *   255 (0 = presumably not copied, 255 = presumably copied)
	 */
	cv::Mat_<uchar> sameAffineTransformationSelection(int minSameTransformPairsTh = 100,
													  int recomputationTh = 2);
	// --- additional getter-methods ---
	/*!
	 * \return vector of all transformation matrices which were used
	 *   to group points together
	 */
	inline std::vector<cv::Mat> & getTrafoMatrices(void){
		return trans;
	}
	/*!
	 * \return all_hypo_points = point-clouds / clusters
	 */
	inline std::vector<std::vector< std::pair<cv::Point,cv::Point> > > & getPointGroups(void){
		return all_hypo_points;
	}

	// --- useful functions which could be needed by other programs ---
	/// indicates the relative position of 3 points: are they clockwise
	/// rotated, counter-clockwise rotated or colinear?
	static int clockwise(const cv::Point2f idx1, const cv::Point2f idx2, const cv::Point2f idx3);
	/// applies affine transformation to a given point
	static cv::Point2f transformPoint(const cv::Point2f &point, const cv::Mat &transformation);

private:

	/// coordinate correspondence map
	std::vector<cv::Mat_<cv::Point> >  corres_map;
	/// the output map as a member variable
	cv::Mat_<uchar> output;
	/// size of a block (typically 16)
	int block_size;
	/// dimensions in x and y (= width & height)
	int dimx, dimy;
	/// verbosity mode (for debugging purposes)
	int verbosity;
	/// maximum distance between estimated point and correspondent point
	int max_dist;
	/// step size between two pixels / neighboring points
	int step;
	/// do we have individual color channels
	bool color;
	/// all transformation matrices found
	std::vector<cv::Mat> trans;
	/// all marked hypothesis points
	std::vector<std::vector<std::pair<cv::Point,cv::Point> > > all_hypo_points;

	//--- flann specific variables ---
	/// flann-tree yes or no
	bool with_tree;
	/// flann-tree of x,y-positions
	cv::flann::Index *tree;
	/// indices of neighbors
	cv::Mat_<int> indices;
	/// distances of neighbors
	cv::Mat_<float> dists;
	/// the actual data matrix (Nx2)
	cv::Mat_<float> data;
	/// the data_array of the data-matrix
	float* data_array;
	/// builds the flann-tree
	void buildTree(void);

	/// run sats normal
	void sats(int minSameTransformPairsTh, int recomputationTh);
	/// runs sats with the flann-tree as neighboring search
	void satsWithTree(int minSameTransformPairsTh, int recomputationTh);

	/// find three suitable points
	bool findThreePoints(int x, int y, cv::Point2f *src,
						 cv::Point2f *dst);

	/// adds the nearest neighbors to the todo-queue
	void addNeighbors(int point_index, std::queue<int> &todo) const;
	/// adds the nearest neighbors searched in the kd-tree to the todo-queue
	void addNeighbors(const cv::Point2f &p, std::queue<cv::Point2f> &todo) const;

	// returns the location (in image coordinates) of the closest block in feature space
	/// get correspondence of point x,y and saves it in result, returns true if it exists
	bool getCorres(int x, int y, int channel, cv::Point2f & result ) const;

	/** computes transformation matrix from a number of pairs; prefers
		spatially extremal points for numerically more stable calculations
	*/
	cv::Mat recomputeTransformation( const cv::Point2f *extremalPoints,
									 const cv::Point2f *extremalPoints_corres,
									 std::vector<cv::Point2f> &hypothesis_points) const;

	/** patches the set of extremal points, such that no three colinear
		points are chosen (points may not be colinear if the transform
		should be computed)
	*/
	void updateExtremalPoints(cv::Point2f *extremalPoints,
							  cv::Point2f *extremalPoints_corres,
							  cv::Point2f newPoint,
							  cv::Point2f newPoint_corres) const;

	/// checks if two points are within maxDist
	bool wrongDist(const cv::Point2f &p1, const cv::Point2f &p2, double max_dist) const;
	/// test the points in the test_points vector if it has a correspondence and if hasn't been already looked up
	bool testPoints(const cv::Point2f & base,
					const cv::Point2f & corres,
					const std::vector<cv::Point2f> & test_points,
					cv::Point2f & result_point,
					cv::Point2f & result_corres) const;

	// --- Helper Functions ---
	struct smaller{
		bool operator() (const cv::Point2f & a, const cv::Point2f & b) const{
			return (a.x < b.x || (a.x == b.x && a.y < b.y)) ? true : false;
		}
	} small;
	struct similar{
		bool operator() (const cv::Point2f & a, const cv::Point2f & b) const {
			return a.x == b.x && a.y == b.y;
		}
	} same;
};
#endif

