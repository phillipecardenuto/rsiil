/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef VERIFICATION_H
#define VERIFICATION_H

#include <map>
#include <vector>

class Verification {
public:
	/*! Constructor
   * \param width the width of our image
   * \param height the height of our image
   * \param decompLvl decomposition Lvl
	*/
	Verification(int blockSize,
				 int width,
				 int height,
				 std::vector<cv::Mat_<cv::Point> > & corres,
				 int decompLvl);

	/// destructor
	~Verification(void);

	/// RANSAC method to estimate transformation matrix
	cv::Mat ransac(const std::string & estimateMethod,
				   int num_iterations,
				   double reprojTh,
				   bool normalize);

	/// RANSAC method to estimate transformation matrix (matrices)
	/// between each cluster
	std::vector<cv::Mat> ransac(const std::string & estimateMethod,
								std::vector<std::vector<std::pair<cv::Point,
									cv::Point> > > clusters,
								 int num_iterations,
								 double reprojTh,
								 bool normalize);
	/** \brief apply actual ransac algorithm to estimate transformation matrix
		\return # inliers
	 */
	static int ransac(const std::string & estimateMethod,
					  const cv::Mat &from,
					  const cv::Mat &to,
					  bool normalize,
					  cv::Mat & trafo,
					  int num_iterations,
					  double reprojTh,
					  std::vector<uchar> & inliers);

	/** \brief computes the 3x3 affine homogenous transformation matrix
			as described in the gold-standard algorithm by Hartley p.130 (2nd ed)
		\param[in] from, to matrices of points
		\param[out] trafo resulting transformation matrix
	 */
	static void affinity(const cv::Mat & from,
						 const cv::Mat & to,
						 cv::Mat &trafo);

	/// verifies the map by testing the shift vectors
	/*! if numMax > 0 then the numMax-maximum of the shiftvectors will be searched and all blockpairs
	*	which < minSameShift to this max will be erased.
	*	\param minSameShift the minimum amount of blocks with the same shift vector
	*/
	void verificateShift(size_t minSameShift, int numMax, int shiftVariance);

	/*! marks points if its neighborhood is forged too i
   */
	//	void verificateNeighbors(int minFirstNeighbor, int minSecondNeighbor);

	/// verifies the map by testing the neighbourhoodArea of a block
	/*!
	*	\param number number of neighbouring blocks which will be analyzed
	*	\param maxDistance maximum distance of pixels from the analyzed one
	*	\param areaPercentage threshold of how many blocks have to fullfill the similarity condition in the
	*		near area (defined by distance-param)
	*	\param step the step size between two blocks
	*/
	//	void verificateArea(int number, int maxDistance, double areaPercentage, int step);

	/*!
   * \return a matrix of similar blocks
   */
	cv::Mat & getSimilarityMatrix(void);
	/*!
   * \return the set of similar block pairs
   */
	inline std::vector< std::multimap<std::pair<int,int>,
			std::pair<cv::Point,cv::Point> > > & getShift()
	{
		return shiftvec;
	}

	/*! \brief applies normalization to the clusters such that the data
		gets zero-mean and the avg-distance to the center is sqrt(2).
		\param src[in,out] matrix of of points which will be transformed
		\param t[out] transformation matrix
	*/
	static void normalizeCluster(cv::Mat & src, cv::Mat & t);
	/// renormalize the transformation matrices as described in Hartley
	/*! renormalization has to be done to the output homography matrix
		i.e. we actually need that transformation in matrix form to get the
		actually homography matrix H from the H_{norm} which was obtained
		when using the normalized points:H =  T'{-1} H_{norm} T
		s. Hartley ~p.110
	*/
	static void denormalizeH( const cv::Mat & t1, const cv::Mat & t2, cv::Mat & trafo );

	/// merge similar transformations and recompute the transformation matrix
	static void mergeTrafos(std::vector<cv::Mat> & all_trafos,
							std::vector<std::vector<cv::Vec2f> > & all_inliers,
							/* ransac parameters */
							const std::string & estimateMethod,
							bool normalize,
							int num_iterations,
							double reprojTh);

	/// root mean square error between a and b
	static double rmse( const cv::Mat & a, const cv::Mat & b );
	static std::string printTrafo(const cv::Mat & trafo);
	void printSet() const;
private:
	/// convert correspondences to shift-vectors
	void convertCorresToShift();

	/// can we skip the trafo due to abstruse scale?
	bool skip(const cv::Mat & _trafo) const;

	std::vector<cv::Mat_<cv::Point> > & corres;
	int block_size;
	/// dimension in x direction of the matrix ( == img.width)
	int dimx;
	/// dimension in y direction of the matrix ( == img.height)
	int dimy;
	/// map from shift vector to Blockpair, so we can easily get all Blocks with a special shift vector
	std::vector<std::multimap<std::pair<int, int>, std::pair<cv::Point,cv::Point> > > shiftvec;
	/// matrix (size = dimx*dimy = img_width*img_height)
	/// used for the marking step
	cv::Mat matrix;
	/// decomposition level (wvelet or gaussian pyramid)
	int decomposeLvl;
	static bool cmp(int a, int b){
		return a > b;
	}
};
#endif
