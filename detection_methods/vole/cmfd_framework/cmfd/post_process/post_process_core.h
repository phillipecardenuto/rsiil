#ifndef POST_PROCESS_CORE_H
#define POST_PROCESS_CORE_H

#include "post_process_config.h"
#include <opencv2/core/core.hpp>

/*! \brief post processing of cmfd output
 * This class is the core-class for all post-processing applications
 */
class PostProcess
{
public:
	PostProcess(const PPConfig & config, const cv::Mat & cmfd_img);
	~PostProcess();

	/*! \brief apply an area threshold of a cmfd output image
	 * If a detected region of a copy-move forgery detection (cmfd)
	 * is larger than an area threshold (config.areaThreshold)
	 * than this area will be kept, otherwise it will be deleted
	 *	\param matrix a possible matrix of feature matches of a cmfd
	 *		output, if this matrix is not empty it will be used to
	 *		verify the path, if no match occurs on the path it is rejected
	 */
	void applyAreaThreshold(const cv::Mat & matrix,
							std::vector<cv::Mat_<cv::Point> > * corres = NULL);

	/*! \brief applies morphological filter to the cmfd img
	 * possible options: DILATE, ERODE, OPEN, CLOSE, etc (see opencv-doku)
	 */
	void applyMorphFilter();
	/*!	fill holes in the cmfd-map by extracting the contours, taking
	 *	only the outer contours and filling them
	 *  \param[in,out] contours extracted contours
	 *	\param[in,out] hierarchy hierarchy of contours
	 */
	void fillHolesByContour(std::vector<std::vector<cv::Point> > & contours,
							std::vector<cv::Vec4i> & hierarchy);

	/// returns the (possibly) manipulated cmfd-map
	inline const cv::Mat & getMatrix() const{
		return cmfd_img;
	}

private:
	/*! \brief helper class for applyAreaThreshold
	 * gives neighbors of a point
	 * \param img plane of cmfd_img
	 * \param p point around which neighbors are gained
	 * \return neighbors of the point 
	 */
	std::vector<cv::Point> getNeighbors(const cv::Mat &img, const cv::Point &p) const;
	/// unmark a path in the cmfd_img
	void unmark(const std::vector<cv::Point> & path, int c);
	/// the post processing config
	const PPConfig & config;
	/// the cmfd output map
	cv::Mat cmfd_img;
};

#endif
