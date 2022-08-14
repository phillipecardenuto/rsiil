/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef MARK_H
#define MARK_H

#include "cv.h"
#include "copymoveframework_config.h"
#include <map>
#include <vector>
#include <string>

// Forward declarations
class Matching;
/*! \brief marks somehow the copied block pairs in the image
 */
class Mark
{
public:
	Mark(cv::Mat & _img,
		 cv::Mat & _origimg,
		 struct mark_cfg &_config,
		 const std::vector<cv::Mat_<cv::Point> > & _corres,
		 int _blockSize);
	void markIt(cv::Mat& matrix);
	void markIt(std::vector< std::multimap<std::pair<int,int>,
				std::pair<cv::Point, cv::Point> > > shift);
	/// marks each region w. different color and shift-vecs if whished with green
	void markRegions(std::vector<std::vector<std::pair<cv::Point,cv::Point> > > & clusters,
					 bool draw_shift_vecs);
	void markContour(const cv::Mat& cmfd_map,
					 std::vector<std::vector<cv::Point> > * contours,
					 std::vector<cv::Vec4i> * hierarchy);
	/// marks the convex hull
	void markHull(const std::vector<std::vector<std::pair<cv::Point,cv::Point> > > & clusters);
	/// write pairs - marks pairs in GT differently
	void writePairs(const std::string & groundTruthFile, bool shift);
	/// write chromaticity matrix of positions
	void writeChrom(void);
	/// computes correlation map from transformations
	std::vector<cv::Mat_<uchar> > computeCorrelation(const std::vector<cv::Mat> & trans,
													 bool write);
	/// just writes out the matrix given, nothing fancy here
	static void writeMatrix(const cv::Mat & matrix);
	/// dumps a matrix with the unique identifier
	static void dumpMatrix(std::string identifier, const cv::Mat & dump_matrix);
	// Mark Detected region with an unique ID
    void labelRegions(
            std::vector<std::vector<std::pair<cv::Point, cv::Point> > > &clusters,
            const cv::Mat &post);
private:
std::vector<cv::Point> getNeighbors(const cv::Mat &img, const cv::Point &p) const
{
	int radius = 2;
	std::vector<cv::Point> neighbors;
	for(int x = std::max(p.x-radius,0); x < std::min(p.x+radius+1, img.cols); ++x) 
			for(int y = std::max(p.y-radius,0); y < std::min(p.y+radius+1, img.rows); ++y) {
				if(x==0 & y==0)
					continue;
				if (x>= img.cols || x < 0)
					continue;
				if (y>= img.rows || y < 0)
					continue;
				neighbors.push_back(cv::Point(x,y));
			}

	/*/ add the neighbors
	if ( p.y+1 < img.rows )
		neighbors.push_back(cv::Point(p.x,	 p.y+1));	// mid below
	if ( p.y-1 >= 0)
		neighbors.push_back(cv::Point(p.x,	 p.y-1));	// mid above
	if ( p.x-1 >= 0)
		neighbors.push_back(cv::Point(p.x-1, p.y));		// left mid
	if ( p.x+1 < img.cols)
		neighbors.push_back(cv::Point(p.x+1, p.y));		// right mid

	// for eight-neighborhood:
	if ( p.x-1 >= 0 && p.y+1 < img.rows )
		neighbors.push_back(cv::Point(p.x-1, p.y+1));	// left below
	if ( p.x+1 < img.cols && p.y+1 < img.rows )
		neighbors.push_back(cv::Point(p.x+1, p.y+1));	// right below
	if ( p.x-1 >= 0 && p.y-1 >= 0 )
		neighbors.push_back(cv::Point(p.x-1, p.y-1));	// left above
    if ( p.x+1 < img.cols && p.y-1 >= 0 )
        neighbors.push_back(cv::Point(p.x+1, p.y-1));	// right above
		*/

	return neighbors;
};

	cv::Mat_<float> computeCorrelMap(const cv::Mat_<uchar> &plane,
									 const cv::Mat &trans,
									 bool inverse = false);
	int block_size;
	cv::Mat & img;
	cv::Mat & origimg;
	struct mark_cfg config;
	const std::vector<cv::Mat_<cv::Point> > & corres;
	/// png-parameters
	std::vector<int> params;
};

#endif // MARK_H
