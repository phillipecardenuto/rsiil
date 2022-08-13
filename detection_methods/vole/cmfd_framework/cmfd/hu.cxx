/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "hu.h"
#include <vector>
#include <cmath>
#include "cv.h"
#include "blur.h"

Hu :: Hu(struct ft_cfg _config, const BlockHandling &b,int threadNum)
	:	Feature(_config, b, threadNum)
{
	CV_Assert(config.numHu > 0 && config.numHu < 8);
	allocateFeature(config.numHu);
}

void Hu :: computeOne(const cv::Mat & cur, int row)
{
	cv::Mat block = cur;

	// mask
	if (config.circle){
		// make deep copy
		block = cur.clone();

		cv::Mat_<uchar> mask(cur.rows, cur.cols, static_cast<uchar>(0));
		cv::circle(mask, cv::Point(cur.cols/2, cur.rows/2),
				   static_cast<uchar>(cur.rows/2), cv::Scalar(1,0), -1);
		for (int y = 0; y < cur.rows; y++){
			for (int x = 0; x < cur.cols; x++){
				if (!mask(y,x))
					block.at<uchar>(y,x) = 0.0;
			}
		}
	}
	// compute moments to degree 3
	cv::Moments mom = cv::moments(block);
	// compute 7 Hu-moments
	double hu[7];
	cv::HuMoments(mom, hu);

	// use as many as adjusted
	for (int n = 0; n < config.numHu; n++){
		feature(row, n) = hu[n];
	}
}
