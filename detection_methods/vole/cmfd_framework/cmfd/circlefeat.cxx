/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "circlefeat.h"
#include "execution.h"
Circle :: Circle(struct ft_cfg _config, const BlockHandling & blocks, int numThreads)
	:	Feature(_config, blocks, numThreads)
{
	allocateFeature(block_size / 2);
}

void Circle :: computeOne(const cv::Mat & current, int row)
{
	// take as many features as blocksize allows
	int num_radius = block_size / 2;
	for (int radius = 1; radius <= num_radius; radius++){
		// create mask
		cv::Mat_<uchar> mask(current.rows, current.cols, static_cast<uchar>(0));
		cv::circle(mask, cv::Point(num_radius, num_radius),
				radius, cv::Scalar(1,0), -1);

		// mean of the disc
		feature(row, radius-1) = cv::mean(current, mask)[0];
	}
}
