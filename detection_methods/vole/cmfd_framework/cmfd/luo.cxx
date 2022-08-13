/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include <iostream>
#include <cmath>

#include "luo.h"
#include "execution.h"

Luo :: Luo(struct ft_cfg _config, const BlockHandling &b,
	int numThreads)
	:	Feature(_config, b, numThreads)
{
	allocateFeature(7);
}

void Luo :: computeOne(const cv::Mat & cur, int row)
{
	if ( cur.channels() != 3 ) {
		throw std::runtime_error("Luo::computeOne: You need to use all three channels, i.e. don't use grayscale-conversion");
	}

	// Averages of every channel (typically rgb)
	cv::Scalar m = cv::mean(cur);
	for (int nr = 0; nr < 3; nr++){
		feature(row, nr) = m[nr];
	}


	double sum = cv::sum(cur)[0]; // compute sum
	// prevent division by zero if sum == 0
	if(sum == 0.0){
		for (int nr = 3; nr < feature.cols; nr++){
			feature(row, nr) = 0.0;
		}
		return;
	}

	// the sump = sum(part1 / sum(part1 + part2))
	// 	= sum(part1 / sum)
	// TODO: probably these 8 for-loops may be merged together
	double sump = 0.0;
	for (int y = 0; y < static_cast<int>(cur.rows / 2); y++){
		for (int x = 0; x < cur.cols; x++){
			sump += cur.at<uchar>(y,x);
		}
	}
	feature(row, 3) = sump / sum;

	sump = 0.0;
	for(int y = 0; y < cur.rows; y++){
		for(int x = 0; x < static_cast<int>(cur.cols / 2); x++){
			sump += cur.at<uchar>(y,x);
		}
	}
	feature(row, 4) = sump / sum;

	sump = 0.0;
	for(int y = 0; y < cur.rows; y++){
		for(int x = y; x < cur.cols; x++){
			sump += cur.at<uchar>(y,x) ;
		}
	}
	feature(row, 5) = sump / sum;

	sump = 0.0;
	for(int y = 0; y < cur.rows; y++){
		for(int x = 0; x < (cur.cols - y); x++){
			sump += cur.at<uchar>(y,x) ;
		}
	}
	feature(row, 6) = sump / sum;
}
