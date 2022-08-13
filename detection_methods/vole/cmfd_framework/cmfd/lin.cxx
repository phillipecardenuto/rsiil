/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include <iostream>
#include <cmath>

#include "lin.h"

Lin :: Lin(struct ft_cfg _config,
	const BlockHandling & blocks,
	int numThreads)
	:	Feature(_config, blocks, numThreads)
{
	allocateFeature(9);
}

void Lin :: computeOne(const cv::Mat & current, int row)
{
	// compute f{1} = mean of complete block
	double ave = cv::mean(current)[0];

	// compute mean of sub blocks
	int w = current.cols;
	int h = current.rows;
	double avg[4];
	avg[0] = cv::mean(current(cv::Range(0, h/2), cv::Range(0, w/2)))[0];
	avg[1] = cv::mean(current(cv::Range(0, h/2), cv::Range(w/2, w)))[0];
	avg[2] = cv::mean(current(cv::Range(h/2, h), cv::Range(0, w)))[0];
	avg[3] = cv::mean(current(cv::Range(h/2, h), cv::Range(w/2, w)))[0];

	//get min and max
	cv::Scalar subavg; // pre state of feature 6-9
	double min = 256.0;
	double max = -256.0;
	for (int i = 0; i < 4; i++){
		// subavg[i] = f{5-9} of the paper
		subavg[i] = avg[i] - ave;
		if (subavg[i] < min){
			min = subavg[i];
		}
		if (subavg[i] > max){
			max = subavg[i];
		}
	}

	// Average of channel
	feature(row, 0) = floor(ave);
	for (int i = 1; i < 9; i++) {
		// features 2-5
		if ( i < 5){
			feature(row, i) = floor(255*(avg[i-1] / (ave == 0 ? 1 : ave)));
		}
		// features 6-9
		else {
			double maxMinusMin = max - min;
			double m2 = min;
			if (maxMinusMin == 0.0){
				maxMinusMin = 1.0;
				m2 = 0.0;
			}
			feature(row, i) = floor( 255 * ((subavg[i-5] - m2)
											/ (maxMinusMin == 0 ? 1 : maxMinusMin) ));
		}
	}
}
