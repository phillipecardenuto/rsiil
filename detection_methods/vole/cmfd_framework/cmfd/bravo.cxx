/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include <iostream>
#include <cmath>

#include "bravo.h"
#include "execution.h"

Bravo :: Bravo(struct ft_cfg _config, const BlockHandling &b,
			   int numThreads)
	:	Feature(_config, b, numThreads)
{
	allocateFeature(4);
}

void Bravo :: computeOne(const cv::Mat & cur, int row)
{
	if ( cur.channels() != 3 ) {
		throw std::runtime_error("Bravo::computeOne: You need to use all three channels, i.e. don't use grayscale-conversion");
	}

	// create mask
	cv::Mat_<uchar> mask(cur.rows, cur.cols, static_cast<uchar>(0));
	cv::circle(mask, cv::Point(cur.cols/2, cur.rows/2),
			static_cast<int>(cur.rows/2), cv::Scalar(1,0), -1);

	// mean of the disc
	cv::Scalar m = cv::mean(cur, mask);

	// compute luminance Y
	cv::Mat xyz;	
	cv::cvtColor(cur, xyz, CV_BGR2XYZ);

	// get probabilities of lumi
	int	chan[] = {1};// take only the Y - channel
	int histSize[] = {256}; // take 256 bins -> should we change that?
	float lumirange[] = {0, 256};
	const float* ranges[] = {lumirange};
	cv::MatND hist; // output histogram
	cv::calcHist(&xyz, 1, chan, mask, // use mask
			hist, 1, histSize, ranges, 
			true, // histogram is uniform
			false); // don't akkumulate

	// compute entropy of the disk
	double entropy = 0.0;
	for (int y = 0; y < 256; y++){
		// compute luminance probability 
		float p = hist.at<float>(y) / (cv::sum(mask)[0]);	
		if (p == 0.0) 
			continue;
		// compute entropy with log to basis log_2
		entropy -= p * (log(p) / log(2));
	}

	if ( Execution::verbosity >=3  ){
		std::cerr << "entropy of block(" << row << ") = " << entropy << std::endl;
	}

	feature(row, 0) = m[0];
	feature(row, 1) = m[1];
	feature(row, 2) = m[2];
	feature(row, 3) = entropy;
}
