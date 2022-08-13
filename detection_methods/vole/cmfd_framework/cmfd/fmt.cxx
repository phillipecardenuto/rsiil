/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "fmt.h"

#include <cmath>
#include "highgui.h"

Fmt :: Fmt(struct ft_cfg _config, const BlockHandling &b, int numThreads)
	:	Feature(_config, b, numThreads)
{
	feature_size = 45;
	dstwidth = config.dstwidth;
	optiM = Feature::optiM(dstwidth);

	allocateFeature(feature_size);
}
void Fmt :: computeOne(const cv::Mat & current, int row)
{
	cv::Mat_<double> A = current; // convert to double

	//--- Discrete fourier transform to ensure translation invariance
	cv::Mat dest;
	cv::dft(A, dest, cv::DFT_COMPLEX_OUTPUT);

	// compute the complex magnitude
	cv::Mat thedest(dest, cv::Range::all(), cv::Range(0, (dest.cols/2) + 1));
	std::vector<cv::Mat> realcomp;
	cv::split(thedest, realcomp);
	cv::Mat real = realcomp[0];
	cv::Mat comp = realcomp[1];

	cv::Mat mag;
	cv::magnitude(real, comp, mag);

	//--- Log polar transformation
	cv::Mat_<float> mag2 = mag;
	CvMat val(mag2);
	//cv::Mat_<double> logpolar(mag.rows, mag.cols);
	cv::Mat_<float> logpolar(180, dstwidth);
	//			cv::Mat_<float> logpolar(90, dstwidth);
	CvMat logp(logpolar);

	// chose the center -> could perhaps give
	// better results with rotated regions?
	if (config.modified)
		cvLogPolar( &val, &logp, cvPoint2D32f(mag.cols/2,mag.rows/2), optiM, CV_WARP_FILL_OUTLIERS);
	else
		cvLogPolar( &val, &logp, cvPoint2D32f(0,0), optiM, CV_WARP_FILL_OUTLIERS);

//		cvLogPolar( &val, &logp, cvPoint2D32f(0,0), optiM, CV_WARP_FILL_OUTLIERS + CV_INTER_LINEAR);
	
	// iterate over every second angle and additionally take angle + 45\degree
	for (int col = 0; col < feature_size; col++) {
		feature(row,col) = 0.0;
		for (int x = 0; x < logpolar.cols; x++){
			double l1 = fabs(logpolar(col+2,x));
//			double l1 = fabs(logpolar(y,x));
			if (l1 == 0.0)
				l1 = 1;
			double l2 = fabs(logpolar(col+2 + 2*feature_size,x));
//			double l2 = fabs(logpolar(y + feature_size,x));
			if (l2 == 0.0)
				l2 = 1;

			feature(row,col) += log(l1) + log(l2);

//			rows(0,y) += log(fabs(logpolar(y,x)) + 1)
//					+ log(fabs(logpolar(y+feature_size,x)) + 1)
//					+ log(fabs(logpolar(y+2*feature_size,x)) + 1)
//					+ log(fabs(logpolar(y+3*feature_size,x)) + 1);

		}
	}

}
