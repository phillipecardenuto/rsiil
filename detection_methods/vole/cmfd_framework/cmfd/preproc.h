/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef PREPROC_H
#define PREPROC_H

#include "cv.h"
#include "copymoveframework.h"
class Preproc {
	public:
		/// Constructor
		// uses _config to process image
		Preproc( cv::Mat &i , struct pr_cfg config);
		/// computes gaussian pyramid and gives back the img with the 
		/// level specified in the config
		cv::Mat computeGaussianPyramid(void);
		/// computes haar wavelet to specified level and gives back the img
		/// with the level specified in the config
		cv::Mat computeHaarWavelet(void);
    private:
		cv::Mat &img;
		struct pr_cfg config;
};

#endif
