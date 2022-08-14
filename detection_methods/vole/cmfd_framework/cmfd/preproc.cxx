/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "preproc.h"
#include "wavelet.h"
#include "execution.h"

Preproc :: Preproc(cv::Mat &i, struct pr_cfg _config)
	:	img(i)
{
	cv::Mat_<uchar> tmp;
	config = _config;
	if (config.chan == "GRAY"){
		cv::cvtColor(img, tmp, CV_RGB2GRAY);	
		img = tmp.clone();
	}
	else if (img.channels() >= 3){
		std::vector<cv::Mat> planes;
		split(img, planes);
		if (config.chan == "RED"){
			img = planes[2].clone();
		}
		else if (config.chan == "GREEN"){
			img = planes[1].clone();
		}
		else if (config.chan == "BLUE"){
			img = planes[0].clone();
		}
		else if (config.chan == "ALL"){
			return;
		}
		else {
			throw std::invalid_argument(std::string(config.chan + " invalid, allowed: ALL, GRAY, GREEN, BLUE, RED"));
		}
	}
}

cv::Mat Preproc :: computeGaussianPyramid(void)
{
	if (config.pyramidLvl <= 0)
		return img;
	if(Execution::verbosity >= 1 )
		std::cerr << "Gaussian Pyramid w. Level " << config.pyramidLvl << std::endl;
	std::vector<cv::Mat> dst;
	cv::buildPyramid(img, dst, config.pyramidLvl);
	img = dst[config.pyramidLvl];
	return img;
}

cv::Mat Preproc :: computeHaarWavelet(void)
{
	if (config.waveletLvl <= 0)
		return img;
	if(Execution::verbosity >=1 ){
		std::cerr << "DWT w. Level " << config.waveletLvl << std::endl;
	}
	std::vector<cv::Mat_<double> > dst;
	//Wavelet::decomposeHaar(img, dst, config.waveletLvl, false); // false: dont normalize result to 0-255
	Wavelet::decomposeHaar(img, dst, config.waveletLvl, true);

	// convert to uchar
	cv::Mat_<uchar> tmp = dst.back();
	img = tmp;

	return img;
}
