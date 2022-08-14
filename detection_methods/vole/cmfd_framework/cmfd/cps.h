/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef CPS_H
#define CPS_H

#include "cv.h"
#include <vector>
#include "copymoveframework.h"
/*! 
 * tries to find copied regions with help of cross power spectrum
 */
class Cps
{
	public:
		/// constructor, subImg = divider in direction -> subimg_x*subimgy = num of parts
		Cps(cv::Mat &img, cv::Mat &origimg, struct ft_cfg cfg, int waveletLvl);
		~Cps(void){}
		/// divide img into non overlipping sub images
		void computeSubImg(void);
		/// computes the (optionally not) normalized cross power spectrum
		static void crossPowerSpectrum(cv::Mat& src1, cv::Mat& src2, cv::Mat& dest, bool normalize = true);
		/// finds duplicated region with help of normalized CPS
		void findRegion(void);
		/// returns the similarity matrix
		cv::Mat_<cv::Vec3b>& getSimilarityMatrix(void) {
			return matrix;
		}
	private:
		cv::Mat &img;
		cv::Mat &origimg;
		/// wavelet multiplier
		int waveletMult;
		/// parts the image will be divided in
		int divx;
		int divy;
		double correlTh;
		int pixelDiff;	
		/// vector of subimages of the original image
		std::vector<cv::Mat > subimages;
		cv::Mat_<cv::Vec3b> matrix;
};
#endif
