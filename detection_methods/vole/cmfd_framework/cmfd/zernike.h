/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef ZERNIKE_H
#define ZERNIKE_H

#include "feature.h"

class GetZernike2;

class Zernike : public Feature {
	public:
		Zernike(struct ft_cfg _config, const BlockHandling &b, int numThreads);
		static unsigned long long fak(int n);
		double zernikePoly(int n, int m, double rho);
		void computeOne(const cv::Mat & cur, int row);
	private:
		std::vector<cv::Mat_<cv::Vec2d> > zernikeMoments;
		double pi;
		/// size of one feature-vector
		int feature_size;
};
#endif
