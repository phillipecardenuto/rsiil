/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef SVD_H
#define SVD_H

#include "feature.h"

class Svd : public Feature {
	public:
		Svd(struct ft_cfg _config, const BlockHandling &b, int numThreads);
       ~Svd(){}
		void computeOne(const cv::Mat & cur, int row);
	private:
		int computeDim(cv::Mat_<double> &chan);
		void project(int row, cv::SVD &s, int dim);
		size_t globalDim;
		/// maximal dimension if dim == -2
		int max_dim;
};
#endif
