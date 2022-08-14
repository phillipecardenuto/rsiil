/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef PCTFEAT_H
#define PCTFEAT_H

#include "feature.h"

class Pctfeat : public Feature {
	public:
		Pctfeat(struct ft_cfg _config, const BlockHandling &b, int numThreads)
			: Feature(_config, b, numThreads){}
        ~Pctfeat(){}
		void computeFeatureVec(void);
	private:
		void computeMean(int begin, int finish, cv::Mat_<double> & mean);
		cv::Mat_<float> computeCov(cv::Mat_<double> & mean);
		int dim;
};
#endif
