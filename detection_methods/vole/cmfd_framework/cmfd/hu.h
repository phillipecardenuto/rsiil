/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef HU_H
#define HU_H

#include "feature.h"

class Hu : public Feature {
	public:
		Hu(struct ft_cfg _config, const BlockHandling &b,
		   int threadNum);
		~Hu(){}
		void computeOne(const cv::Mat & cur, int row);
	private:
};
#endif
