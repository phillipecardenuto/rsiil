/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef CIRCLE_H
#define CIRCLE_H

#include "feature.h"

class Circle : public Feature {
	public:
		Circle(struct ft_cfg _config, const BlockHandling &b, int numThreads);
        ~Circle(){}
		void computeOne(const cv::Mat & current_block, int row);
	private:
};
#endif
