/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef LUO_H
#define LUO_H

#include "feature.h"

class Luo : public Feature {
	public:
		Luo(struct ft_cfg _config, const BlockHandling &b,
			int numThreads);
        ~Luo(){}
		void computeOne(const cv::Mat & cur, int row);
	private:
};
#endif
