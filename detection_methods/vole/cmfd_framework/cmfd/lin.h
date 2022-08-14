/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef LIN_H
#define LIN_H

#include "feature.h"

class Lin : public Feature {
	public:
		Lin(struct ft_cfg _config,
			const BlockHandling &b,
			int numThreads);
        ~Lin(){}
		void computeOne(const cv::Mat & current_block, int row);
	private:
};
#endif
