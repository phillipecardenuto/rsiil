/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef DWTFEAT_H
#define DWTFEAT_H

#include "feature.h"

class Dwtfeat : public Feature 
{
	public:
		Dwtfeat(struct ft_cfg config,
				const BlockHandling & blocks,
				int numThreads);
        ~Dwtfeat(){}
		void computeOne(const cv::Mat & cur, int row);
	private:

};
#endif
