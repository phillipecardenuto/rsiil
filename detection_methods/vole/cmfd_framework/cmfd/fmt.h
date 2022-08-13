/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef FMT_H
#define FMT_H

#include "feature.h"
/*!
 * uses fourier mellin transformation as features for blocks
 */
class Fmt : public Feature 
{
	public:
		Fmt(struct ft_cfg _config, const BlockHandling &b, int numThreads);
        ~Fmt(){}
		void computeOne(const cv::Mat & current_block, int row);
	private:
		double optiM;
		int dstwidth;
		int feature_size;
};
#endif
