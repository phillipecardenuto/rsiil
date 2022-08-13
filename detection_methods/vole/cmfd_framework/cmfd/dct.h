/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef DCT_H
#define DCT_H

#include "feature.h"
//#define BLK_SIZE 16
class Dct : public Feature {
	public:
		Dct(struct ft_cfg _config, const BlockHandling &b,int threadNum);
        ~Dct(){}
		void computeOne(const cv::Mat & cur, int row);
	private:
//		static const double quantizationMatrix[8][8];
//		double quant16[BLK_SIZE][BLK_SIZE];
		cv::Mat_<double> quant8;
		cv::Mat_<double> quant16;

		void buildQuantizationMatrix(void);
//		double cos_terms[BLK_SIZE][BLK_SIZE];
		cv::Mat_<double> cos_terms;
		double q;
		const double PI;
};
#endif
