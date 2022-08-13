/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef BLUR_H
#define BLUR_H

#include "feature.h"

/// this class uses blur moments as features proposed by Mahdian et al
class Blur : public Feature {
	public:
        Blur(struct ft_cfg config, const BlockHandling & blocks, int num_threads)
			:	Feature(config, blocks, num_threads)
		{}
        Blur(struct ft_cfg config, const BlockHandling & blocks, int num_threads, cv::Mat_<float> &feature)
			:	Feature(config, blocks, num_threads)
		{
			this->feature = feature;
		}
		~Blur(){}
		Blur& clone(void){
			Blur *n = new Blur(config, blocks, num_threads, feature);
			return *n;
		}
		void computeFeatureVec(void);
		void computeRange(int, int);
		void computeOne(const cv::Mat &act, int row);
	
		/// initializes tables with -1
		void initializeTables(void);

		/// computes the centroids for centralized moments
		/// also computes mue00, mue01, mue10
		void computeCentroids(cv::Mat_<double> &chan);

		/// computes centralized moments for a matrix
		double mue(unsigned int p, unsigned int q, cv::Mat_<double> &act);
	private:
		int binom(unsigned int n, unsigned int k) const;
		double blur(unsigned int p, unsigned int q, cv::Mat_<double> &act);
		double mue(unsigned int p, unsigned int q){
			return mue(p,q, chan);
		}
		cv::Mat_<double> chan;
		/// centroid in x dir
		double xt;
		/// centroid in y dir
		double yt;
		// speed up tables
		double blurTable[7][7];
		/// table for faster look up for mues
		double mueTable[7][7];
};
#endif
