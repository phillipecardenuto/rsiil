/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef PCA_H
#define PCA_H

#include "cv.h"
#include "feature.h"

class PCA {
	public:
		/// eps will be used to set new dimension (s. reduce) if no dimension is set
		PCA(Feature & feat, double eps);
		~PCA(){}

		/*! performs pca and reduces the dimension of the feature vectors
		 * \param dim new dimension of featurespace (= featuresize)
		 *  if dim isnt set or set to 0, the new dimension will be set, so it satisfies
		 *  1 - eps = sum[0,dim](evals) / sum[0,ftSize](evals)
		 */
		void reduce(int dim, int minFeatureSize, int maxFeatureSize);
		/*! check if evecs contains nan, if so, we remove the whole row
		 * from the evecs-matrix and returns new dimension
		 */
		static int checkNan(cv::Mat & evecs);
		/// returns dimension with given error, between [minFeatureSize,maxFEatureSize]
		static int getDimension(const cv::Mat & evals,
								double error,
								int minFeatureSize,
								int maxFeatureSize);
	private:
		Feature & feat;
		/// reference to the blocks
		cv::Mat_<float> & features;
		/// the accuracy of diagonalization
		double eps;
};


#endif
