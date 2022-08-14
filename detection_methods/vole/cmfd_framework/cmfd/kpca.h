/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef KPCA_H
#define KPCA_H

#include "feature.h"

class Kpca : public Feature {
	public:
		/// eps will be used to set new dimension (s. reduce) if no dimension is set
		Kpca(struct ft_cfg _config, const BlockHandling &b, int numThreads);
		~Kpca(void){}

		/*! performs Kpca and reduces the dimension of the feature vectors
		 * \param dim new dimension of featurespace (= featuresize)
		 *  if dim isnt set or set to 0, the new dimension will be set, so it satisfies
		 *  1 - eps = sum[0,dim](evals) / sum[0,ftSize](evals)
		 */
		void computeFeatureVec(void);
	private:
		cv::Mat_<float> selectTrainingSamples(void);
		cv::Mat_<float> getTrainingKernel(const cv::Mat_<float> & data);
		cv::Mat_<float> getTestKernel(const cv::Mat_<float> & A_trn);
		/// the kernel function (Gaussian Kernel)
		float kernel(const cv::Mat & i, const cv::Mat & j);
		cv::Mat_<float> centerKernel(const cv::Mat_<float> & K);
		// center testing Kernel
		cv::Mat_<float> centerTestKernel(const cv::Mat_<float> & K_tst,
										 const cv::Mat_<float> &K_trn);
		int getDimension(const cv::Mat & evecs);
		/// the accuracy of diagonalization
		double eps;
		/// dimension of features
		int dim;
		/// variance of Gaussian Kernel
		double sigma;
		/// number of training samples
		int num_samples;
};

#endif
