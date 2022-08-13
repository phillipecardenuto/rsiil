/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "pca.h"
#include "execution.h"

PCA::PCA(Feature & _feat, double eps)
	: feat(_feat), features(feat.getFeature())
{
	this->eps = eps;
}

int PCA :: checkNan(cv::Mat & evecs)
{
	while(1){
		cv::Point pt;
		if ( cv::checkRange(evecs, true, &pt, -FLT_MAX, FLT_MAX) ){
			break;
		}
		if (Execution::verbosity >= 2)
			std::cout << " found nan at " << pt.x << " , " << pt.y
					  << " evec(pt): " << evecs.at<float>(pt.y, pt.x)
					  << " -> remove row" << std::endl;
//		std::cerr << evecs;
		// as the evals are sorted we say that after we got 'nan' all other
		// rows won't give additional information thus remove them
		evecs = evecs.rowRange(0, pt.y).clone();
		if (evecs.empty() || evecs.rows == 0){

			throw std::runtime_error("KPCA: ComputeFeatureVec(): evecs are empty after removing nan");
		}
		if (Execution::verbosity >= 1)
			std::cout << " removed nan values, new evecs.rows: " << evecs.rows << std::endl;
	}
	return evecs.rows;
}

int PCA :: getDimension(const cv::Mat & evals, double error,
						int minFeatureSize, int maxFeatureSize)
{
	double sum_evals = cv::sum(evals)[0];
	int err = 1 - error;
	double sum_trunc = 0.0;
	int num_dim = 0;
	for (int i = 0; i < evals.rows; i++){
		sum_trunc += evals.at<float>(i,0);
		num_dim++;
		if ( (sum_trunc / sum_evals) >= err ){
			break;
		}
	}
	if (num_dim > maxFeatureSize){
		num_dim = maxFeatureSize;
	}
	if (num_dim < minFeatureSize){
		num_dim = minFeatureSize;
	}
	if (Execution::verbosity >= 2){
		std::cout << "\tcomputed new dimension: " << num_dim << std::endl;
	}
	return num_dim;
}

void PCA :: reduce(int dim, int minFeatureSize, int maxFeatureSize)
{
	cv::PCA p(features, cv::Mat(), CV_PCA_DATA_AS_ROW, dim);
	cv::Mat evals = p.eigenvalues;

	// get new dimension, if it's set to 0
	if (dim <= 0){
		getDimension(evals, eps, minFeatureSize, maxFeatureSize);

		// change eigenvecs & eigenvalues if dim changes
		p.eigenvalues = p.eigenvalues.rowRange(0,dim).clone();
		p.eigenvectors = p.eigenvectors.rowRange(0,dim).clone();
	}
	dim = checkNan(p.eigenvectors);

	if (Execution::verbosity >=2){
		std::cerr << " project to the new space with dimension: " << dim << std::endl;
	}

	cv::Mat_<float> new_feat = p.project(features);

	// free old feature-matrix
	feat.free();
	// assign the new one
	features = new_feat;
}
