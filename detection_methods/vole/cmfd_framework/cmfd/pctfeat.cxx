/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "pctfeat.h"
#include "execution.h"
#include "pca.h"
#include <cassert>

#include <boost/thread.hpp>
#include <boost/bind.hpp>

void Pctfeat::computeMean(int begin, int finish, cv::Mat_<double> & mean)
{
	for (int i = begin; i < finish; i++){

		cv::Mat block = blocks[i];

		cv::Mat_<double> tmp;
		if (config.circle){	
			cv::Mat_<double> tmpsrc = block;
			cv::Mat_<double> mask(block.rows, block.cols, 0.0);
			cv::circle(mask, cv::Point(block.cols/2, block.rows/2),
					   static_cast<int>(block.rows/2), cv::Scalar(1.0,0.0), -1);
			double newSize = cv::sum(mask)[0];
			tmp = cv::Mat_<double>(1,newSize);
			int cnt = 0;
			for (int y = 0; y < block.rows; y++){
				for (int x = 0; x < block.cols; x++){
					if (mask(y,x) != 0)
						tmp(0,cnt++) = static_cast<double>(tmpsrc(y,x));
				}
			}
		}
		else {
			tmp = block;
			tmp = tmp.reshape(1,1);
		}
		
		mean += tmp;
	}
}

// compute Covariance matrix
// note: this is indeed a strange version of computing the Covar
// more sense: cov(matrix of all blocks)
cv::Mat_<float> Pctfeat::computeCov(cv::Mat_<double> & mean)
{
	cv::Mat_<double> tmpdst;
	cv::Mat_<double> cov(mean.cols, mean.cols, 0.0);

	for (size_t i = 0; i < blocks.size();  i++){
		cv::Mat_<double> block = blocks[i];
		block = block.reshape(1,1);

		// computes (x-mean)^T .* (x-mean)
		// not (x-mean) .* (x-mean)^T as the features are row-wise
		// in the matrix, not column-wise
		cv::mulTransposed(block, tmpdst, true, mean);

		cov += tmpdst;
	}
	cv::Mat_<float> cov_flt = cov / blocks.size();
	return cov_flt;
}

void Pctfeat::computeFeatureVec(void)
{
	// preperation for threading
	std::vector<cv::Mat_<double> > means;
	means.resize(num_threads);

	// the real mean and covar matrix
	cv::Mat_<double> mean = cv::Mat_<double>::zeros(1, blocks[0].rows *  // 0 matrix
													blocks[0].cols);// *

	// to compute pca only from a circular block
	if (config.circle){
		cv::Mat block = blocks[0];
		cv::Mat_<double> mask(block.rows, block.cols, 0.0);
		cv::circle(mask, cv::Point(block.cols/2, block.rows/2),
				   static_cast<int>(block.rows/2), cv::Scalar(1.0,0.0), -1);
		int newSize = static_cast<int>(cv::sum(mask)[0]);
		mean = cv::Mat_<double>::zeros(1, newSize);
	}

	for (int i = 0; i < num_threads; i++) {
		means[i] = mean.clone();
	}

	// boundings of threads
	unsigned int vectorsize = blocks.size();
	int len = vectorsize / num_threads;
	int rest = vectorsize % num_threads;
	int begin = len + rest;

	if (Execution::verbosity >= 2){
		std::cout << " compute mean\n";
	}

	// compute mean with threads
	{
		boost::thread_group thrd_grp;
		// first thread computes the overhead addiotionally
		boost::thread *t = new boost::thread(boost::bind(
												 &Pctfeat::computeMean, this, 0, begin, boost::ref(means[0]) ));
		thrd_grp.add_thread(t);
		for (int i = 1; i < num_threads; i++){
			t = new boost::thread(boost::bind(
									  &Pctfeat::computeMean, this, begin + (i-1)*len, begin + i*len, boost::ref(means[i])));
			thrd_grp.add_thread(t);
		}
		thrd_grp.join_all();
	}

	// compute overall mean
	for (int i = 0; i < num_threads; i++) {
		mean += means[i];
	}
	mean /= blocks.size();

	if (Execution::verbosity >= 2){
		std::cout << " compute cov\n";
	}

	cv::Mat_<float> cov = computeCov(mean);

	// compute eigenvectors and eigenvalues
	if (Execution::verbosity >=2){
		std::cout << " compute evals and evecs\n";
	}
	cv::Mat evecs;
	cv::Mat evals;
	cv::eigen(cov, evals, evecs);

	if (Execution::verbosity >=2){
		for (int y = 0; y < evals.rows; y++) {
			for (int x = 0; x < evals.cols; x++) {
				std::cout << evals.at<double>(y,x) << " ";
			}
			std::cout << std::endl;
		}
	}
	// get new dimension, if it's set to 0
	dim = config.dim;
	if (dim == 0){
		dim = PCA::getDimension(evals, config.eps, config.minFeatureSize, config.maxFeatureSize);
	}
	if (dim < evals.rows){
		// change eigenvecs & eigenvalues if dim changes
		evals = evals.rowRange(0,dim).clone();
		evecs = evecs.rowRange(0,dim).clone();
	}

	// check if we produced nan, if so, we remove the whole row
	// from the evecs-matrix
	dim = PCA::checkNan(evecs);

	allocateFeature(dim);

	// project blocks to the new space
	if (Execution::verbosity >= 1){
		std::cout << " project to the new space with dimension: " << dim << std::endl;
	}

	for (unsigned int i = 0; i < vectorsize; i++){
		// get correct channel in double
		cv::Mat_<float> block = blocks[i];
		block = block.reshape(1,1);
		cv::Mat_<float> mean_flt = mean;
		block -= mean_flt; // make data zero-mean

		for (int d = 0; d < dim; d++){
			feature(i,d) = block.dot(evecs.row(d));
		}
	}
}
