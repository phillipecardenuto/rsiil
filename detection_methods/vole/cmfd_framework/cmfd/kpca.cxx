/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "kpca.h"
#include "execution.h"
#include "pca.h"

#include <stdlib.h>
#include <time.h>

Kpca :: Kpca(struct ft_cfg config, const BlockHandling &b, int numThreads)
	: Feature(config, b, numThreads)
{
	eps = config.eps;
	sigma = config.sigma;
	num_samples = config.numSamples;
	dim = config.dim;
	cv::Mat bl = blocks[0];
}

cv::Mat_<float> Kpca :: selectTrainingSamples(void)
{
	cv::Mat_<float> A_trn(num_samples, block_size*block_size);
	// we use training samples at equidistant points
	double distance = blocks.size() / num_samples;

	for (int i = 0; i < num_samples; i++){
		int r = static_cast<int>( i*distance );
		// convert to float

		cv::Mat_<float> rand_block = blocks[r];
		// matrix -> vector
		rand_block = rand_block.reshape(1,1);
		// copy it to our training matrix
		cv::Mat_<float> A_row = A_trn.row(i);
		rand_block.copyTo(A_row);
	}
	int non_zero = cv::countNonZero(A_trn);
	if (non_zero == 0)
		std::cerr << "WARNING: Kpca :: selectTrainingSamples: no non_zero values\n";
	return A_trn;
}

// Gaussian Kernel
// K(xi, xj) = exp(-||(xi-xj)||^2 / 2*sigma^2)
float Kpca :: kernel(const cv::Mat & i, const cv::Mat & j)
{
	double norm = cv::norm(i, j);
	return static_cast<float>( std::exp(-norm / (2*sigma*sigma)) );
}

cv::Mat_<float> Kpca :: getTrainingKernel(const cv::Mat_<float> & data)
{
	int rows = data.rows;
	// create kernel matrix
	cv::Mat_<float> K(rows, rows);
	for (int i = 0; i < rows; i++) {
		for (int j = i; j < rows; j++) {
			float k = kernel(data.row(i), data.row(j));
			K(i, j) = k;
			// Kernel matrix is symmetric
			K(j, i) = k;
		}
	}
	return K;
}


cv::Mat_<float> Kpca :: centerKernel(const cv::Mat_<float> & K)
{
	cv::Mat_<float> U(K.rows, K.cols, 1.0/K.rows);

	return K - U*K - K*U + U*K*U;
}

// K_ij = kernel(test(i),train(j))
cv::Mat_<float> Kpca :: getTestKernel(const cv::Mat_<float> & A_trn)
{
	vout(1) << "Create Test Kernel matrix of dimension: "
				  << blocks.size() << " x " << num_samples << std::endl;

	// create Kernel-matrix of size (N x M), where N = all_blocks,
	// and M = number of training samples
	cv::Mat_<float> K(blocks.size(), num_samples);

	for (int y = 0; y < K.rows; y++) {
		for (int x = 0; x < K.cols; x++) {
			cv::Mat_<float> tmp = blocks[y];
			tmp = tmp.reshape(1,1);
			K(y, x) = kernel( tmp, A_trn.row(x) );
		}
	}
	return K;
}

cv::Mat_<float> Kpca :: centerTestKernel(const cv::Mat_<float> & K_tst,
                                         const cv::Mat_<float> & K_trn)
{
    // 1_M'
    float *data = new float[K_tst.rows*K_tst.cols];
    cv::Mat_<float> U_tst(K_tst.rows, K_tst.cols, data);
	for (int y = 0; y < U_tst.rows; ++y)
		for (int x = 0; x < U_tst.cols; ++x)
			U_tst[y][x] = 1.0/num_samples;

//	memset(data, 1.0/num_samples, K_tst.rows*K_tst.cols*sizeof(float));
    // 1_M
    cv::Mat_<float> U_trn(num_samples, num_samples, 1.0/num_samples);

    cv::Mat result = K_tst - U_tst*K_trn - K_tst*U_trn + U_tst*K_trn*U_trn;
    delete[] data;
    return result;
}


void Kpca :: computeFeatureVec()
{
	// don't need to make the input-data zero-mean as we do
	// that in feature space

	// establish training-matrix
	cv::Mat_<float> A_trn = selectTrainingSamples();

	// compute Training-Kernel-matrix
	// K_trn = num_samples x num_samples
	cv::Mat_<float> K_trn = getTrainingKernel(A_trn);
	// center Training-Kernel
	K_trn = centerKernel(K_trn);

	// diagonalize K -> eigenvalues, eigenvectors
	cv::Mat evecs;
	cv::Mat evals;
	cv::eigen(K_trn, evals, evecs);

	if (Execution::verbosity >= 3){
		std::cout << "Eigenvalues: ";
		for (int i = 0; i < evals.rows; i++){
			std::cout << evals.at<float>(i,0) << " ";
		}
		std::cout << std::endl;
	}

	// reduce dimensionality if possible
	if (dim == 0){
		dim = PCA::getDimension(evals, eps, config.minFeatureSize, config.maxFeatureSize);
	}
	if (dim < evals.rows){
		evals = evals.rowRange(0,dim).clone();
		evecs = evecs.rowRange(0,dim).clone();
	}

	// normalize eigenvectors
	for (int i = 0; i < evecs.rows; i++){
		cv::Mat r = evecs.row(i);
		r /= std::sqrt(evals.at<float>(i,0));
	}

	// check if we produced nan, if so, we remove the whole row
	// from the evecs-matrix
	dim = PCA::checkNan(evecs);
	// --------------------------------

	// create test gram (or kernel) matrix using A_tst and A_trn
	// A_tst = num-test x num-samples
	cv::Mat_<float> K_tst = getTestKernel(A_trn);

	// center K_tst
	K_tst = centerTestKernel(K_tst, K_trn);

	// project all vectors of K_tst to the eigenvectors (kernel principal
	// components)
	allocateFeature(dim);
	feature.setTo(0.0);

	for (size_t t = 0; t < blocks.size(); t++){
		for (int k = 0; k < dim; k++){
			for (int i = 0; i < num_samples; i++)
			{
				feature(t,k) += K_tst(t,i) * evecs.at<float>(k,i);

				// test of nan -> if nan: flann of opencv231 takes forever
//				if (K_tst(t,i) != K_tst(t,i)){
//					std::cerr << "K_tst(" << t << "," << i << ") is nan" << std::endl;
//					std::cerr << " sigma was: " << sigma << " and num_samples: " << num_samples << std::endl;
//					exit(1);
//				}
//				if (evecs.at<float>(k,i) != evecs.at<float>(k,i)){
//					std::cerr << "evecs(" << k << "," << i << ") is nan" << std::endl;
//					std::cerr << "Evecs( " << k << ", ...): ";
//					for (int e = 0; e < evecs.cols; e++){
//						std::cerr << evecs.at<float>(k,e) << " ";
//					}
//					std::cerr << std::endl;
//					std::cerr << " sigma was: " << sigma << " and num_samples: " << num_samples << std::endl;
//					exit(1);
//				}
//				if (feature(t,k) != feature(t,k)){
//					std::cerr << "feature(" << t << "," << k << ") is nan" << std::endl;
//					std::cerr << " sigma was: " << sigma << " and num_samples: " << num_samples << std::endl;
//					exit(1);
//				}
			}
		}
	}

	// this is the same as the for loops above but w.o. own memory-handling
//	feature = K_tst * evecs.t();

}
