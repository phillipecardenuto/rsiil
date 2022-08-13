/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "svd.h"
#include "execution.h"
#include <cassert>
#include "cv.h"

#include <boost/thread.hpp>
#include <boost/bind.hpp>

Svd :: Svd(struct ft_cfg _config, const BlockHandling & blocks, int numThreads)
    :	Feature(_config, blocks, numThreads)
{
    // if we don't wanna variable dimension and nobody has set the dimension
    // compute one from the first block - this is not recommended
    if (config.dim == 0){
        cv::Mat_<float> cur = blocks[0];

        if (config.circle){
            cv::Mat_<double> mask(cur.rows, cur.cols, 0.0);
            cv::circle(mask, cv::Point(cur.cols/2, cur.rows/2),
                       static_cast<int>(cur.rows/2), cv::Scalar(1.0,0.0), -1);
            cur *= mask;
        }

        // if dimension set to 0 we chose the dimension as maximal fraction of the SVs
        cv::SVD s(cur, cv::SVD::NO_UV); // dont need u and v here
        cv::Mat_<double> w = s.w;
        globalDim = computeDim(w);

        if ( Execution::verbosity >= 1 ){
            std::cerr << "SVD: new dimension (first block-first channel) = " << globalDim << std::endl;
        }
        allocateFeature(globalDim);
    }
    // if dim  == -1 take the average dimension
    else if (config.dim == -1) {
        size_t dim = 0;
        for (size_t i = 0; i < blocks.size(); i++)
        {
            cv::Mat cur = blocks[i];
            cv::Mat_<double> chan = cur; // convert to double

            if (config.circle){
                cv::Mat_<double> mask(chan.rows, chan.cols, 0.0);
                cv::circle(mask, cv::Point(chan.cols/2, chan.rows/2),
                           static_cast<int>(chan.rows/2), cv::Scalar(1.0,0.0), -1);
                chan *= mask;
            }

            cv::SVD s(chan, cv::SVD::NO_UV); // dont need u and v here
            cv::Mat_<double> w = s.w;

            dim += computeDim(w);
        }
        globalDim = dim / blocks.size();
        allocateFeature(globalDim);
        if ( Execution::verbosity >= 1 ){
            std::cerr << "SVD: new dimension (avg dim) = " << globalDim << std::endl;
        }
    }
    else if (config.dim == -2){
        // this is completly arbitrary
        max_dim = 32;
        allocateFeature(max_dim);
        feature.setTo(0.0);
    }
    else{
        globalDim = config.dim;
        allocateFeature(globalDim);
        if ( Execution::verbosity >= 1 )
            std::cerr << "SVD: new dimension is computed individually for every block or is set by dim\n";
    }
}


// get new dimension, if it's set to 0
int Svd :: computeDim(cv::Mat_<double> &w)
{
    // new dimension of block i
    int dimi;
    // sum over all evals
    double evalsum = cv::sum(w)[0];

    double crit = 1.0 - config.eps;
    double sum = 0.0;
    int cnt = 0;

    // determine new dimension
    for (int y = 0; y < w.rows; y++) {
        cnt++;
        double tmp = w.at<double>(y,0);
        if (tmp == 0.0){
            break;
        }
        sum += tmp;
        if (crit <= sum / evalsum){
            break;
        }
    }

    dimi = cnt;

    if (dimi == w.rows && Execution::verbosity >=2){
        std::cerr << "Warning: SVD failed: New Dimension == Old Dimension - No useful dimension could be found with eps: " << config.eps << std::endl;
    }

    return dimi;
}

void Svd :: project(int row, cv::SVD &s, int dim)
{
    // redue the rank
    cv::Mat_<double> tmpu = s.u(cv::Range(0,dim), cv::Range::all());
    cv::Mat_<double> wdiag = cv::Mat_<double>::zeros(s.u.cols,s.u.rows);
    cv::Mat_<double> tmpw = s.w;
    for (int p = 0; p < dim; p++){
        wdiag(p,p) = tmpw(p,0);
        //	feature(0,p) = tmpw(p,0);
    }
    cv::Mat_<double> tmpvt = s.vt(cv::Range::all(), cv::Range(0, dim));

    // reduced rank approximation
    cv::Mat_<double> dec = tmpu * wdiag * tmpvt;

    cv::Mat_<float> tmp;
    // for the approach of dwt + svd -> no second SVD
    if (config.modified){
		//	dec = dec.reshape(1,1);
        //	std::cerr << feature.cols << " " << feature.rows << " dim: " << dim<< std::endl;
        //	feature = dec;
		tmp = tmpw;
    }
    else {
        // second SVD computation -> features = singular values
        cv::SVD svd(dec, cv::SVD::FULL_UV);
        tmp = svd.w;
    }
	tmp = tmp.reshape(1,1);
    for (int i = 0; i < dim; i++){
        feature(row, i) = tmp(0,i);
    }
}

void Svd :: computeOne(const cv::Mat & cur, int row)
{
    int dim;
    cv::Mat_<double> chan = cur; // convert to double
    if (config.circle){
        cv::Mat_<double> mask(chan.rows, chan.cols, 0.0);
        cv::circle(mask, cv::Point(chan.cols/2, chan.rows/2),
                   static_cast<int>(chan.rows/2), cv::Scalar(1.0,0.0), -1);
        chan *= mask;
    }

    cv::SVD s(chan, cv::SVD::FULL_UV);

    if (config.dim == -2){
        cv::Mat_<double> w = s.w;
        dim = computeDim(w);
        dim = std::min(dim, max_dim);
    } else {
        dim = globalDim;
    }

    project(row, s, dim);

}

