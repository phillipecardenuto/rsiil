/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include <iostream>
#include <cmath>

#include <boost/thread.hpp>
#include <boost/bind.hpp>

#include "blur.h"
static const int fak[] = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320};

// we only need binom coefficients up to 7, 
// so do it this way:
int Blur::binom(unsigned int n, unsigned int k) const {
	if (k > n)
		return 0;
	if (k == 0 || k == n) return 1;

	return fak[n] / (fak[k] * fak[n-k]);
}


double Blur::mue(unsigned int p, unsigned int q, cv::Mat_<double> &act) {
	// speed-up - look in table if already calculated
	if (mueTable[p][q] != -1){
		return mueTable[p][q];
	}
	double m = 0;
	cv::Mat_<uchar> mask(act.rows, act.cols, static_cast<uchar>(0));
	if (config.circle){ 
		cv::circle(mask, cv::Point(act.cols/2, act.rows/2), 
				static_cast<uchar>(act.rows/2), cv::Scalar(1,0), -1);
	}
	for (int j = 1; j <= act.rows; j++){
		for (int i = 1; i <= act.cols; i++){
			if (config.circle && mask(j,i) == 0)
				continue;
			m += pow((double)(i - xt),(double) p) 
				* pow( (double)(j - yt), (double)q)
				* act(j-1, i-1);
		}
	}
	mueTable[p][q] = m;
	return m;
}

double Blur::blur(unsigned int p, unsigned int q, cv::Mat_<double> &act) {
	// !modified -> axis symmetric -> 24 moments
	
	// if order == even 
	// (Flusser, Degraded Image Analysis: An invariant approach, Appendix A
	if (config.modified) {
		if ( (p + q) % 2 == 0) {
			return 0;
		}
		// if p == q then the blur invariant is equal to a lower order one
		if ( p == q ) {
			return 0;
		}
	}
	else {
		if ( (p == q ) && (p % 2 == 0) ) {
			return 0;
		}
	}

	if ( ((p == 0) && (q == 1)) || ((p == 1) && (q == 0)) ){
		return 0;
	}
	// another speed-up -> look in table if already computed:
	if (blurTable[p][q] != -1){
		return blurTable[p][q];
	}
	unsigned long long b = 0;
	for (unsigned int n = 0; n <= p; n++){
		for (unsigned int m = 0; m <= q; m++){
			if (!config.modified){
				if ( (n + m) % 2 == 0){
					continue;
				}
			}
			if ( (n + m <= 0) || (n + m >= p +q) ){
				continue;
			}
			double rek = blur(p - n, q - m, act);
			if (rek == 0) { // little speed-up
				continue;
			}
			b += static_cast<double>(binom(p, n)) * 
				static_cast<double>(binom(q, m)) *
				rek *
				mue(n, m, act);
		}
	}

	double mue00 = mue(0, 0, act);
	double c = mue(p, q, act);
	
	// this will never happen, 
	// cause we dont use 28 invariants (would be for radially symmetric)
	if (!config.modified && (p % 2 == 0) && (q % 2 == 0)){
		c -= mue(q, p, act);
	}

	if (mue00 != 0.0){ // prevent possible division by zero
		c -= (b / mue00);
	}
	else {
		c -= b;
	}

	cv::Mat_<uchar> mask(act.rows, act.cols, static_cast<uchar>(0));
	if (config.circle){ 
		cv::circle(mask, cv::Point(act.cols/2, act.rows/2), 
				static_cast<uchar>(act.rows/2), cv::Scalar(1,0), -1);
	}

	// normalize // blockSize * blockSize better? contradiction of mahdian with flusser here
	double n;
	if (config.circle)
		n = pow( (double)(cv::sum(mask)[0] * 0.5), (double)(p+q) ) * mue00;
	else
		n = pow( ((double)block_size * 0.5), (double)p + q ) * mue00;
	if (n != 0.0){ // prevent possible division by zero
		c /= n;
	}

	blurTable[p][q] = c;

	return c;
}

void Blur::computeCentroids(cv::Mat_<double> &chan){
	double m00 = 0.0;
	double m01 = 0.0;
	double m10 = 0.0;
	cv::Mat_<uchar> mask(chan.rows, chan.cols, static_cast<uchar>(0));
	if (config.circle){
		// create mask
		cv::circle(mask, cv::Point(chan.cols/2, chan.rows/2), 
				static_cast<uchar>(chan.rows/2), cv::Scalar(1,0), -1);
		for (int y = 0; y < chan.rows; y++){
			for (int x = 0; x < chan.cols; x++){
				m00 += chan(y,x);
			}
		}
	}
	else
		m00 = static_cast<double>(cv::sum(chan)[0]); // compute sum over our channel

	for (int y = 1; y <= chan.rows; y++){
		for (int x = 1; x <= chan.cols; x++){
			if (config.circle && mask(y,x) == 0)
				continue;
			double val = chan(y - 1, x - 1);
			m01 += y * val;
			m10 += x * val;
		}
	}
	
	// we need these!
	mueTable[0][0] = m00;
//	cv::Moments m = moments(chan);

	mueTable[0][1] = 0.;
	mueTable[1][0] = 0.;
	

	if (m00 != 0.0){ // prevent division by zero
		xt = m10 / m00;
		yt = m01 / m00;
	}
	else {
		xt = 0.0;
		yt = 0.0;
	}
}

void Blur :: initializeTables(void){
	for (int y = 0; y < 7; y++){
		for (int x = 0; x < 7; x++){
			blurTable[x][y] = -1;
			mueTable[x][y] = -1;
		}
	}
}

void Blur :: computeRange(int start, int end)
{
	for (int i = start; i < end; i++) {
		computeOne(blocks[i], i);
	}
}

void Blur :: computeFeatureVec(void)
{
	if (config.modified){
		allocateFeature(18);
	}
	else { 
		allocateFeature(24);
	}
	// boundings of threads
	int len = blocks.size() / num_threads;
	int rest = blocks.size() % num_threads;
	int begin = len + rest; 

	boost::thread_group thrd_grp;
	// first thread computes the overhead additionally
	boost::thread *t = new boost::thread(boost::bind(
				&Blur::computeRange, this->clone(), 0, begin));
	thrd_grp.add_thread(t);

	for (int i = 1; i < num_threads; i++){
		t = new boost::thread(boost::bind(
					&Blur::computeRange, this->clone(), begin + (i-1)*len, begin + i*len ));
		thrd_grp.add_thread(t);
	}
	thrd_grp.join_all();
}

void Blur :: computeOne(const cv::Mat & block, int row)
{
	CV_Assert(feature.data != NULL);
	CV_Assert(block.channels() == 1);
	// convert to double
	cv::Mat_<double> cur = block;

	initializeTables();
	computeCentroids(cur);

	if (config.modified){
		feature(row, 0) = blur(3, 0, cur);
		feature(row, 1) = blur(2, 1, cur);
		feature(row, 2) = blur(1, 2, cur);
		feature(row, 3) = blur(0, 3, cur);

		feature(row, 4) = blur(5, 0, cur);
		feature(row, 5) = blur(4, 1, cur);
		feature(row, 6) = blur(3, 2, cur);
		feature(row, 7) = blur(2, 3, cur);
		feature(row, 8) = blur(1, 4, cur);
		feature(row, 9) = blur(0, 5, cur);

		feature(row, 10) = blur(7, 0, cur);
		feature(row, 11) = blur(6, 1, cur);
		feature(row, 12) = blur(5, 2, cur);
		feature(row, 13) = blur(4, 3, cur);
		feature(row, 14) = blur(3, 4, cur);
		feature(row, 15) = blur(2, 5, cur);
		feature(row, 16) = blur(1, 6, cur);
		feature(row, 17) = blur(0, 7, cur);
	}
	else 
	{
		feature(row, 0) = mue(1, 1);

		feature(row, 1) = blur(1, 2, cur);
		feature(row, 2) = blur(2, 1, cur);
		feature(row, 3) = blur(0, 3, cur);
		feature(row, 4) = blur(3, 0, cur);

		feature(row, 5) = blur(1, 3, cur);
		feature(row, 6) = blur(3, 1, cur);
		feature(row, 7) = blur(3, 2, cur);
		feature(row, 8) = blur(2, 3, cur);
		feature(row, 9) = blur(4, 1, cur);
		feature(row, 10) = blur(1, 4, cur);

		feature(row, 11) = blur(0, 5, cur);
		feature(row, 12) = blur(5, 0, cur);
		feature(row, 13) = blur(3, 3, cur);
		feature(row, 14) = blur(1, 5, cur);
		feature(row, 15) = blur(5, 1, cur);
		feature(row, 16) = blur(0, 7, cur);
		feature(row, 17) = blur(1, 6, cur);
		feature(row, 18) = blur(2, 5, cur);
		feature(row, 19) = blur(3, 4, cur);
		feature(row, 20) = blur(7, 0, cur);
		feature(row, 21) = blur(6, 1, cur);
		feature(row, 22) = blur(5, 2, cur);
		feature(row, 23) = blur(3, 4, cur);
	}  
}
