/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include <boost/thread.hpp>
#include <boost/bind.hpp>

#include "execution.h"
#include "feature.h"

Feature::Feature(struct ft_cfg cfg,
				 const BlockHandling & _blocks,
				 int num_threads)
	: blocks(_blocks)
{
	CV_Assert(!blocks.empty());
	block_size = blocks.getBlockSize();
	config = cfg;
	this->num_threads = num_threads;
	feature_data = NULL;
}

Feature :: ~Feature(void){
	// delete internal used data for feature-mattrix
	if (feature_data != NULL){
		delete[] feature_data;
		feature_data = NULL;
	}
}

void Feature :: computeOne(const cv::Mat & cur, int row)
{
	// convert to float
	cv::Mat_<float> cur_f = cur;
	cur_f = cur_f.reshape(1,1);
	cv::Mat_<float> r = feature.row(row);
	cur_f.copyTo(r);
}

void Feature :: computeRange(int start, int end)
{
	for (int i = start; i < end; i++) {
		computeOne(blocks[i], i);
	//	if (i % 1000 == 0)
	//		std::cerr << ".";
	}
}

void Feature :: computeFeatureVec()
{
	// if NOTHING is chosen, really nothing happens (usefull for CPS)
	if (config.featvec == NOTHING) {
		return;
	}
	if (config.featvec == NO){
		allocateFeature(block_size*block_size);
	}
	if (feature.empty()){
		std::cerr << "WARNING: feature hasn't been allocated!"
				  << "-> Call allocateFeature()! You could get strange results\n";
	}
	if (Execution::verbosity >=1){
		std::cout << "->Compute Feature vectors\n";
	}
	// boundings of threads
	int len = blocks.size() / num_threads;
	int rest = blocks.size() % num_threads;
	int begin = len + rest; 

	boost::thread_group thrd_grp;
	// first thread computes the overhead additionally
	boost::thread *t = new boost::thread(boost::bind(
				&Feature::computeRange, this, 0, begin));
	thrd_grp.add_thread(t);

	for (int i = 1; i < num_threads; i++){
		t = new boost::thread(boost::bind(
					&Feature::computeRange, this, begin + (i-1)*len, begin + i*len ));
		thrd_grp.add_thread(t);
	}
	thrd_grp.join_all();
	if (Execution::verbosity >=1){
		std::cout << "<-Compute Feature vectors\n";
	}
}

void Feature :: allocateFeature(int cols)
{
	if (Execution::verbosity >= 2){
		std::cout << "Allocate feature with " << cols << " cols"
				  << " and " << blocks.size() << " rows" << std::endl;
	}
	if (!feature.empty() || blocks.empty() || cols < 0){
		std::cerr << "WARNING: Feature::allocateFeature(): !feature.empty() || blocks.empty() || cols < 0\n";
		std::cerr << "-->return";
		return;
	}
	// if we already have features and want to free the memory
	// call allocateFeature with 0
	if (cols == 0 && feature_data != NULL){
		delete[] feature_data;
		feature_data = NULL;
		return;
	}
	//feature = cv::Mat_<float>(blocks.size(), cols);
	try {
		feature_data = new float[blocks.size() * cols];
	}
	// just give some more output if bad_alloc exception happens,
	// but throw it further
	catch(std::bad_alloc ex){
		std::cerr << "Feature::allocateFeature: bad alloc error\n";
		throw std::bad_alloc(ex);
	}
	feature = cv::Mat_<float>(blocks.size(), cols, feature_data);
}

void Feature :: free(void) {
	allocateFeature(0);
}

void Feature :: printFeature(int row) const
{
	printFeature(this->feature, row);
}

void Feature :: printFeature(const cv::Mat_<float> &feature, int row)
{
	for (int x = 0; x < feature.cols; x++){
		std::cout << feature(row, x) << " ";
	}
	std::cout << std::endl;
}

void Feature :: normalize()
{
	normalize(this->feature);
}

void Feature :: normalize(cv::Mat & feature)
{
	for (int i = 0; i < feature.rows; i++){
		cv::Mat row = feature.row(i);
		cv::Mat dest;
		// normalize with the L2-norm
		cv::normalize(row, dest);
		dest.copyTo(row);
	}
}

// quantize our features
void Feature :: quantize()
{
	int qStart = config.qStart;
	int qEnd = config.qEnd;
	double q = static_cast<double>(config.qNum);

	if (qStart == qEnd) {
		if (Execution::verbosity >= 1 )
			std::cerr << "No quantization chosen\n";
		return;
	}

	if (Execution::verbosity >= 1 ){
		std::cerr << "quantize Now\n";
	}

	if (qEnd == -1){
		if (Execution::verbosity >= 1 )
			std::cerr << "quantize(): rows x cols " << feature.rows
					  << " x " << feature.cols << std::endl;
		qEnd = feature.cols;
	}
	if (qEnd < qStart) {
		throw std::invalid_argument("Feature::quantize: qEnd < qStart");
	}
	if ( feature.cols < qStart || qStart < 0) {
		throw std::invalid_argument("Feature::quantize qStart is out of range");
	}
	if (feature.cols < (qEnd - qStart)){
		throw std::invalid_argument("Feature::quantize: act->featureSize < num qEnd - qStart");
	}

	for (int y = 0; y < feature.rows; y++){
		for (int x = qStart; x < qEnd; x++){
			double val = feature(y, x);
			// take absolute value
			if (config.qAbs) {
				val = fabs(val);
			}
			// round it down
			if (config.qRound){
				val = floor(val);
			}
			// quantize
			val /= q;
			// round it down
			if (config.qRound2){
				val = floor(val);
			}
			feature(y, x) = val;
		}
	}
}
\

double Feature :: optiM(int width)
{
	return ( (width) / ( log(width) - log(sqrt(2)) ) );
}

