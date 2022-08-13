/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include <stdexcept>
#include <iostream>
#include "blockhandling.h"
#include "feature.h"
#include "execution.h"
#include "copymoveframework.h"

#include <boost/thread.hpp>
#include <boost/bind.hpp>
BlockHandling::
BlockHandling(cv::Mat &_img, struct bl_cfg cfg)
	: img(_img)
{
	config = cfg;
	if (config.step < 1){
		throw std::runtime_error("BlockHandling: wrong stepsize");
	}
	if (config.blockSize < 1){
		throw std::runtime_error("BlockHandling: wrong blockSize");
	}
	if (config.roi.size() % 4 != 0) {
		throw std::runtime_error("Blockhandling: wrong amount of ROI-params - 4 for every ROI");
	}
}

BlockHandling::BlockHandling(const BlockHandling &a)
	:	img(a.img)
{	
	config = a.config;
	blocks = a.blocks;
}

// the destructor
BlockHandling :: ~BlockHandling(){}

void BlockHandling::addBlocks(unsigned int startX, unsigned int startY, unsigned int endX, unsigned int endY){
	int dx = endX - startX;
	int dy = endY - startY;
	// number of blocks
	int sx = (dx - config.blockSize + config.step) / config.step;
	int sy = (dy - config.blockSize + config.step) / config.step;
	if (sx <= 0 || sy <= 0){
		std::cout << "WARNING: cant generate enough blocks with"
				  << " block size: " << config.blockSize
				  << " step: " << config.step
				  << " width: " << dx << " height: " << dy << std::endl;
		return;
	}
	if (Execution::verbosity >= 2){
		std::cout << "reserve for " << blocks.size() + sx*sy << " blocks\n";
	}
	// save some time and reserve space
	blocks.reserve(blocks.size() + sx*sy);
	keypoints.reserve(keypoints.size()+ sx*sy);

	for (int y = startY; y + config.blockSize <= endY; y += config.step){ // for every row
		for (int x = startX; x + config.blockSize <= endX; x += config.step){ // for every column
			// cut a special range respectivly a block out of the image
			cv::Mat tmp = img( cv::Range(y, y + config.blockSize),
					cv::Range(x, x + config.blockSize) );
			blocks.push_back(tmp);
			cv::KeyPoint key(cv::Point2f(x,y), config.blockSize);
			// add shift
//			cv::KeyPoint key( cv::Point2f(x+config.blockSize/2, y+config.blockSize/2),
//							  config.blockSize );
			keypoints.push_back(key);
		}
	}
}

void BlockHandling::computeBlocks(void ) {
	unsigned int startX;
	unsigned int startY;
	int endX;
	int endY;

	if (config.roi.empty()){
		startX = 0;
		startY = 0;
		endX = img.cols;
		endY = img.rows;

		addBlocks(startX, startY, endX, endY);
	}
	else {
		for (unsigned int n = 0; n < config.roi.size(); n+=4)
		{
			// some checks of corners of the diagonal 
			// left up and right down
			if (config.roi[n] < config.roi[n+2] && config.roi[n+1] < config.roi[n+3]) {
				startX = config.roi[n];
				startY = config.roi[n+1];
				endX = config.roi[n+2];
				endY = config.roi[n+3];
			}
			// right down and left up
			else if (config.roi[n] > config.roi[n+2] && config.roi[n+3] > config.roi[n+3]) {
				startX = config.roi[n+2];
				startY = config.roi[n+3];
				endX = config.roi[n];
				endY = config.roi[n+1];
			}
			// right up and left down
			else if (config.roi[n] > config.roi[n+2] && config.roi[n+3] < config.roi[n+3]) {
				startX = config.roi[n+2];
				startY = config.roi[n+1];
				endX = config.roi[n];
				endY = config.roi[n+3];
			}
			// left down and right up
			else if (config.roi[n] < config.roi[n+2] && config.roi[n+3] > config.roi[n+3]) {
				startX = config.roi[n];
				startY = config.roi[n+3];
				endX = config.roi[n+2];
				endY = config.roi[n+1];
			}
			else {
				throw std::invalid_argument("Blockhandling: wrong config.roi-params");
			}
			// prevent blocks out of img
			if (endX >= img.cols){
				endX = img.cols;
			}
			if (endY >= img.rows){
				endY = img.rows;
			}

			addBlocks(startX, startY, endX, endY);
		}
	}
}

BlockHandling& BlockHandling::operator=(const BlockHandling &a)
{
	if (this != &a) {
		img = a.img;
		blocks = a.blocks;
		config = a.config;
	}
	return *this;
}

// a has to be single chan!
cv::Mat BlockHandling :: logPolar(const cv::Mat &a, bool center, int deg, int dstwidth) 
{
	if (a.channels() != 1) 
		throw std::runtime_error("BlockHandling::logPolar a.channels != 1");
	
	// determine optimal M (= manipulator for logpolar coord trafo)
	double optiM = Feature::optiM(dstwidth);
	// prepare data
	cv::Mat_<float> dat = a;  // convert to float
//	dat = dat.reshape(a.data.rows);
	CvMat val(dat);
	cv::Mat_<float> logpolar(deg, dstwidth);
	CvMat logp(logpolar);

	// compute log polar coords
	if (center)
		cvLogPolar( &val, &logp, 
				cvPoint2D32f(dat.cols/2, dat.rows/2), 
				optiM, CV_WARP_FILL_OUTLIERS);
	else
		cvLogPolar( &val, &logp, 
				cvPoint2D32f(0,0),
				optiM, CV_WARP_FILL_OUTLIERS);

	// 1-D representation
	cv::Mat dst;
	cv::reduce(logpolar, dst, 0, CV_REDUCE_SUM);

	return dst;
}

double BlockHandling :: getEntropy(cv::Mat mat, bool circle) const
{
	// circular blocks -> create mask
	cv::Mat_<uchar> mask;
	if (circle){
		mask = cv::Mat_<uchar>(mat.rows, mat.cols, static_cast<uchar>(0));
		cv::circle(mask, cv::Point(mat.cols/2, mat.rows/2), 
				static_cast<int>(mat.rows/2), cv::Scalar(1,0), -1);
	}

	// get probabilities of lumi
	int chan[1]; // take only the Y -channel - have only this one
	chan[0] = 0;
	int histSize[] = {256}; // take 256 bins -> should we change that?
	float lumirange[] = {0, 256};
	const float* ranges[] = {lumirange};
	cv::MatND hist; // output histogram
	cv::calcHist(&mat, 1, chan, mask, // use mask
			hist, 1, histSize, ranges, 
			true, // histogram is uniform
			false); // don't akkumulate

	double entropy = 0.0;
	float div;
	if (circle)
		div = cv::sum(mask)[0];
	else 
		div = mat.rows * mat.cols;
	for (int y = 0; y < 256; y++){
		// compute luminance probability 
		float p = hist.at<float>(y) / static_cast<float>(div); 		
		if (p == 0.0) 
			continue;
		// compute entropy with log to basis log_2
		entropy -= p * (log(p) / log(2));
	}
	return entropy;
}

void BlockHandling :: checkEntropyRange(int start, int end, bool circle)
{
	for (int i = start; i < end; i++){
		cv::Mat_<float> gray;
		if (blocks[i].channels() == 1)
			gray = blocks[i];
		else {
			CV_Assert(blocks[i].channels() == 3);
			cv::cvtColor(blocks[i], gray, CV_BGR2GRAY);
		}
		entropies[i] = getEntropy(gray, circle);
	}
}

// checks entropy in multiple threads
void BlockHandling :: checkEntropy(int threadnum, double entropyTh, bool circle) 
{
	if (Execution::verbosity >= 2){
		std::cout << "-> BlockHandling::checkEntropy\n";
	}
	entropies.resize(blocks.size());
	
	// boundings of threads
	int len = blocks.size() / threadnum;
	int rest = blocks.size() % threadnum;
	int begin = len + rest; 

	boost::thread_group thrd_grp;
	// first thread computes the overhead additionally
	boost::thread *t = new boost::thread(boost::bind(
				&BlockHandling::checkEntropyRange, this, 0, begin, circle));
	thrd_grp.add_thread(t);

	for (unsigned int i = 1; i < (unsigned int) threadnum; i++){
		t = new boost::thread(boost::bind(
					&BlockHandling::checkEntropyRange, this, begin + (i-1)*len, begin + i*len, circle ));
		thrd_grp.add_thread(t);
	}
	thrd_grp.join_all();

	// remove blocks with too low entropy
	// this doesn't work anymore
	//lowEntropy.entropy = entropyTh;
	//blocks.erase(std::remove_if(blocks.begin(), blocks.end(), lowEntropy), blocks.end());

	// why is that so much slower than the line above?
	for( int i = 0; i < (int)blocks.size(); i++ ){
		if (entropies[i] < entropyTh){
			blocks.erase(blocks.begin() + i);
			keypoints.erase(keypoints.begin() + i);
		}
	}
	if (Execution::verbosity >= 2){
		std::cout << "<- BlockHandling::checkEntropy\n";
	}
}

void BlockHandling::printBlock(const cv::Mat & a) {
	std::cerr << "--------Data of Block-----------\n";
	std::cerr << "x = " << (getPosi(a)).x << std::endl;
	std::cerr << "y = " << (getPosi(a)).y << std::endl;
	std::cerr << "channels = " << a.channels() << std::endl;
	std::cerr << "size = " << a.cols << "x" << a.rows << std::endl;
	if (a.empty()){
		std::cerr << "data is empty\n";
		return;
	}
	
	for (int y = 0; y < a.rows; y ++){
		for (int x = 0; x < a.cols; x ++){
			std::vector<cv::Mat> planes;
			split (a, planes);
			for (int c = 0; c < a.channels(); c ++){
				std::cerr << ((int) planes[c].at<uchar>(y,x));
				if (c < a.channels() - 1)
					std::cerr << " ";
			}
			std::cerr << "|";
		}
		std::cerr << std::endl;
	}
	std::cerr << std::endl;
}
