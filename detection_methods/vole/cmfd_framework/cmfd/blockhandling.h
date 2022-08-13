/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef BLOCKHANDLING_H
#define BLOCKHANDLING_H

#include <vector>
#include "cv.h"
#include "copymoveframework_config.h"

class BlockHandling {
	public:
		BlockHandling(cv::Mat& i, struct bl_cfg cfg);
		/// Copy Constructor // flat copy of underlying matrices;
		BlockHandling(const BlockHandling &a);
		~BlockHandling();
		/// sets the blocks from the image
		void computeBlocks(void);
		/// return reference to the internal vector where the blocks are saved
		const std::vector<cv::Mat> & getBlocks(void) const{
			return blocks;
		}
		/// keypoints
		const std::vector<cv::KeyPoint> & getKeypoints(void) const {
			return keypoints;
		}
		/// returns true if there are no blocks in the vector
		inline bool empty(void) const {
			return blocks.empty();
		}
		/// returns number of blocks
		inline size_t size(void) const {
			return blocks.size();
		}
		/// returns the blocksize of one block
		inline int getBlockSize(void) const {
			return config.blockSize;
		}
		inline const cv::Mat & operator[](unsigned int i) const{
			return blocks[i];
		}

		/// assigns Blockhandling class to another
		BlockHandling& operator=(const BlockHandling &a);
		
		/// returns the position of a block in the image
		static cv::Point getPosi(const cv::Mat & a){
			cv::Point p;
			cv::Size s;
			a.locateROI(s, p);
			return p;
		}
	
		/// computes a log polar representation of a block
		/// and sums it up to a 1-D representation
		static cv::Mat logPolar(const cv::Mat & a, bool center, int deg, int dstwidth);
		static void printBlock(const cv::Mat & a);
		/// checks entropy of the blocks and removes them if they are lower than entropyTh
		void checkEntropy(int threadnum, double entropyTh, bool circle);
	private:	
		//struct cmp{
		//bool operator()(const cv::Mat& a) const {
		//} cmpEntropy;
		double getEntropy(cv::Mat mat, bool circle) const;
		void checkEntropyRange(int start, int end, bool circle);
		/// vector of pointers to blocks
		std::vector<cv::Mat> blocks;
		std::vector<cv::KeyPoint> keypoints;
		/// vectors of entropies
		std::vector<double> entropies;
		/// reference to our image
		cv::Mat& img;
		/// our config (step, blockSize, roi)
		struct bl_cfg config;
		/// helper functio for adding blocks to the block vector
		void addBlocks(unsigned int startX, unsigned int startY, unsigned int endX, unsigned int endY);
};
#endif
