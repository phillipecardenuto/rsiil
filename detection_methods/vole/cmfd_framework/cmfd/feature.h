/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef FEATURE_H
#define FEATURE_H

#include "opencv2/opencv.hpp"
#include "copymoveframework_config.h"
#include "blockhandling.h"

/// This is the feature extraction basis class
class Feature {
	public:
		//Feature(){}
		/// Constructor 
        Feature(struct ft_cfg cfg, const BlockHandling & _blocks, int num_threads = -1);
		virtual ~Feature(void);
		virtual void computeOne(const cv::Mat & current_block, int row);
		virtual void computeRange(int start, int end);
		virtual void computeFeatureVec(void);
		/// quantize the featurevectors in bins
		void quantize(void);
        // needed by FMT, TODO: maybe move fnct to that class
		static double optiM(int width);
		// Note: would work if we dont have maxFeatureSize
//		const cv::Mat_<float> & getFeature() const{
//			return feature;
//		}
		cv::Mat_<float> & getFeature(){
			if (config.maxFeatureSize > 0 && feature.cols > config.maxFeatureSize){
				// here we change from self-handled data to opencv-handled memory
				// as we assume that the feature-size will be much smaller now
				feature = feature(cv::Range::all(), cv::Range(0, config.maxFeatureSize)).clone();
				delete feature_data;
				feature_data = NULL;
				return feature;
			}
			else {
				return feature;
			}
			//return feature;
		}
		/// normalizes features by the L2-norm
		void normalize(void);
		/// static variant of normaliing feature with the L2 norm
		static void normalize(cv::Mat & feature);
		/// print a feature vector
		void printFeature(int row) const;
		/// static variant of printing a feature vector
		static void printFeature(const cv::Mat_<float> & feature, int row);
		/// only a subclass of Feature may call 'allocateFeature'
		/// but you may free here the internally used data
		/// !use with care!
		void free(void);
	protected:
		/// you have to call that function!
		void allocateFeature(int cols);
		/// matrix of feature vectors
		cv::Mat_<float> feature;
		/// our feature struct
		struct ft_cfg config;
		/// reference to our blocks
		const BlockHandling & blocks;
		/// size of every block
		int block_size;
		/// number of threads
		int num_threads;
	private:
		/// internal used data-array for the feature-matrix
		float *feature_data;
};

#endif
