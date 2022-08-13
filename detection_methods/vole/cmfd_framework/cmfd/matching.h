/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef MATCHING_H
#define MATCHING_H

#include <stdexcept>
#include <map>
#include <vector>
#include <set>
#include "cv.h"

/*! 
 * this class tries to match similar feature-vectors 
 */

class Matching {
	public:
		/*!
		 * \param r our feature vectors which will be analyzed
		 */
		Matching(const std::vector<cv::KeyPoint> & _keypoints,
				 const cv::Mat & _feature,
				 const cv::Size & _img_size,
				 struct mat_cfg config);
		/// destructor
		~Matching(void);

		/// copmutes the similar Blocks and stores them in vector of shift vecs and correspondence matrix 
		void computeMatching();
		
		/// returns the index of flann
		inline cv::flann::Index* getIndex(void) const { return ind; }

		/// checks distance and similarity criterion as specified in the config 
		/// returns true if distance smaller and similarity ok
		bool checkDistAndCrit(const cv::Mat_<float> & first, const cv::Mat_<float> & second,
								cv::Point & pa, cv::Point & pb,
								int *dxp, int *dyp);
		
		/// Matching Criterion: euclidian Distance
		bool euclidianDistanceCriterion(const cv::Mat_<float> & first, const cv::Mat_<float> & second) const;
		/// Matching: Criterion: Luo (= vector of thresholds)
		bool LuoCriterion(const cv::Mat_<float> & first,const  cv::Mat_<float> & second) const;
		/// use the maximum of the phase correlation as criterion
		bool cpsCriterion(const cv::Mat_<float> & first,const  cv::Mat_<float> & second) const;
		/// criterions of bravo-Solorio
		bool BravoCriterion(const cv::Mat_<float> & first,const  cv::Mat_<float> & second) const;
		/// use the correlation coefficient
		bool correlCriterion(const cv::Mat_<float> & first,const  cv::Mat_<float> & second) const;
	
		/// computes Chebyshev distance between too block positions and compares with the threshold th
		bool chebyshevDistance(int dx, int dy) const;
		/// computes Euclidian distance between too block positions and compares with the threshold th
		bool euclidianDistance(int dx, int dy) const;


		// -- Helper functions for the correspondence map
		inline cv::Point2f getCorres(int x, int y, int c) const{
			return cv::Point2f(corres[c](y,x).x, corres[c](y,x).y);
		}
		inline std::vector<cv::Mat_<cv::Point> > & getCorresMatrix(void){ 
			return corres;
		}
		static cv::Point2f getCorres(const cv::Mat_<cv::Point> & cor, int x, int y){
			return cv::Point2f(cor(y,x).x, cor(y,x).y);
		}
		static cv::Point2f getCorres(const std::vector<cv::Mat_<cv::Point> > & cor, int x, int y, int c){
			return cv::Point2f(cor[c](y,x).x, cor[c](y,x).y);
		}
		static cv::Point getCorres(const std::vector<cv::Mat_<cv::Point> > & cor, const cv::Point & p, int c){
			return cv::Point2f(cor[c](p.y,p.x).x, cor[c](p.y,p.x).y);
		}
	private:
		/// inserts the corresponding blockpositions in a matrix
		void putBlocksInPointMatrix(void);
		// mini class to compare features' row indices
		class LessThanIdx {
			public:
				LessThanIdx(const cv::Mat_<float> & _features)
					: features(_features){}
				bool operator()(int row_a, int row_b) const{
					for (int x = 0; x < features.cols; x++){
						if (features(row_a, x) < features(row_b, x))
							return true;
						if (features(row_a, x) > features(row_b, x))
							return false;
					}
					// if both are completly equal, has to be false - strictly weak order!
					return false; 
				}
				const cv::Mat_<float> & features;
		};
		/// matching ann with kd tree - helper function
		void computeMatchingKD(void);
		/// matching with lexicographic sorting
		void computeMatchingLexi(void);
		/// matching nearest neighbor with brute-force
		void computeNN(void);
		/// checks the minDistance and the similarity, if ok it puts in the blocks in bset and shiftvec
		void checkAndPut(const cv::Mat_<float> & first, const cv::Mat_<float> & second, 
						cv::Point p_first, cv::Point p_second);
		/// checks entropy of the feature matrix at a certain column
		/// Note: atm only used for Bravo-Solorio et al
		void checkEntropy(int column);
		
		/// the keypoints = positions of features
		const std::vector<cv::KeyPoint> & keypoints;
		/// matrix of features, every row is a feature
		const cv::Mat & features;
		const cv::Size & img_size;
		/// corresponds to the indices matrix - if row is marked as false, the feature and thus
		/// its position will be ignored
		cv::Mat_<uchar> mask;
		/// configuration struct (view copymoveframework.h)
		struct mat_cfg config;
		/// euclidian threshold
		double euc_th;
		/// Index of flann
		cv::flann::Index *ind;
		/// the criterions in one flag variable
		int crit;
		/// datastructures for fast access to corresponding points
		std::vector<cv::Mat_<cv::Point> > corres;
		std::map<int,int> critToId;
};
#endif
