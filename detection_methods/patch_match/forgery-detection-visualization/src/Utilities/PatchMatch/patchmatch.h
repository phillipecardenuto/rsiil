/*
 * Original work Copyright (c) 2017, Thibaud Ehret <ehret.thibaud@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PATCHMATCH_HEADER_H
#define PATCHMATCH_HEADER_H

#include <vector>
#include <algorithm> 
#include <random>
#include <unordered_map>
#include "../FeatManager/featManager.h"
#include "../LibImages.h"
#include "../Utilities.h"
#include "../parameters.h"
#ifdef _OPENMP
#include <omp.h>
#endif

class Patchmatch
{
	public:
		/**
		 * @brief Construct a PatchMatch object on which matching will be done
		 *
		 * @param fm = FeatureManager on which PatchMatch will be applied
		 matrix;
		 * @param th = the (spatial) distance threshold below which the matching are discarded 
		 **/
		Patchmatch(FeatManager& fm, int th);

		/**
		 * @brief Initialize the displacement map at random (must still be valid matches)
		 *
		 * @return 0 if no problem, -1 otherwise
		 **/
		int initializeRandom();

		/**
		 * @brief Do one step of PatchMatch propagation
		 **/
		void propagate();

		/**
		 * @brief Do one step of PatchMatch backpropagation
		 **/
		void backpropagate();

		/**
		 * @brief Do N steps of PatchMatch propagations (alternating between actual propagations and backpropagations) with a random search between each propagations
		 * @param N = the number of propagations and backpropagations to be made
		 **/
		void propagateNtimes(int N);

		/**
		 * @brief Extract the displacement maps (both along the x and y axis) associated to the matches computed by PatchMatch
		 *
		 * @param dispX = for each patch, correspond to the displacement along the x axis between the patch and its match
		 * @param dispY = for each patch, correspond to the displacement along the y axis between the patch and its match
		 **/
		void extract(std::vector<int>& dispX, std::vector<int>& dispY);

		/**
		 * @brief Copy the distance associated to the matches
		 *
		 * @param dists = vector in which the result is stored
		 **/
		void dists(std::vector<float>& dists) 
		{
			dists.resize(imSize.wh);
			for(int i = 0; i < imSize.wh; ++i)
				dists[i] = distances[i];
		};

		/**
		 * @brief Apply a random search PatchMatch step to the current displacement map
		 **/
		void randomSearch();

		/*
		* @assign the value of matches to a vector
		*/
		void get_detected_matches(std::vector<int>& copy, std::vector<bool>&detectionMask){
				copy.assign(matches.begin(),matches.end());
				for(int i=0; i<copy.size();i++){
					if(!detectionMask[i]){
						copy[i] =0;
					}
				}
		}
	private:
		/**
		 * @brief Propagate to currentId the match from direction using a zero-order propagation
		 *
		 * @param currentId = index of the patch where the propagation is done
		 * @param direction = direction in which the propagation is done:		 
		 * 	0 = left
		 * 	1 = top
		 * 	2 = top-left
		 * 	3 = top-right
		 **/
		inline void propToId(unsigned currentId, int direction);

		/**
		 * @brief Propagate to currentId the match from direction using a first-order propagation
		 *
		 * @param currentId = index of the patch where the propagation is done
		 * @param direction = direction in which the propagation is done:		 
		 * 	0 = left
		 * 	1 = top
		 * 	2 = top-left
		 * 	3 = top-right
		 **/
		inline void propToId1(unsigned currentId, int direction);

		/**
		 * @brief Symmetrize the match of id if it is better
		 *
		 * @param id = index of the patch which match is to be symmetrized
		 **/
		inline void symmetrize(unsigned id);

		/// Manager for the features to do the patchmatch on
		FeatManager* fm;

		/// match map for the different nearest neighbors
		std::vector<int> matches;
		/// Distances associated to the matches
		std::vector<float> distances;

		ImageSize imSize;

		int th;
};
#endif
