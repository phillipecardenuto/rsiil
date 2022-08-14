/*
 * Original work Copyright (c) 2016, Thibaud Ehret <ehret.thibaud@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef FEATMANAGER_HEADER_H 
#define FEATMANAGER_HEADER_H 

#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include "../LibImages.h"

/**
 *
 * Basis for different feature manager
 *
 */

class FeatManager 
{
	public:
		/** 
		 * @brief Compute the squared euclidean distance between the first feature and the second feature
		 * 
		 * @param id1: Index of the first feature
		 * @param id2: index of the second feature
		 *
		 * @return distance between the two patches
		 */
		virtual float distance(unsigned id1, unsigned id2) = 0;

		/**
		 * @brief Compute the index of the neighboring patch in the specified direction
		 *
		 * @param feature: Query feature
		 * @param direction: direction for the new patch. 0: x+1. 1: y+1. 2: x-1. 3: y-1.
		 *
		 * @return the index of the resulting feature
		 */
		virtual unsigned getNeighboringFeat(unsigned feat, int direction) = 0;

		/**
		 * @brief Compute the index of the neighboring patch in the specified direction (in the inverse direction compared to the previous function)
		 *
		 * @param feat: Query feature
		 * @param direction: direction for the new feature. 0: x-1. 1: y-1. 2: x+1. 3: y+1.
		 *
		 * @return the index of the resulting feature
		 */
		virtual unsigned getInverseNeighboringFeat(unsigned feat, int direction) = 0;

		/**
		 * @brief Accessor for the size of the feature map
		 *
		 * @return the size of the feature map
		 */
		virtual ImageSize sz() = 0;

	private:
};
#endif
