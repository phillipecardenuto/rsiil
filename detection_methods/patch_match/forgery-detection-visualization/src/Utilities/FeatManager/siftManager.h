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

#ifndef SIFTPM_HEADER_H
#define SIFTPM_HEADER_H

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <string>
#include "featManager.h"
#include "../LibImages.h"
extern "C" {
    #include "vl/generic.h"
    #include "vl/dsift.h"
    #include "vl/imopv.h"
}

/**
 *
 * Create a manager for sift features
 *
 */

class SiftManager: public FeatManager {
	public:
		/**
		 * @brief Basic constructor
		 *
		 * @param descr: The feature map
		 * @param sz: Size of the feature map 
		 */
		SiftManager(const std::vector<float>& image, ImageSize dsz, int binSize, bool flip);

		/**
		 * check featManager.h for information on these functions
		 * @{
		 */
		float distance(unsigned id1, unsigned id2);

		unsigned getNeighboringFeat(unsigned feat, int type);
		unsigned getInverseNeighboringFeat(unsigned feat, int type);
		ImageSize sz() {return descSize;};
	private:
		/// The feature map on which the manager is based
		std::vector<float> features1;
		/// The feature map on which the manager is based for the flipped version
		std::vector<float> features2;
		/// Size of the feature map
		ImageSize descSize;
};
#endif
