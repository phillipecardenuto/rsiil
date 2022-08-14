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

#ifndef ORIPM_HEADER_H
#define ORIPM_HEADER_H

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <complex>
#include <cmath>
#include <iostream>
#include "../parameters.h"
#include "featManager.h"
#include "../LibImages.h"

#define PI 3.1415926535

/**
 *
 * Create a manager for Zernike moments
 *
 */

class ZMManager: public FeatManager {
	public:
		/**
		 * @brief Basic constructor
		 *
		 * @param image: The image
		 * @param dsz: Size of the image 
		 * @param order: Maximum order for the moments 
		 * @param numRho: Number of radial sampling 
		 * @param numTheta: Number of angle sampling 
		 * @param zm: Type of Zernike moments to compute 
		 * @param flip: Test again the flipped image 
		 */
		ZMManager(const std::vector<float>& image, ImageSize dsz, int sp, int order, int numRho, int numTheta, int zm, bool flip);


		/**
		 * @brief Compute the Zernike moments for a given image
		 *
		 * @param image: The image
		 * @param dsz: Size of the image 
		 * @param order: Maximum order for the moments 
		 * @param numRho: Number of radial sampling 
		 * @param numTheta: Number of angle sampling 
		 * @param zm: Type of Zernike moments to compute 
		 */
		void computeZM(const std::vector<float>& image, ImageSize dsz, int sp, int order, int numRho, int numTheta, int zm, std::vector<std::vector<float> >& features);

		/**
		 * check featManager.h for information on these functions
		 * @{
		 */
		float distance(unsigned id1, unsigned id2);

		unsigned getNeighboringFeat(unsigned feat, int type);
		unsigned getInverseNeighboringFeat(unsigned feat, int type);
		ImageSize sz() {return descSize;};

		std::vector<std::vector<float> > exportFeatures() {return features1;};
		
		~ZMManager();
	private:
		/// The feature map on which the manager is based
		std::vector<std::vector<float> > features1;
		/// The feature map on which the manager is based for the flipped version
		std::vector<std::vector<float> > features2;
		/// Size of the feature map
		ImageSize descSize;
};
#endif
