/*
	Copyright(c) 2012 Christian Riess <christian.riess@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef VOLE_GROUND_TRUTH_SUBIMAGE_H
#define VOLE_GROUND_TRUTH_SUBIMAGE_H

#include "ground_truth_common.h"

/** Encapsulates an image or a snippet as a combination of the color
 * information and the transparency information. Note that the transparency
 * mask is processed with morphologic operators that are thresholding at
 * opaque_intensity. Values higher than opaque_intensity are considered opaque
 * pixels, the others as transparent pixels.
 */
class SubImage {
public:

	SubImage(const cv::Mat_<cv::Vec3b> colors, const cv::Mat_<unsigned char> mask, int opaque_intensity);
	void show() const;


	cv::Mat_<cv::Vec3b> colors;
	cv::Mat_<unsigned char> mask, thrmask;
	cv::Mat_<cv::Vec3b> maskBGR;
	int opaque_intensity;
};

#endif // VOLE_GROUND_TRUTH_SUBIMAGE_H
