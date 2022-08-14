/*
	Copyright(c) 2012 Christian Riess <christian.riess@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef VOLE_GROUND_TRUTH_COMMON_H
#define VOLE_GROUND_TRUTH_COMMON_H
#include <cv.h>

#define EPSILON 0.00001

typedef cv::Vec3f BGR;
typedef cv::Mat_<BGR> MatBGR;
typedef cv::MatConstIterator_<BGR> ItBGR;

#endif // VOLE_GROUND_TRUTH_COMMON_H
