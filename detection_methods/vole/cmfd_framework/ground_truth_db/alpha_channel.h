/*
	Copyright(c) 2012 Christian Riess <christian.riess@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef VOLE_GT_ALPHA_CHANNEL_H
#define VOLE_GT_ALPHA_CHANNEL_H

#include "cv.h"

namespace vole {
	namespace cmfdgt {
		// computes a bounding box around all not fully transparent pixels,
		// from the image boundary to inside of the image
		cv::Rect computeBoundingBox(cv::Mat_<unsigned char> full_alpha_image);
	}
}

#endif // VOLE_GT_ALPHA_CHANNEL_H
