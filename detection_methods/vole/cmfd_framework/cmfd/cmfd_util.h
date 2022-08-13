/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef CMFD_UTIL_H
#define CMFD_UTIL_H

#include <string>

namespace cmfd
{
	// computes the basename of the image
	// Note: qt-library has sth similar too
	std::string getBasename(std::string inputfile);
}

#endif
