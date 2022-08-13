/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "cmfd_util.h"
namespace cmfd
{
	std::string getBasename(std::string inputfile)
	{
		std::string image_basename;
		size_t min_index = inputfile.find_last_of("/\\");
		size_t max_index = inputfile.find_last_of('.');
		if (min_index == std::string::npos) { min_index = 0; } else { min_index++; }
		if (max_index == std::string::npos) {
			image_basename = inputfile.substr(min_index, std::string::npos);
		} else {
			if (max_index == std::string::npos)
				image_basename = inputfile.substr(min_index, std::string::npos);
			else
				image_basename = inputfile.substr(min_index, max_index-min_index);
		}

		return image_basename;
	}
}
