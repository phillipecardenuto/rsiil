/*
	Copyright(c) 2012 Christian Riess <christian.riess@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef FILESYSTEM_HELPERS_H
#define FILESYSTEM_HELPERS_H

#include <sys/stat.h>
#include <iostream>

class FilesystemHelpers {
public:

	static bool file_exists(std::string filename);

	/* FIXME double-check: function is maybe not portable */
	static std::string strip_last_filename_component(std::string file);

	static bool recursive_mkdir(std::string path);

	// TODO function looks buggy; better test it
	static std::string basename_without_extension(std::string file);

};




#endif // FILESYSTEM_HELPERS_H
