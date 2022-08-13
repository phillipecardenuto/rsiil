/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "log.h"
#include <iostream>

#include <opencv/cv.h>

Log :: Log()
{
	start_time = static_cast<double>(cv::getTickCount());
	filename = "/tmp/cmfd.log";
	log_level = 1;
}

Log :: Log(const std::string & filename,
		   int log_level,
		   bool log_time)
{
	set(filename, log_level, log_time);
}

bool Log :: fexists(const char *filename)
{
	std::ifstream ifile(filename);
	return ifile;
}

Log :: ~Log()
{
	if (log_level > 0){
		std::ofstream file;
		// as we have boost::filesystem anyway as lib we can use it to check if the file
		// already exists, then append
		if ( fexists(filename.c_str()) ){
			file.open( filename.c_str(), std::ios::out | std::ios::app );
		}
		else {
			file.open( filename.c_str() );
		}
		file << os.str();
		if (log_time){
			double t = static_cast<double>(cv::getTickCount());
			file << "Elapsed time: " << (t-start_time)/cv::getTickFrequency() << "\n";
		}
		file.close();
	}
}

void Log :: set(const std::string &_filename,
				int _log_level,
				bool _log_time)
{
	filename = _filename;
	log_level = _log_level;
	log_time = _log_time;
}

std::ostringstream & Log :: operator()(int mode)
{
	if (mode <= log_level){
		if (log_time){
			double t = static_cast<double>(cv::getTickCount());
			// print elapsed time since start
			os << (t - start_time) / cv::getTickFrequency() << ": ";
		}
		return os;
	}
	else {
		dummy.clear();
		// probably there exist sth more clever
		return dummy;
	}
}
