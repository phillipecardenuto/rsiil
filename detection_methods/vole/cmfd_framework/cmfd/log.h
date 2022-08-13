/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef LOG_H
#define LOG_H

#include <iosfwd>
#include <fstream>
#include <sstream>

// very very simple log class
// not threadsafe!

class Log{
public:
	Log();
	Log(const std::string & filename,
		int log_level = 1,
		bool log_time = false);
	/// writes out the stream
	~Log();
	void set(const std::string & filename,
			 int log_level = 1,
			 bool log_time = false);
	/// gives stream which is used for logging
	std::ostringstream & operator()(int mode = 1);
private:
	/// check if file already exists
	bool fexists(const char * filename);
	/// file-name to log into
	std::string filename;
	/// only calls with modes <= log_level will
	/// be logged
	int log_level;
	/// the starting time
	double start_time;
	/// gives the elapsed time in seconds
	bool log_time;
	/// dummy stringstream will be cleared every time
	std::ostringstream dummy;
	/// real stringstream, will be written to file
	/// when the destructor is called
	std::ostringstream os;
};

#endif // LOG_H
