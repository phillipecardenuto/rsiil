/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef COPYMOVEFRAMEWORK_H
#define COPYMOVEFRAMEWORK_H

#include <iostream>
#include "command.h"
#include "copymoveframework_config.h"
#include <opencv2/core/core.hpp>

class CopyMoveFramework : public vole::Command {
public:
	CopyMoveFramework(void);
	~CopyMoveFramework(void) { /*delete i;*/ }
	/// execution with saving
	int execute(void);
	/// execution without saving
	cv::Mat* execute_headless(void);
	
	void printShortHelp(void) const;
	void printHelp(void) const;

	void printConfig(void);
	
	// this is our config struct. please keep it public, baby!
	CmfdConfig config;
	// temporary image 
    // TODO: pointer to Mat is unnecessarily complicated.
    // (consider inlining 'execute_headless')...
	cv::Mat *i_img;
private:
	/// safed old time
	double old_time;
	/// updates the time and returns relative time
	inline double relativeTime(){
		double t = static_cast<double>(cv::getTickCount());
		double elapsed = (t - old_time) / cv::getTickFrequency();
		old_time = t;
		return elapsed;
	}
};

#endif // COPYMOVEFRAMEWORK_H
