/*
	Copyright(c) 2012 Christian Riess <christian.riess@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef VOLE_GROUND_TRUTH_H__
#define VOLE_GROUND_TRUTH_H__

#include <iostream>
#include "command.h"
#include "ground_truth_config.h"

namespace vole {

class GroundTruth : public Command
{
public:
	GroundTruth();

	int execute();

	void printHelp() const;
	void printShortHelp() const;

	GroundTruthConfig config;
private:

};

}

#endif // VOLE_GROUND_TRUTH_H__
