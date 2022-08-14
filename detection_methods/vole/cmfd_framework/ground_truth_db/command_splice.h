/*
	Copyright(c) 2012 Christian Riess <christian.riess@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef VOLE_COMMAND_SPLICE_H__
#define VOLE_COMMAND_SPLICE_H__

#include <iostream>
#include "command.h"
#include "splice_config.h"
#include "cv.h"

namespace vole { namespace cmfd_gt {

/** Splices an image and a number of snippets and generates the attached ground
 * truth. It applies additionally rotation, scaling, Gaussian noise and JPEG
 * artifacts, depending on the configuration values.
 */
class CommandSplice : public vole::Command
{
public:
	CommandSplice();

	int execute();

	void printHelp() const;
	void printShortHelp() const;

	SpliceConfig config;

private:
//	void parseOptions(void);
	std::pair<std::string, std::string> get_output_file_names();
	void correctInsertPositions(std::vector<int> &source_positions, std::vector<int> & insert_positions, cv::Mat &orig);

	std::vector<int> insert_positions;
	std::vector<cv::Mat_<cv::Vec3b> > snippets;
	std::vector<cv::Mat_<uchar> > snippets_alpha;
};

} }

#endif // VOLE_COMMAND_SPLICE_H__
