/*
	Copyright(c) 2012 Christian Riess <christian.riess@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef VOLE_CMFD_GT_POST_PROCESSING_H
#define VOLE_CMFD_GT_POST_PROCESSING_H

#include "command.h"
#include "gt_post_processing_config.h"

namespace vole { namespace cmfdgt {

class CommandGtPostProcessing : public Command
{
public:
	CommandGtPostProcessing();

	int execute();

	void printHelp() const;
	void printShortHelp() const;

	GtPostProcessingConfig config;
};

}}

#endif // VOLE_CMFD_GT_POST_PROCESSING_H
