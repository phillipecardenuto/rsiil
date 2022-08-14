/*
	Copyright(c) 2012 Christian Riess <christian.riess@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef VOLE_CMFD_COMMAND_REVERSE_ENGINEER_SNIPPET_H
#define VOLE_CMFD_COMMAND_REVERSE_ENGINEER_SNIPPET_H

#include "command.h"
#include "reverse_engineer_snippet_config.h"

namespace vole {
namespace cmfdgt {

class CommandReverseEngineerSnippet : public vole::Command {
public:
	CommandReverseEngineerSnippet();

	int execute();

	void printHelp() const;
	void printShortHelp() const;

	ReverseEngineerSnippetConfig config;
private:

};

}
}
#endif // VOLE_CMFD_REVERSE_ENGINEER_SNIPPET_H
