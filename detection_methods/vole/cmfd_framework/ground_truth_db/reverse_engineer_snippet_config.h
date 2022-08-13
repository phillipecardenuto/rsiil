/*
	Copyright(c) 2012 Christian Riess <christian.riess@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef VOLE_CMFD_REVERSE_ENGINEER_SNIPPET_CONFIG_H
#define VOLE_CMFD_REVERSE_ENGINEER_SNIPPET_CONFIG_H

#include <iostream>
#include "vole_config.h"

namespace vole { namespace cmfdgt {

class ReverseEngineerSnippetConfig : public vole::Config {
public:	
	ReverseEngineerSnippetConfig(const std::string &prefix = std::string());

	//! input file name
	std::string orig_file;
	//! file name of full image including the tampered region
	std::string copy_file;
	//! file name of manipulated region
	std::string snippet_file;
	//! alpha threshold for the overlay
	std::string snippet_alpha_file;
	//! alpha threshold for the overlay (whole image)
	std::string full_snippet_alpha_file;

	virtual std::string getString() const;

	protected:
	#ifdef WITH_BOOST_PROGRAM_OPTIONS
		virtual void initBoostOptions();
	#endif // WITH_BOOST_PROGRAM_OPTIONS
};

} }

#endif // VOLE_CMFD_REVERSE_ENGINEER_SNIPPET_CONFIG_H
