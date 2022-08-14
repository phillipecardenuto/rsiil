/*
	Copyright(c) 2012 Christian Riess <christian.riess@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "reverse_engineer_snippet_config.h"

using namespace boost::program_options;

namespace vole { namespace cmfdgt {

ReverseEngineerSnippetConfig::ReverseEngineerSnippetConfig(const std::string &prefix) : Config(prefix) {
	#ifdef WITH_BOOST_PROGRAM_OPTIONS
		initBoostOptions();
	#endif // WITH_BOOST_PROGRAM_OPTIONS
}


std::string ReverseEngineerSnippetConfig::getString() const {
	std::stringstream s;

	if (prefix_enabled)
		s << "[" << prefix << "]" << std::endl; 

	s	<< "orig=" << orig_file << " # Image to process" << std::endl
		<< "copy=" << copy_file << 
			   " # image containing the copy-move forgery" << std::endl
		<< "snippet=" << snippet_file << 
			   " # manipulated region" << std::endl
		<< "alpha_snippet=" << snippet_alpha_file << 
			   " # snippet alpha channel file" << std::endl
		<< "full_alpha_snippet=" << full_snippet_alpha_file << 
			   " # snippet alpha channel file in the whole image dimensions" << std::endl
		;
	return s.str();
}

#ifdef WITH_BOOST_PROGRAM_OPTIONS
void ReverseEngineerSnippetConfig::initBoostOptions() {
	options.add_options()
		(key("orig"), value(&orig_file)->default_value(""),
			"Image to process")
		(key("copy"), value(&copy_file)->default_value(""),
			"image containing the copy-move forgery")
		(key("snippet"), value(&snippet_file)->default_value(""),
			"manipulated region")
		(key("alpha_snippet"), value(&snippet_alpha_file)->default_value(""),
			"snippet alpha channel file")
		(key("full_alpha_snippet"), value(&full_snippet_alpha_file)->default_value(""),
			   "snippet alpha channel file in the whole image dimensions")
	;
}
#endif // WITH_BOOST_PROGRAM_OPTIONS



} }
