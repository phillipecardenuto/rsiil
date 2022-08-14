/*
	Copyright(c) 2012 Christian Riess <christian.riess@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef VOLE_GT_POST_PROCESSING_CONFIG_H
#define VOLE_GT_POST_PROCESSING_CONFIG_H

#include "vole_config.h"

namespace vole { namespace cmfdgt {

class GtPostProcessingConfig : public Config {
public:
	GtPostProcessingConfig(const std::string &prefix = std::string());

	//! ground truth file for post-processing
	std::string ground_truth_file;

	//! post-processed ground truth file (i.e. the output)
	std::string pp_ground_truth_file;

	int support_size;
	int support_ratio;

	virtual std::string getString() const;

	protected:
	#ifdef WITH_BOOST
		virtual void initBoostOptions();
	#endif // WITH_BOOST


};

}}

#endif // VOLE_GT_POST_PROCESSING_CONFIG_H

