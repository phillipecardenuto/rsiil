/*
	Copyright(c) 2012 Christian Riess <christian.riess@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef VOLE_SPLICE_CONFIG_RANDOM_BLOCK_H
#define VOLE_SPLICE_CONFIG_RANDOM_BLOCK_H

#include "vole_config.h"

#include <iostream>
#include <vector>


/**
 * Configuration parameters for the ground truth generation.
 */
namespace vole {

class SpliceConfigRandomBlock : public Config {
public:	
	SpliceConfigRandomBlock(const std::string &prefix = std::string());

	//! Original (untampered) image file
	std::string orig_file;
	//! directory for all intermediate files
	std::string output_directory;
	//! prefix for the output file name
	std::string file_prefix;
	//! Output image shall be compressed as jpeg
	bool output_is_jpeg;
	/// region of interest of input file
	std::string rio;

	// number of blocks to paste: rand[min,max]
	int num_pastes_min;
	int num_pastes_max;
	/// block area size = block_size * block_size
	int block_size;

	bool no_ground_truth;
	bool ground_truth_is_snippet;
	bool ground_truth_is_cmfd;

	//! add noise only on snippet
	bool local_noise;
	//! add noise afterwards
	bool global_noise;

	//! add rotation
	bool rot;
	//! add scaling
	bool scale;
    bool global_cont_scale;
    bool global_discrete_scale;

	//! add jpeg artifacts per snippet
	bool local_jpeg;
	//! add jpeg artifacts on final image
	bool global_jpeg;

	//! noise level; Standard deviation of Gaussian noise (for snippet)
	int local_noise_std_dev; // divided by 1000
	int global_noise_std_dev; // divided by 1000

	// rotation section
	int angle;  // in degrees

	// scaling section
	int scaling; // in per mille (1/1000)
	int global_cont_scale_factor;
	int global_discrete_scale_factor;

	// jpeg section
	int local_jpeg_quality;
	int global_jpeg_quality;

	// correct the insert positions for global scaling
	bool correct_insert_positions;
	
	virtual std::string getString() const;

	protected:
	#ifdef WITH_BOOST_PROGRAM_OPTIONS
		virtual void initBoostOptions();
	#endif // WITH_BOOST

};

}


#endif // VOLE_SPLICE_CONFIG_H
