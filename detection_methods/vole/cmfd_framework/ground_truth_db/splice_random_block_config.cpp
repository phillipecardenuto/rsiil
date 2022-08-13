/*
	Copyright(c) 2012 Christian Riess <christian.riess@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "splice_random_block_config.h"

using namespace boost::program_options;

namespace vole {

SpliceConfigRandomBlock::SpliceConfigRandomBlock(const std::string &prefix) : Config(prefix) {
	#ifdef WITH_BOOST_PROGRAM_OPTIONS
		initBoostOptions();
	#endif // WITH_BOOST
}

std::string SpliceConfigRandomBlock::getString() const {
	std::stringstream s;

	if (prefix_enabled)
		s << "[" << prefix << "]" << std::endl; 

	s	<< "local_noise=" << (local_noise ? 1:0) << " # add noise only to the snippet" << std::endl
		<< "global_noise=" << (global_noise ? 1:0) << " # add noise to the whole image" << std::endl
		<< "rot=" << (rot ? 1:0) << " # add rotation" << std::endl
		<< "scale=" << (scale ? 1:0) << " # add scaling" << std::endl
		<< "local_jpeg=" << (local_jpeg ? 1:0) << " # add JPEG compression per snippet" << std::endl
		<< "global_jpeg=" << (global_jpeg ? 1:0) << " # add JPEG compression on the final image" << std::endl
		<< "orig=" << orig_file << " # Original (untampered) image file" << std::endl
		<< "output=" << output_directory << " # Working directory" << std::endl
		<< "output_is_jpeg=" << (output_is_jpeg ? 1:0) << " # Add JPEG compression to final image"
		<< "file_prefix=" << file_prefix << " # prefix for the output file name" << std::endl
		<< "num_blocks_min=" << num_pastes_min << " # lower bound of numbers of blocks"
		<< "num_blocks_max=" << num_pastes_max << " # upper bound of numbers of blocks"
		<< "block_size=" << block_size << " # length of one edge of the block -> area = block_size * block_size"
		<< "rio=" << rio << " # region of interest (e.g. 0 0 800 600)"
		<< "no_ground_truth=" << (no_ground_truth ? 1:0) << " # only add postprocessing to original images, without splicing snippets" << std::endl
		<< "ground_truth_is_snippet=" << (ground_truth_is_snippet ? 1:0) << " # use only snippet for ground truth, not source region" << std::endl
		<< "ground_truth_is_cmfd=" << (ground_truth_is_cmfd ? 1:0) << " # use snippet and source region for ground truth" << std::endl
		<< "local_noise_std_dev=" << local_noise_std_dev << 
			   " # Standard deviation of Gaussian noise if noise should be added to the snippet" << std::endl
		<< "global_noise_std_dev=" << global_noise_std_dev << 
			   " # Standard deviation of Gaussian noise if noise should be added to the whole_image" << std::endl
		<< "correct_insert_positions=" << (correct_insert_positions ? 1:0) << " correct the insert position " << std::endl
		<< "angle=" << angle << " # angle for rotation of the manipulated region" << std::endl

		// params
 		<< "lnoise_stddev=" << local_noise_std_dev <<
			" # Standard deviation in per mille (1/1000) of Gaussian noise (added to the snippet)" << std::endl
 		<< "gnoise_stddev=" << global_noise_std_dev <<
			" # Standard deviation in per mille (1/1000) of Gaussian noise (added to the whole image)" << std::endl
		<< "scaling=" << scaling <<
			" # scaling factor of the snippet, per mille (1/1000)" << std::endl
		<< "cont_scale_factor=" << global_cont_scale_factor <<
			" # global continous scaling factor, per mille (1/1000)" << std::endl
		<< "discrete_scale_factor=" << global_discrete_scale_factor <<
			" # global discrete scaling factor, per mille (1/1000)" << std::endl
		<< "angle=" << angle <<
			" # angle for rotation of the manipulated region" << std::endl
		<< "gjquality=" << global_jpeg_quality <<
			" # quality factor for jpeg compression on whole image" << std::endl
		<< "ljquality=" << local_jpeg_quality <<
			" # quality factor for jpeg compression per snippet" << std::endl
		;
	return s.str();
}


#ifdef WITH_BOOST_PROGRAM_OPTIONS
void SpliceConfigRandomBlock::initBoostOptions() {
	options.add_options()
		// operations
		(key("local_noise"), bool_switch(&local_noise)->default_value(false), "add noise only to the snippet")
		(key("global_noise"), bool_switch(&global_noise)->default_value(false), "add noise to the whole image")
		(key("rot"), bool_switch(&rot)->default_value(false), "add rotation")
		(key("scale"), bool_switch(&scale)->default_value(false), "add scaling")
		(key("global_cont_scale"), bool_switch(&global_cont_scale)->default_value(false), "add global continous scaling")
		(key("global_discrete_scale"), bool_switch(&global_discrete_scale)->default_value(false), "add global discrete scaling")
		(key("local_jpeg"), bool_switch(&local_jpeg)->default_value(false), "add JPEG compression per snippet")
		(key("global_jpeg"), bool_switch(&global_jpeg)->default_value(false), "add JPEG compression on whole image")

		// files
		(key("orig"), value(&orig_file)->default_value(""), "Original (untampered) image file")
		(key("output,O"), value(&output_directory)->default_value("/tmp/"), "Directory where the image is written to ")
		(key("output_is_jpeg"), bool_switch(&output_is_jpeg)->default_value(false), "Add JPEG compression to final image")
		(key("file_prefix"), value(&file_prefix)->default_value(""), "Prefix for the output file name")
		(key("no_ground_truth"), bool_switch(&no_ground_truth)->default_value(false),
			"only add postprocessing to original images, without splicing snippets")
		(key("ground_truth_is_snippet"), bool_switch(&ground_truth_is_snippet)->default_value(false),
			" use only snippet for ground truth, not source region")
		(key("ground_truth_is_cmfd"), bool_switch(&ground_truth_is_cmfd)->default_value(false),
			" use snippet and source region for ground truth")
		(key("rio"), value(&rio)->default_value(""),
			"region of interest of input file")

		// snippet source and insertion positions
		(key("num_pastes_min"), value(&num_pastes_min)->default_value(2),
			 "number of pastes of one random block - lower bound (num_blocks = rand[min,max])")
		(key("num_pastes_max"), value(&num_pastes_max)->default_value(2),
			 "number of pastes of one random block - upper bound (num_blocks = rand[min,max])")
		(key("block_size"), value(&block_size)->default_value(64),
			 "one dimension of the block (area of block = block_size*block_size")

		(key("correct_insert_positions"), bool_switch(&correct_insert_positions)->default_value(false),
			"Correct the insert positions")

		// params
 		(key("lnoise_stddev"), value(&local_noise_std_dev)->default_value(10),
			"Standard deviation in per mille (1/1000) of Gaussian noise (added to the snippet)")
 		(key("gnoise_stddev"), value(&global_noise_std_dev)->default_value(10),
			"Standard deviation in per mille (1/1000) of Gaussian noise (added to the whole image)")
		(key("scaling"), value(&scaling)->default_value(10),
			"scaling factor of the snippet, per mille (1/1000)")
		(key("cont_scale_factor"), value(&global_cont_scale_factor)->default_value(500),
			"global scaling factor for continuous scaling (= linear interpolation), per mille (1/1000)")
		(key("discrete_scale_factor"), value(&global_discrete_scale_factor)->default_value(500),
			"global scaling factor for discrete scaling (= nearest neighbor interpolation), per mille (1/1000)")
		(key("angle"), value(&angle)->default_value(10), // e.g. 0
			"angle for rotation of the manipulated region")
		(key("gjquality"), value(&global_jpeg_quality)->default_value(90), // e.g. 0
			"quality factor for jpeg compression on whole image")
		(key("ljquality"), value(&local_jpeg_quality)->default_value(90), // e.g. 0
			"quality factor for jpeg compression per snippet")
	;
}
#endif // WITH_BOOST_PROGRAM_OPTIONS

}

