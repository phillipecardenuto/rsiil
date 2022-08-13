/*
	Copyright(c) 2012 Christian Riess <christian.riess@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "command_splice_random_block.h"

#include "highgui.h"
#include "string_conversion.h"
#include "splice_core.h"

#include <fstream>
#include <sstream>

using namespace boost::program_options;

namespace vole { namespace cmfd_gt {

CommandSpliceRandomBlock::CommandSpliceRandomBlock()
	: Command(
		"splice_rb",
		config,
		"Christian Riess",
		"christian.riess@informatik.uni-erlangen.de")
{
}

int CommandSpliceRandomBlock::execute()
{
	if (config.orig_file.length() < 1) { std::cerr << "Add original image (--orig)!" << std::endl; printHelp(); return 1; }

	if (config.file_prefix.length() < 1) { std::cout << "WARNING: Add file prefix (--file_prefix) to better distinguish files!" << std::endl; }

	if (!(config.ground_truth_is_snippet || config.ground_truth_is_cmfd || config.no_ground_truth)) {
		std::cerr << "either --ground_truth_is_snippet, --ground_truth_is_cmfd or --no_ground_truth must be specified!" << std::endl;
		return 1;
	}

	cv::Mat_<cv::Vec3b> orig = cv::imread(config.orig_file, CV_LOAD_IMAGE_COLOR);
    if (orig.empty()) { std::cerr << "error loading image " << config.orig_file << std::endl; return 1; }


	if (!config.rio.empty()){
		// e.g. 0 0 800 600
		std::vector<int> rio = vole::StringConversion::stringToIntList(config.rio, ',');
		if (rio.size() != 4){
			std::cerr << "if you give a region of interest you have to give 4 point\n";
			return 1;
		}
		orig = orig(cv::Range(rio[1] ,rio[3]), // rows
					cv::Range(rio[0], rio[2])).clone(); // cols
	}
	// initialize randomizer
	srand(time(NULL));


	std::vector<int> source_positions;
	std::vector<int> insert_positions;
	std::vector<cv::Mat_<cv::Vec3b> > snippets;
	std::vector<cv::Mat_<unsigned char> > alpha_snippets;
	// + 1 as upper bound is inclusive
	int num_pastes = rand() % (config.num_pastes_max-config.num_pastes_min + 1) + config.num_pastes_min;
	if (num_pastes == 0){
		std::cerr << "num_blocks = 0 -> wrong min/max ? \n";
		return 1;
	}
	std::cout << "insert " << num_pastes
			  << " random blocks with edge length: "
			  << config.block_size << std::endl;
	// determine random source and insert position
	int s_x = rand() % (orig.cols-config.block_size);
	int s_y = rand() % (orig.rows-config.block_size);
	for (int i = 0; i < num_pastes; i++){
		int i_x = rand() % (orig.cols-config.block_size);
		int i_y = rand() % (orig.rows-config.block_size);
		// prevent same position
		if (s_x == i_x && s_y == s_y){
			i--;
			continue;
		}
		cv::Mat_<cv::Vec3b> snippet = orig(cv::Range(s_y, s_y + config.block_size),
									   cv::Range(s_x, s_x + config.block_size)).clone();
		cv::Mat_<uchar> alpha_snippet(config.block_size, config.block_size, (uchar)255);
		source_positions.push_back(s_x);
		source_positions.push_back(s_y);
		insert_positions.push_back(i_x);
		insert_positions.push_back(i_y);
		snippets.push_back(snippet);
		alpha_snippets.push_back(alpha_snippet);
	}


	// correct insert position
	if (config.correct_insert_positions){
		correctInsertPositions(source_positions, insert_positions, orig);
	}

	SpliceCore splicing(orig, snippets, alpha_snippets, source_positions, insert_positions, 200 /* FIXME make configurable*/);

	if (!config.no_ground_truth) { // perform all per-snippet operations

		if (config.local_noise) {
			for (int i = 0; i < (int) snippets.size(); ++i) {
				splicing.add_snippet_noise(i, config.local_noise_std_dev/1000.0);
			}
		}

		if (config.rot) {
			for (int i = 0; i < (int) snippets.size(); ++i) {
				splicing.add_rotation(i, config.angle, false);
			}
		}

		if (config.scale) {
			for (int i = 0; i < (int) snippets.size(); ++i) {
				splicing.add_scaling(i, config.scaling / 1000.0, false);
			}
		}

		if (config.local_jpeg) {
			for (int i = 0; i < (int) snippets.size(); ++i) {
				splicing.add_snippet_jpeg(i, config.local_jpeg_quality);
			}
		}

		splicing.splice(config.ground_truth_is_cmfd);
	}

	
	// global scaling in 2 different ways (linear vs nearest neighbor interpolation)
	if (config.global_cont_scale) {
			// with linear interpolation
		splicing.add_global_cont_scale(config.global_cont_scale_factor / 1000.0);
	} 
	if (config.global_discrete_scale){
		splicing.add_global_discrete_scale(config.global_discrete_scale_factor / 1000.0);
	}

	if (config.global_noise) {
		splicing.add_global_noise(config.global_noise_std_dev/1000.0);
	}
	
	if (config.global_jpeg && !config.output_is_jpeg) {
		splicing.add_global_jpeg(config.global_jpeg_quality);
	}

	cv::Mat_<cv::Vec3b> tampered = splicing.getTamperedImage();
	cv::Mat_<unsigned char> ground_truth = splicing.getAlphaImage();

	std::pair<std::string, std::string> output_files = get_output_file_names();

	if (!config.output_is_jpeg) {
		cv::imwrite(output_files.first, tampered);
	} else {
		std::vector<int> params;
		params.push_back(CV_IMWRITE_JPEG_QUALITY);
		params.push_back(config.global_jpeg_quality);
		cv::imwrite(output_files.first, tampered, params);
	}

	cv::imwrite(output_files.second, ground_truth);

/*
	cv::Mat_<cv::Vec3b> sanity;
	tampered.copyTo(sanity);
	for (int y = 0; y < sanity.rows; ++y) {
		for (int x = 0; x < sanity.cols; ++x) {
			sanity[y][x][2] = ground_truth[y][x];
		}
	}
	cv::imwrite("/tmp/sanity.png", sanity);
*/
	if (config.verbosity > 0) {
		std::cout << "wrote files." << std::endl
			<< "copy: " << output_files.first << std::endl
			<< "gt:   " << output_files.second << std::endl
			<< "sanity:   " << "/tmp/sanity.png" << std::endl;
	}
	return 0;
}

std::pair<std::string, std::string> CommandSpliceRandomBlock::get_output_file_names() {
	std::stringstream tampered, gt;
	tampered << config.output_directory << "/";
	gt << config.output_directory << "/";
	if (config.file_prefix.length() > 0) {
		tampered << config.file_prefix;
		gt << config.file_prefix;
	} else {
		tampered << "anonymous_forgery";
		gt << "anonymous_forgery_";
	}

	if (!config.no_ground_truth) {
		tampered << "_copy";
	}
	gt << "_gt";

	if (config.num_pastes_min == config.num_pastes_max){
		tampered << "_rb" << config.num_pastes_max;
		gt << "_rb" << config.num_pastes_max;
	}
	else {
		tampered << "_rb" << config.num_pastes_min << "-" << config.num_pastes_max;
		gt << "_rb" << config.num_pastes_min << "-" << config.num_pastes_max;
	}

	if (config.local_noise) {
		tampered << "_ln" << config.local_noise_std_dev;
		gt << "_ln" << config.local_noise_std_dev;
	}

	if (config.rot) {
		tampered << "_r" << config.angle;
		gt << "_r" << config.angle;
	}

	if (config.scale) {
		tampered << "_s" << config.scaling;
		gt << "_s" << config.scaling;
	}
	if (config.global_cont_scale){
		tampered << "_gcs" << config.global_cont_scale_factor;
		gt << "_gcs" << config.global_cont_scale_factor;
	}
	if (config.global_discrete_scale){
		tampered << "_gds" << config.global_discrete_scale_factor;
		gt << "_gds" << config.global_discrete_scale_factor;
	}
	if (config.local_jpeg) {
		tampered << "_lj" << config.local_jpeg_quality;
		gt << "_lj" << config.local_jpeg_quality;
	}

	if (config.global_noise) {
		tampered << "_gn" << config.global_noise_std_dev;
		gt << "_gn" << config.global_noise_std_dev;
	}

	if (config.global_jpeg) {
		tampered << "_gj" << config.global_jpeg_quality;
		gt << "_gj" << config.global_jpeg_quality;
	}

	if (config.output_is_jpeg)
		tampered << ".jpg";
	else
		tampered << ".png";

	gt << ".png";

	return std::pair<std::string, std::string>(tampered.str(), gt.str());
}

void CommandSpliceRandomBlock::correctInsertPositions(std::vector<int> &source_positions, std::vector<int> & insert_positions, cv::Mat& orig){
	int scaling_factor = 1000;
	if (config.global_cont_scale) {
		scaling_factor = config.global_cont_scale_factor;
	}
	if (config.global_discrete_scale) {
		scaling_factor = config.global_discrete_scale_factor;
	}
	if (scaling_factor != 1000){
		double factor = scaling_factor / 1000.0;
		if (factor < 1.0)
			factor = 1.0 / factor;
		// todo: works atm only with scaling factors of precision of 2
		int factor_int = static_cast<int>(factor);
		if ((factor - static_cast<double>(factor_int)) > 0.0)
			factor_int *= 10;
		for (size_t i = 0; i < insert_positions.size(); i++){
			int source_pos = source_positions[i];
			int insert_pos = insert_positions[i];
			int source_rest = source_pos % factor_int;
			int insert_rest = insert_pos % factor_int;
			// fits already -> skip
			if (source_rest == insert_rest)	
				continue;
			int new_insert_pos = -1;
			
			// try to adjust upwards
			if ( std::abs(insert_rest - source_rest) > (factor_int/2) ) {
				int end_iter = ((i % 2) == 0 
						? std::min(orig.cols, insert_pos+factor_int) 
						: std::min(orig.rows, insert_pos+factor_int));
				for (int k = insert_pos+1; k < end_iter; k++){
					if (k == source_pos)
						continue;
					int new_rest = k % factor_int;
					if (source_rest == new_rest){
						new_insert_pos = k;
						break;
					}
				}
			}
			// try to adjust downwards
			else {
				int end_iter = std::max(0, insert_pos - factor_int);
				for (int k = insert_pos-1; k >= end_iter; k--){
					if (k == source_pos)
						continue;
					int new_rest = k % factor_int;
					if (source_rest == new_rest){
						new_insert_pos = k;
						break;
					}
				}
			}
			if (new_insert_pos != -1){ // found a better insertion position
				std::cout << "correct insert position from " << insert_positions[i] 
					<< " to " << new_insert_pos << std::endl;
				insert_positions[i] = new_insert_pos;
			}
		} // end for every insert_position
	} // end if factor != 1000
}		
void CommandSpliceRandomBlock::printHelp() const {
	std::cout << "splices images, creates ground truth automatically" << std::endl;
	std::cout << std::endl;
	std::cout << "This command combines an input image and random blocks" << std::endl;
	std::cout << "to a forgery. Optionally, noise, JPEG compression" << std::endl;
	std::cout << "artifacts, rotation and scaling can be applied. Noise and" << std::endl;
	std::cout << "JPEG artifacts can be applied to the block alone, or to" << std::endl;
	std::cout << "the whole combined image (or to both). Note that all" << std::endl;
	std::cout << "transformations must currently be the same for all inserted" << std::endl;
	std::cout << "blocks. However, this can be changed in a straightforward" << std::endl;
	std::cout << "manner by extending the parameter set to per-block parameters." << std::endl;
}

void CommandSpliceRandomBlock::printShortHelp() const {
	std::cout << "splices images with random blocks, creates ground truth automatically" << std::endl;
}

} } // end namespace vole, cmfd_gt
