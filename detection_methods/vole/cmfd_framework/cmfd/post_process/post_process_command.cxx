#include "post_process_command.h"
#include "cmfd_util.h"
#include "post_process_core.h"
// OpenCV
#include <opencv2/highgui/highgui.hpp>

#include <sstream>
using namespace boost::program_options;

PostProcessCommand :: PostProcessCommand(void)
 : Command(
		"pp", //post processing
		config,
		"Vincent Christlein",
		"vincent.christlein@stud.informatik.uni-erlangen.de")
{
}

int PostProcessCommand :: execute(void) 
{
	if (config.verbosity >= 1){
		printConfig();
		std::cout << "---EXECUTE---\n";
	}
	cv::Mat *img = execute_headless();;
	if (config.graphical) {
		// show what we painted using opencv
		cv::imshow("Forgeries", *img);
		cv::waitKey();
	}
	std::stringstream ss;
	ss << config.outputdir << "/" << cmfd::getBasename(config.inputfile) << "_pp.png";
	if (config.verbosity >= 1){
		//std::cout << "save image to " << config.outputdir << config.f.featvec << "_out.png" << std::endl;
		std::cout << "save image to " << ss.str() << std::endl;
	}
	cv::imwrite(ss.str(), *img);

	delete i;
    return 0;
}

cv::Mat* PostProcessCommand :: execute_headless(void) 
{
	// check basename
	std::string image_basename = cmfd::getBasename(config.inputfile);
	if (image_basename.length() < 1) {
		throw std::runtime_error("ERROR: input image basename is empty?!");
	} 
	if (config.verbosity > 1) {
		std::cout << "image basename: " << image_basename << std::endl;
		std::cout << "load now image\n";
	}

	// load images
	cv::Mat img;
	if ( (img = cv::imread(config.inputfile)).empty() )
		throw std::runtime_error("PostProcessCommand::execute_headless: error in loading image");
	if ( img.channels() != 3 && img.channels() != 1 )
		throw std::runtime_error("PostProcessCommand::execute_headless: error in loading image - unsupported # of channels");
	if ( config.graphical ) {
		cv::imshow("original", img); 
	}
	cv::Mat matrix;
	if ( !config.matrix.empty() ){
		if ( (matrix = cv::imread(config.matrix)).empty() )
			throw std::runtime_error("PostProcessCommand::execute_headless: error in loading matrix-image");
	}
	cv::Mat orig;
	if ( !config.orig.empty() ){
		if ( (orig = cv::imread(config.orig)).empty() )
			throw std::runtime_error("PostProcessCommand::execute_headless: error in loading orig-image");
	}

	// the actual preprocessing is handled by the core-class
	PostProcess pp(config, img);
	if ( config.rmSmallAreas ){
		pp.applyAreaThreshold(matrix);
	}
	if ( !config.morphOp.empty() ){
		pp.applyMorphFilter();
	}
	if ( config.fillHolesByContour ){
		std::vector<std::vector<cv::Point> > contour;
		std::vector<cv::Vec4i> hierarchy;
		pp.fillHolesByContour(contour, hierarchy);
	}

	i = new cv::Mat(img.rows, img.cols, img.type()); // important, cause execution class exists only local
	*i = img.clone();
	return i;
}


void PostProcessCommand :: printShortHelp(void) const 
{
	std::cout << "post processing of detection maps generated by forgery detection methods, e.g. cmfd" << std::endl;
}

void PostProcessCommand :: printHelp(void)  const 
{
	printShortHelp();
	// TODO
}

void PostProcessCommand :: printConfig(void) 
{
//	config.printConfig();
	std::cout << "---PARAMETERS---\n"
			  << "inputImage: " << config.inputfile << std::endl
			  << "matrix: " << config.matrix << std::endl
			  << "orig: " << config.orig << std::endl
			  << "outputpath: " << config.outputdir << std::endl;
	std::cout << "\tareaThreshold: " << config.areaThreshold << std::endl;
	std::cout << "\tmorphOp: ";
	for (size_t i = 0; i < config.morphOp.size() ; i++)
		std::cout << config.morphOp[i] << " ";
	std::cout << std::endl;
	std::cout << "\tmorphIter: ";
	for (size_t i = 0; i < config.morphIter.size() ; i++)
		std::cout << config.morphIter[i] << " ";
	std::cout << std::endl;
	// TODO: rest
}