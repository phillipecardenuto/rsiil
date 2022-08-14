/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "copymoveframework.h"
#include "execution.h"
#include "log.h"
// OpenCV
#include "cv.h"
#include "highgui.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <sstream>
#include <fstream>

using namespace boost::program_options;

// in the constructor we set up all parameters the user may configure
CopyMoveFramework::CopyMoveFramework(void)
 : Command(
		"cmfd", //copy move forgery detection
		config,
		"Vincent Christlein",
		"vincent.christlein@stud.informatik.uni-erlangen.de")
{
}


int CopyMoveFramework::execute(void) {
	if (config.verbosity >= 1){
		printConfig();
		std::cout << "---EXECUTE---\n";
	}
	cv::Mat *img = execute_headless();
	if (img != NULL && config.graphical) {
        // show what we painted using opencv
		cv::imshow("Forgeries", *img);
		cv::waitKey();
	}
	if (img != NULL && config.mk.writeOutput && config.loadCorresFile.empty()){
		std::stringstream ss;
		ss << config.outputdir << "/" << Execution::image_basename << "_output" << Execution::suffix << ".png";
		if (Execution::verbosity >= 1)
			std::cout << "save image to " << ss.str() << std::endl;
		cv::imwrite(ss.str(), *img);
	}
	if (img != NULL)
		delete img;
    return 0;
}

cv::Mat* CopyMoveFramework :: execute_headless(void)
{
	Execution ex(config);
	cv::Mat img;	
	std::vector<cv::Mat> corres_as_img;

	Log log;

	double start_time = (double)cv::getTickCount();
	old_time = start_time;

	// no corres-file specified to load
	if (config.loadCorresFile.empty()){
		if (Execution::image_basename.length() < 1) {
			throw std::invalid_argument("ERROR: image_basename is empty?!");
		}
		vout(1) << "image basename: " << Execution::image_basename << std::endl;

		if (config.log_times && config.log_file.empty()){
			config.log_file = std::string(config.outputdir + '/' + Execution::image_basename
										  /*+ config.suffix*/ + "_times" +  ".log");
		}
		log.set(config.log_file , config.log_times ? 1 : 0);

		ex.loadImage();

		img = ex.getImage();
		if (config.graphical) {
			cv::imshow("original", img); 
		}
		int iterations = 1;
		std::vector<std::string> chans;
		// FIXME: this currently only works if we use keypoint-based features
		if (config.p.chan == "ALL" and config.f.descriptorType != "" ){
			iterations = 3;
			chans.push_back("RED");
			chans.push_back("GREEN");
			chans.push_back("BLUE");
		}
		else {
			chans.push_back(config.p.chan);
		}

		for (int i = 0; i < iterations; i++){
			config.p.chan = chans[i];
			ex.preprocess();
			ex.generateBlocks();
			log() << "BC: " << ex.getBlockCount() << "\n"
				  << "PT: " << relativeTime() << "\n";

			ex.computeFeatures();
			double t = relativeTime();
			if (Execution::verbosity >= 2)
				std::cout << "Feature computation time needed: " << t << std::endl;
			log() << "FT: " << t << "\n";

			ex.computeMatching();
			log() << "MT: " << relativeTime() << "\n";

			std::vector<cv::Mat> corres = ex.corresToImg();
			// insert vector at end of overall vector
			corres_as_img.insert(corres_as_img.end(), corres.begin(), corres.end());
		}
	}
	else {
		// 1. load serialized infos from file
		std::ifstream ifs(config.loadCorresFile.c_str());
		if (!ifs) {
			throw std::invalid_argument("CopyMoveFramework::execute_headless: "
										"correspondence-file could not be loaded");
		}
		boost::archive::text_iarchive ia(ifs);
		ia >> ex;

		if (config.log_times && config.log_file.empty()){
			config.log_file = std::string(config.outputdir + '/' + Execution::image_basename
										  /*+ config.suffix*/ + "_times" +  ".log");
		}
		log.set(config.log_file , config.log_times ? 1 : 0);

		// 2. now load all correspondence maps saved as imgs 
		int i = 0;
		while(true){
			std::stringstream ss;
			if (Execution::image_basename.empty()){
				std::cerr << "WARNING: image_basename is empty, try to construct from ser-file"
						  << std::endl;
			}
			else if (Execution::verbosity >= 2)
				std::cout << "image basename: " << Execution::image_basename << std::endl;

			// corresfile may be somewhere else than the output-directory is
			size_t found = config.loadCorresFile.find_last_of('/');
			if ( found != std::string::npos ){
				ss << config.loadCorresFile.substr(0, found+1);
			}
			ss << Execution::image_basename /*<< Execution::suffix */<< "_corr_" << i << ".png";

			cv::Mat corres_img = cv::imread(ss.str(), -1);
			if ( corres_img.empty() ){
				if ( i == 0 ){
					throw std::runtime_error(std::string("CopyMoveFramework::execute_headless: "
														 + ss.str() + " couldn't be read"));
				} else {
					break;
				}
			}
			if ( corres_img.type() != CV_16UC(3) ) {
				throw std::runtime_error(std::string("CopyMoveFramework::execute_headless: "
													 + ss.str() + " couldn't be read - wrong format"));
			}
			if (Execution::verbosity >= 1)
				std::cout << "loaded " << ss.str() << std::endl
						  << "image basename: " << Execution::image_basename << std::endl;
			corres_as_img.push_back(corres_img);
			i++;
		}
		ex.imgToCorres(corres_as_img);
		// load img if one is given (need one for correlation map to warp the image)
		if ( ! config.inputfile.empty() ) {
			ex.loadImage();
			ex.preprocess();
		}
	}

	ex.verifySimilarity();
	// normally we only want to log until we wrote the corresfile
	if (!config.writeCorresFile)
		log() << "VT: " << relativeTime() << "\n";

	ex.markIt();
	if (!config.writeCorresFile)
		// log correlation + marking time
		log() << "CT: " << relativeTime() << "\n";

	// write serialized infos to a serialization file and the correspondence map(s) to png
	if (config.writeCorresFile)
	{
		if (Execution::verbosity >= 1)
			std::cout << "write corres file\n";

		// 1. save serialization
		//ex.corresMatToVec();
		std::stringstream ss;
		ss << Execution::outputdir << "/" << Execution::image_basename << "_info.ser";
		std::ofstream ofs(ss.str().c_str());
		if (!ofs){
			throw std::runtime_error(std::string("CopyMoveFramework::execute_headless: " + ss.str() + " - file couldn't be saved"));
		}
		boost::archive::text_oarchive oa(ofs);
		oa << ex;

		// 2. now save correspondence map(s) as png
		std::vector<int> params;
		params.push_back(CV_IMWRITE_PNG_COMPRESSION);
		params.push_back(9); // highest compression
		for (size_t i = 0; i < corres_as_img.size(); i++){
			std::stringstream ss2;
			ss2 << Execution::outputdir << "/" << Execution::image_basename << "_corr_" << i << ".png";
			if (corres_as_img[i].type() != CV_16UC(3)){
				std::cerr << "WARNING: You are saving the correspondence map not as 16bit unsgigned short img!\n";
			}
			if ( ! cv::imwrite(ss2.str(), corres_as_img[i], params) ){
				throw std::runtime_error(std::string("CopyMoveFramework::execute_headless: " + ss2.str() + " couldn't be saved"));
			}
		}
	}
	if (!config.writeCorresFile){
		// overall-time & block-count
		log() << "OT: " << (cv::getTickCount() - start_time) / cv::getTickFrequency() << "\n"
			  << "MC: " << ex.getMatchCount() << "\n";
	}
	// annoying stuff
	if (config.loadCorresFile.empty()){
		i_img = new cv::Mat(img.rows, img.cols, img.type()); // important, cause execution class exists only local
		*i_img = img.clone();
		return i_img;
	}
	return NULL;
}


void CopyMoveFramework::printShortHelp(void) const {
	std::cout << "Copy-Move Forgery Detection Framework" << std::endl;
}


void CopyMoveFramework::printHelp(void)  const {
	std::cout << "The Copy Move Framework tries to detect forgeries by ";
	std::cout << "checking if parts of the image are copied to another place in the picture" << std::endl;
}

void CopyMoveFramework::printConfig(void) {
	std::cout << "---PARAMETERS---\n"
			  << "log: ";
	if (config.log_times){
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << "\n\tlog_file: " << config.log_file << std::endl;
	std::cout << "verbosity: " << config.verbosity << std::endl
		<< "inputImage: " << config.inputfile << std::endl
		<< "outputpath: " << config.outputdir << std::endl
		<< "suffix: " << config.suffix << std::endl
		<< "groundTruthFile: " << config.groundTruthFile << std::endl
		<< "numThreads: " << config.numThreads << std::endl
		<< "loadCorresFile: " << config.loadCorresFile << std::endl;
	std::cout << "writeCorresFile: ";
	if (config.writeCorresFile){
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << "\ncheckEntropy: ";
	if (config.checkEntropy){
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << std::endl;
	std::cout << "PREPROCESSING\n";
	std::cout << "\tchan: " << config.p.chan << std::endl;
	std::cout << "\tpyramidLvl: " << config.p.pyramidLvl << std::endl;
	std::cout << "\twaveletLvl: " << config.p.waveletLvl << std::endl;

	std::cout << "BLOCKS\n"
		<< "\tblocksize: " << config.b.blockSize << std::endl
		<< "\tstep: " << config.b.step << std::endl
		<< "\troi: ";
	for (unsigned int i = 0; i < config.b.roi.size(); i++){
		std::cout << config.b.roi[i] << " ";
	}
	std::cout << std::endl
		<< "FEATURE\n"
		<< "\tdetectorType: " << config.f.detectorType << std::endl
		<< "\tdescriptorType: " << config.f.descriptorType << std::endl
		<< "\tnormalize: ";
	if (config.f.normalize){
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << std::endl
		<< "\tfeaturevec: " << config.f.featvec << std::endl
		<< "\tmaxFeatureSize: " << config.f.maxFeatureSize << std::endl
		<< "\tminFeatureSize: " << config.f.minFeatureSize << std::endl
		<< "\tmodified: ";
	if (config.f.modified){
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << std::endl
		<< "\tcircle: ";
	if (config.f.circle){
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << std::endl
		<< "\tqf: " << config.f.qf << std::endl
		<< "\tnumHu: " << config.f.numHu << std::endl
		<< "\tdecomposeLvl: " << config.f.decomposeLvl << std::endl
		<< "\tdegZernike: " << config.f.degZernike << std::endl
		<< "\tusePCA: ";
	if (config.f.usePCA){
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << std::endl;
	std::cout << "\tdim: " << config.f.dim << std::endl
		<< "\teps: " << config.f.eps << std::endl
		<< "\tsigma: " << config.f.sigma << std::endl
		<< "\tnumSamples: " << config.f.numSamples << std::endl
		<< "\tdstwidth: " << config.f.dstwidth << std::endl

		<< "\tQUANTIZATION\n"
		<< "\t\tqStart: " << config.f.qStart << std::endl
		<< "\t\tqEnd: " << config.f.qEnd << std::endl
		<< "\t\tqNum: " << config.f.qNum << std::endl

		<< "\t\tqRound: ";
		if (config.f.qRound){
			std::cout << "true";
		}
		else {
			std::cout << "false";
		}
		std::cout << std::endl;
		std::cout << "\t\tqRound2: ";
		if (config.f.qRound2){
			std::cout << "true";
		}
		else {
			std::cout << "false";
		}
		std::cout << std::endl;

		std::cout << "\t\tqAbs: ";
		if (config.f.qAbs){
			std::cout << "true";
		}
		else {
			std::cout << "false";
		}

	std::cout << "\n\tuseCps: ";
	if (config.f.useCps){
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << "\n\t\tpartsInX: " << config.f.partsInX << std::endl;
	std::cout << "\t\tpartsInY: " << config.f.partsInY << std::endl;
	std::cout << "\tcorrelTh: " << config.f.correlTh << std::endl;
	std::cout << "\tpixelDiff: " << config.f.pixelDiff << std::endl;	

	std::cout << "\nMATCHING\n"
		<< "\tminDistEuclidian: " << config.m.minDistEuclidian << std::endl
		<< "\tminDistCheby: " << config.m.minDistCheby << std::endl
		<< "\tcrit: ";
	for (size_t i = 0; i < config.m.crit.size(); i++){
		std::cout << config.m.crit[i] << " ";
	}
	std::cout << "\n\tth: ";
	for (size_t i = 0; i < config.m.th.size(); i++){
		std::cout << config.m.th[i] << " ";
	}
	
	std::cout << "\n\teuclidianTh: ";
	if (config.m.euclidianTh) {
		std::cout << "true";
	}
	else {
		std::cout << "false";
	}
	std::cout << "\n\tkdsort: ";
	if (config.m.kdsort) {
		std::cout << "true";
	}
	else {
		std::cout << "false";
	}
	std::cout << "\n\tnnBrute: ";
	if (config.m.nnBrute) {
		std::cout << "true";
	}
	else {
		std::cout << "false";
	}

	std::cout << std::endl
		<< "\tnumRows: " << config.m.numRows << std::endl;
	std::cout << "\tcorrelTh: " << config.m.correlTh << std::endl;
	std::cout << "\tentropyTh: " << config.m.entropyTh << std::endl;

	std::cout	<< "VERIFICATION\n"
				<< "\tcluster: ";
	if (config.v.cluster) {
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << "\n\t\tcutTh: " << config.v.cutTh;
	std::cout << "\n\t\tnormalizeCluster: ";
	if (config.v.normalizeCluster) {
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << "\n\testimateMethod: " << config.v.estimateMethod;
	std::cout << std::endl;
	std::cout << "\t\tranReprojTh: " << config.v.ranReprojTh
		<< "\n\t\tranIterations: " << config.v.ranIterations
		<< std::endl;
	std::cout << "\n\tfastsats: ";
	if (config.v.fastsats) {
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << "\n\t\twithTree: ";
	if (config.v.withTree) {
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << "\n\t\trecomputationTh: " << config.v.recomputationTh;
	std::cout << "\n\tcheckShift: ";
	if (config.v.checkShift) {
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << std::endl
		<< "\tminSameShift: " << config.v.minSameShift << std::endl;
	std::cout << "\t\tnumMax: " << config.v.numMax << std::endl;
	std::cout << "\t\tshiftVariance: " << config.v.shiftVariance << std::endl;
	std::cout << std::endl
		<< "\t\tmaxDist: " << config.v.maxDist << std::endl;

	std::cout << "\n\tcomputeCorrelMap: ";
	if (config.mk.computeCorrelMap) {
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << "\nMARKING\n"
	<< "\tmarkRegions: ";
	if (config.mk.markRegions) {
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << std::endl
			  << "\tmarkShift: ";
	if (config.mk.markShift) {
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << std::endl
			  << "\tmarkContour: ";
	if (config.mk.markContour) {
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << std::endl
			  << "\twriteChrom: ";
	if (config.mk.writeChrom) {
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << std::endl;
	std::cout << "\twritePairs: ";
	if (config.mk.writePairs) {
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << std::endl;
	std::cout << "\twriteMatrix: ";
	if (config.mk.writeMatrix) {
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << std::endl;
	std::cout << "\twriteOutput: ";
	if (config.mk.writeOutput) {
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << std::endl;
	std::cout << "\twriteCluster: ";
	if (config.mk.writeCluster) {
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << std::endl;
	std::cout << "\twriteCorrelation: ";
	if (config.mk.writeCorrelation) {
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << std::endl;
	std::cout << "\twritePost: ";
	if (config.mk.writePost) {
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << "\n\tscaleUp: ";
	if (config.mk.scaleUp) {
		std::cout << "true";
	} else {
		std::cout << "false";
	}
	std::cout << "\n\tuseOrig: ";
	if (config.mk.useOrig) {
		std::cout << "true";
	} else {
		std::cout << "false";
	}	

	std::cout << std::endl;
	std::cout << "\tminVal: " << config.mk.minVal << std::endl;
	std::cout << "\tbinaryTh: " << config.mk.binaryTh << std::endl;
	std::cout << config.post.getString();
}

