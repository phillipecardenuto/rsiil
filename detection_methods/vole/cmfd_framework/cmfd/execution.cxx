/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "execution.h"

#include "pca.h"
#include "kpca.h"
#include "wavelet.h"
#include "fastsats.h"
#include "buildcluster.h"
#include "cmfd_util.h"
#include "preproc.h"
#include "blockhandling.h"
#include "featfactory.h"
#include "matching.h"
#include "verification.h"
#include "mark.h"
#include "cps.h"
#include "log.h"
#include "post_process/post_process_core.h"
// for brisk we need the modules brisk and agast
// furthermore we have to compile with ssse enabled (see brisk/CMakeList.txt)
//#include <brisk.h>

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/nonfree/features2d.hpp>
#if (CV_MAJOR_VERSION >= 2 && CV_MINOR_VERSION > 3)
	#include <opencv2/nonfree/nonfree.hpp>
#endif
//#include <opencv2/features2d/features2d.hpp>

#include <boost/thread.hpp>
//#include <boost/bind.hpp>

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <time.h>

//using namespace boost::program_options;

// forwards for static elements
int Execution::verbosity;
std::string Execution::outputdir;
std::string Execution::image_basename;
std::string Execution::suffix;

Execution::Execution(CmfdConfig &_config)
	: config(_config)
{
    imgBlocks = NULL;
    feature = NULL;
    match = NULL;
    veri = NULL;
    cps = NULL;
    verbosity = config.verbosity;
    outputdir = config.outputdir;
    image_basename = cmfd::getBasename(config.inputfile);
	suffix = config.suffix;
	time = (double)cv::getTickCount();
	old_time = time;

    if (config.f.useCps)
        config.f.featvec = NOTHING;

    // how many threads do we want?
    if (config.numThreads == -1){
        threadnum = boost::thread::hardware_concurrency();
    }
	else {
		threadnum = std::min(config.numThreads, static_cast<int>(boost::thread::hardware_concurrency()));
    }
	if (threadnum == 0){
        threadnum = 1;
    }

    // create keypoint detector and feature extractor
    if ( ! config.f.detectorType.empty() || ! config.f.descriptorType.empty()){
        //		if (config.f.detectorType == "SURF"){
        //			cv::FeatureDetector *fd = new cv::DynamicAdaptedFeatureDetector
        //									(new cv::SurfAdjuster(), 1000, 2000, 50);
        //			feature_detector = cv::Ptr<cv::FeatureDetector>(fd);
        //		}
        if (config.f.detectorType == "SIFT"){
			vout(2) << "MY DETECTORTYPE is SIFT\n";
// think the parameters of SIFT changed or am I wrong?
#if (CV_MAJOR_VERSION <= 2 && CV_MINOR_VERSION <= 3)
            cv::FeatureDetector *fd = new cv::SiftFeatureDetector(cv::SIFT::DetectorParams::GET_DEFAULT_THRESHOLD(),
                                                                  cv::SIFT::DetectorParams::GET_DEFAULT_EDGE_THRESHOLD(),
                                                                  7);
			// Note: this is unfortunately not implemented yet for SIFT in opencv
			//CV_Assert(cv::SIFT::DescriptorParams::DEFAULT_IS_NORMALIZE == true);
			feature_detector = cv::Ptr<cv::FeatureDetector>(fd);

#else
			// lets create it over the standard way the modification of the 7 octave-layers
			// has also only been used by Pan (not Amerini)
			// default: _nOctaveLayers=3,
			vout(2) << "opencv version is: " << CV_MAJOR_VERSION << " " << CV_MINOR_VERSION << " " << CV_SUBMINOR_VERSION << std::endl;
			cv::FeatureDetector *fd = new cv::SiftFeatureDetector();
			feature_detector = cv::Ptr<cv::FeatureDetector>(fd);
			//std::cout << "B\n";
			//feature_detector = cv::FeatureDetector::create( config.f.detectorType );
			//std::cout << "C\n";
#endif
        }
		else if (config.f.detectorType == "ORB"){
			// this should be img-size dependent
			int nr_features = 2000;
			feature_detector = cv::Ptr<cv::FeatureDetector>(
						new cv::OrbFeatureDetector(nr_features) );
		}
		/*
		else if (config.f.detectorType == "BRISK"){
			feature_detector = cv::Ptr<cv::FeatureDetector>(
						new cv::BriskFeatureDetector(60,4) );
		}
		*/
        else {
			// Note: When OpenCV is built statically, dynamically created classes
			// (via Algorithm::create) can fail because linker excludes the "unused"
			//       object files. To avoid this problem, create classes explicitly
            feature_detector = cv::FeatureDetector::create( config.f.detectorType );
        }
		if (feature_detector.empty()){
			throw std::runtime_error("Execution::Execution:feature detector type: " + config.f.detectorType + " is not supported");
		}
		/*
		if (config.f.descriptorType == "BRISK"){
			descriptor_extractor = cv::Ptr<cv::DescriptorExtractor>(
						new cv::BriskDescriptorExtractor());
		}
		else {
		*/
			descriptor_extractor = cv::DescriptorExtractor::create( config.f.descriptorType );
		//}
		if (descriptor_extractor.empty()){
			throw std::runtime_error("Execution::Execution:feature descriptor type: " + config.f.descriptorType + " is not supported");
		}
    }

    if (config.writeCorresFile && config.loadCorresFile.size() > 0)
		throw std::runtime_error("Execution::Execution: you can't load and write at the same time!");
}

Execution::~Execution(){
    delete veri;
    delete cps;
}

void Execution::loadImage(){
	vout(1) << "Load Image " << config.inputfile << std::endl;

	if ( ! std::ifstream(config.inputfile.c_str()) ){
		throw std::runtime_error(std::string("Execution::loadImage: can't read: "
											 + config.inputfile));
	}
	if ( (imgorig = cv::imread(config.inputfile, 1)).empty() ){
		throw std::runtime_error("Exe::loadImage: error in loading");
	}
	img = imgorig.clone(); // working copy

	vout(1) << "Loaded image " << config.inputfile << " with size: "
			<< img.cols << "x" << img.rows << std::endl
			<< "-->image basename: " << Execution::image_basename << std::endl;
}

void Execution::preprocess(){
    if (img.empty()) {
        throw std::runtime_error("Exe::generateBlocks: no img loaded");
    }

    // channels
    Preproc p(img, config.p);

    // Gaussian pyramid
    if (config.p.pyramidLvl > 0) {
        img = p.computeGaussianPyramid();
        if (config.graphical){
            cv::imshow("pyra", img);
        }
    }

    // DWT
    if (config.p.waveletLvl > 0) {
        if (config.graphical){
            cv::Mat dstimg;
            Wavelet::decomposeHaar(img, dstimg, config.p.waveletLvl);
            cv::imshow("wavelet-decomposition", dstimg);
        }
        img = p.computeHaarWavelet();
    }

	if (verbosity >= 1)
		std::cout << "Preprocessing done: img-dimension: " << img.cols << "x" << img.rows << std::endl;

    if (config.graphical){
        cv::imshow("after preprocessing", img);
        cv::waitKey(0);
    }
    dimy = img.rows;
    dimx = img.cols;
	//	Feature::logPolarTransform(img);

	markMatrix = cv::Mat_<uchar>(img.rows, img.cols);
}

void Execution::generateBlocks()
{
	if (verbosity >= 1)
		std::cout << "->Generate Blocks\n";

    if (img.empty())
        throw std::runtime_error("Exe::generateBlocks: no img loaded");

	// detect keypoints
	if ( !feature_detector.empty() ){
		if (verbosity >= 1)
			std::cout << "compute keypoints\n";
		feature_detector->detect( img, keypoints );
		if (verbosity >= 1)
			std::cout << " detectected " << keypoints.size() << " keypoints" << std::endl;
		config.b.blockSize = 1;
	}
    // no img Blocks will be generated atm
	else if ( config.f.useCps ){
        cps = new Cps(img, imgorig, config.f, config.p.waveletLvl);
        cps->findRegion();
    }
    else {
        imgBlocks = new BlockHandling(img, config.b);
		imgBlocks->computeBlocks();
        if (verbosity >=1 ) {
            std::cout << " computed " << imgBlocks->size() << " blocks\n";
            int expected = ( (dimx-config.b.blockSize+config.b.step)
                             * (dimy-config.b.blockSize+config.b.step) )
                    / config.b.step;
            std::cout << " should be: " << expected << std::endl;
        }
// DEBUG
//		cv::Mat bl = (*imgBlocks)[0];
//		imgBlocks->printBlock(bl);
    }
    // remove blocks with too low entropy
    if (config.checkEntropy && imgBlocks != NULL){
        imgBlocks->checkEntropy(threadnum, config.m.entropyTh, config.f.circle);
		if (verbosity >= 1)
			std::cout << " new # blocks: " << imgBlocks->size() << std::endl;
    }
	if (verbosity >= 1)
		std::cout << "<-Generate Blocks\n";
}

void Execution :: computeFeatures()
{
	if (verbosity >= 1)
		std::cout << "Compute now features\n";
	if (imgBlocks != NULL && imgBlocks->empty() && cps == NULL && keypoints.empty()){
		throw std::runtime_error("Exe::computeFeatures: there are no blocks - call generateBlock() before");
    }
	// TODO/FIXME this is a small ugly hack to get the waveletLvl / pyramidLvl also in ft.cfg // fix that perhaps
    // also allows atm either wavelet decomp or gaussian pyramid decompose
	if (config.p.waveletLvl > 0){
        config.f.decomposeLvl = config.p.waveletLvl;
	}
	else if (config.p.pyramidLvl > 0){
        config.f.decomposeLvl = config.p.waveletLvl;
	}

    // keypoints
    if ( ! keypoints.empty() ){
        if (descriptor_extractor.empty()){
            throw std::runtime_error("Execution::computeFeatures(): have keypoints but no descriptor");
        }
		if (verbosity >= 1)
			std::cout << "feature extraction with kepoint based method\n";

		// compute descriptor at keypoints
		descriptor_extractor->compute(img, keypoints, feature_matrix);
		if (config.f.normalize){
			Feature::normalize(feature_matrix);
        }

        if (verbosity >= 2){
            std::cout << "Feature: " << feature_matrix.rows
                      << " x " << feature_matrix.cols << std::endl;
            std::cout << "print first feature\n";
            Feature::printFeature(feature_matrix, 0);
            std::cout << "print last feature: " << feature_matrix.rows - 1 << "\n";
            Feature::printFeature(feature_matrix, feature_matrix.rows-1);
        }
    }
	// block-based
    else {
        feature = FeatureFactory::getInstance(config.f.featvec, config.f, *imgBlocks, threadnum);

        if (feature == NULL){
            throw std::runtime_error("Exe::computeFeatures: feature == NULL -> call initFeature first");
        }

		// block-based
        feature->computeFeatureVec();
        if (config.f.normalize){
			feature->normalize();
        }
        if (verbosity >= 2){
            std::cout << "Feature: " << feature->getFeature().rows
                      << " x " << feature->getFeature().cols << std::endl;
            std::cout << "print first feature\n";
            feature->printFeature(0);
            std::cout << "print last feature: " << imgBlocks->size()-1 << "\n";
            feature->printFeature(imgBlocks->size()-1);
        }
    }

    // reduce the features
	if (config.f.usePCA){
		if (verbosity >= 1)
			std::cout << "\tCompute PCA\n";

		PCA p(*feature, config.f.eps);
		p.reduce(config.f.dim, config.f.minFeatureSize, config.f.maxFeatureSize);
		if (verbosity >= 2){
			std::cout << "Feature: " << feature->getFeature().rows
					  << " x " << feature->getFeature().cols << std::endl;
			std::cout << "print first feature\n";
			feature->printFeature(0);
			std::cout << "print last feature: " << imgBlocks->size()-1 << "\n";
			feature->printFeature(imgBlocks->size()-1);
		}
    }
}

void Execution::computeMatching(){
    if (config.f.useCps)
        return;   
    // TODO in future: change perhaps blocks completly to keypoints
    if (imgBlocks != NULL){
        keypoints = imgBlocks->getKeypoints();
    }

	// add standard threshold for euclidian ratio
	if ( config.m.crit.size() == 1
			&& pow(2.0, (double)(config.m.crit[0]-1)) == EUCLIDIAN_RATIO
			&& config.m.th.empty() )
	{
		config.m.th.push_back(0.6);
	}

	while(1){
		if ( feature_matrix.empty() ) { // block-based
			match = new Matching(keypoints, feature->getFeature(), img.size(), config.m);
		}
		else { // keypoint-based
			match = new Matching(keypoints, feature_matrix, img.size(), config.m);
		}

		if (verbosity >= 1)
			std::cout << "Match the blocks\n";

		match->computeMatching();
		corres_matrix = match->getCorresMatrix();
		int mc = getMatchCount();
		vout(2) << "Matchcount: " << mc << std::endl;

		if (!config.m.crit.empty() && !config.m.th.empty()){
			vout(2) << " crit-size: " << config.m.crit.size() << " th-size: "
					<< config.m.th.size() << " crit-nr: " << config.m.crit[0]
					<< " adjusted crit-nr: " << pow(2.0, (double)(config.m.crit[0]-1))
					<< " th: "  << config.m.th[0] << " ER-nr: " << EUCLIDIAN_RATIO << std::endl;
		}
		else {
			break;
		}
		if (mc < 50
				&& config.m.crit.size() == 1
				&& pow(2.0, (double)(config.m.crit[0]-1)) == EUCLIDIAN_RATIO
				&& config.m.th.size() == 1 && config.m.th[0] <= 0.8 )
		{
			config.m.th[0] += 0.05;
		}
		else if (mc > 1000
				 && config.m.crit.size() == 1
				 && pow(2.0, (double)(config.m.crit[0]-1)) == EUCLIDIAN_RATIO
				 && config.m.th.size() == 1 && config.m.th[0] >= 0.2 )
		{
			config.m.th[0] -= 0.05;
		}
		else {
			break;
		}
		// clean-up
		delete match;
	}
	delete match;

    if (verbosity >= 2){
		std::cout << getMatchCount() << " correspondences after matching\n";
    }
	if (verbosity >= 3){
		cv::Mat_<uchar> all_matches = cv::Mat_<uchar>::zeros( corres_matrix[0].rows,
									 corres_matrix[0].cols );
		for (size_t i = 0; i < corres_matrix.size(); i++){
			for (int y = 0; y < corres_matrix[i].rows; y++) {
				for (int x = 0; x < corres_matrix[i].cols; x++) {
					if (corres_matrix[i](y,x).x != -1){
						all_matches(y,x) = 255;
					}
				}
			}
		}
		Mark::dumpMatrix("_matches", all_matches);
	}

    // we don't need these classes any more
	//delete match;
    if (feature != NULL){
        delete feature;
        feature = NULL;
    }
	if (imgBlocks != NULL){
        delete imgBlocks;
        imgBlocks = NULL;
    }
    match = NULL;
}

void Execution::verifySimilarity(){
	if (config.f.useCps)
		return;
	if (veri != NULL) {
		delete veri;
	}
	if (corres_matrix.empty()) {
		throw std::runtime_error("Exe::initVerification: corres_matrix.empty()");
	}
	if (verbosity >= 1)
		std::cout << "Init Verification with blocksize: "
				  << (config.f.descriptorType.empty() ? config.b.blockSize : 0) << std::endl;

	veri = new Verification(config.f.descriptorType.empty() ? config.b.blockSize : 0,
							dimx,
							dimy,
							corres_matrix,
							config.p.pyramidLvl + config.p.waveletLvl);

	if (veri == NULL){
		throw std::runtime_error("Exe::verifySimilarity: sim = NULL -> call initSimilarity() before");
	}

	if (config.v.checkShift){
		if (verbosity >= 1)
			std::cout  << "Check now the shift vectors\n";
		if ( config.v.fastsats || config.v.cluster){
			throw std::runtime_error("Execution::verifySimilarity: only checkShift should be selected");
		}

		veri->verificateShift(config.v.minSameShift, config.v.numMax, config.v.shiftVariance);
	}

	if (config.v.fastsats){
		if (verbosity >= 1)
			std::cout 	<< "start FastSATS\n"
						<< "corres_matrix.size() = " << corres_matrix.size() << std::endl;
		Fastsats fast_sats(corres_matrix, config.b.blockSize, verbosity,
						   config.v.maxDist, config.b.step,
						   config.p.chan == "ALL",
						   config.v.withTree );
		markMatrix = fast_sats.sameAffineTransformationSelection(config.v.minSameShift,
																 config.v.recomputationTh);
		transformations = fast_sats.getTrafoMatrices();
//		if ( !config.v.estimateMethod.empty() )
		point_groups = fast_sats.getPointGroups();
	}

	if ( config.v.cluster ){
		if (verbosity >= 1)
			std::cout  << "start Clustering\n";
		float cutTh = config.v.cutTh;
		if (config.v.cutTh == 0.0){ // automatic adjustment
			if (getMatchCount() <= 100){
				cutTh = 75;
			}
			else {
				if (config.f.descriptorType == "SIFT")
					cutTh = 25;
				else
					cutTh = 50;
			}
		}
		cluster::BuildCluster cluster(corres_matrix, cutTh);
		if ( config.v.fastsats ){
			cluster.compute(point_groups);
		}
		else {
			cluster.compute();
		}
		point_groups = cluster.getSplitPointGroups();
		//			if (point_groups.size() < 20) {
		//				break;
		//			}
		//			cutTh += 20;
		//			cluster.setCutThreshold(cutTh);
		//		}

		// Ransac onto clusters (as Amerini did it)
		if ( !config.v.estimateMethod.empty() ){
			transformations = veri->ransac( config.v.estimateMethod,
											point_groups,
											config.v.ranIterations,
											config.v.ranReprojTh,
											config.v.normalizeCluster );
		}
		cv::Mat cluster_mark_matrix = cluster.getMarkedMatrix( config.b.blockSize / 2.0 );

		// output
		if ( config.mk.writeCluster ){
			std::stringstream ss;
			ss << Execution::outputdir << "/" << Execution::image_basename << "_cluster"
			   << Execution::suffix << ".png";
			vout(1) << "write " << ss.str() << std::endl;
			cv::imwrite(ss.str(), cluster_mark_matrix);
		}

		// build marking matrix
		std::vector<cv::Mat> planes;
		cv::split(cluster_mark_matrix, planes);
		// make everything to white
		cv::threshold(planes[0], markMatrix, 0.1, 255, cv::THRESH_BINARY);
	}

	// General Ransac onto the correspondence maps (as Pan did it)
	if ( !config.v.estimateMethod.empty() && !config.v.cluster ){
		// not supported atm
		if ( config.v.checkShift ){
			return;
		}
		vout(1)  << "start RANSAC\n";
		if ( config.v.fastsats ){
			// Ransac onto clusters (as Amerini did it)
			// but this time onto sats-clusters
			if ( !config.v.estimateMethod.empty() ){
				transformations = veri->ransac( config.v.estimateMethod,
												point_groups,
												config.v.ranIterations,
												config.v.ranReprojTh,
												config.v.normalizeCluster );
			}
		}
		else {
			// variante of Pan
			transformations.push_back( veri->ransac(config.v.estimateMethod,
													config.v.ranIterations,
													config.v.ranReprojTh,
													config.v.normalizeCluster ) );
		}
	}
}

// REFACTORING needed: this is not really a clean solution
// revise that function and also the Mark-class (mark.cxx/.h)!
void Execution::markIt(){
	if (verbosity >= 1)
		std::cout  << "Mark similar blocks\n";

    if ( ! config.loadCorresFile.empty() && img.empty() && config.mk.writeMatrix ){
		if (config.v.checkShift){
            Mark::writeMatrix(veri->getSimilarityMatrix());
		}
		else {
			Mark::writeMatrix(markMatrix);
		}
		return;
	}

	if (img.empty())
		throw std::runtime_error("Exe::markIt: no img loaded");

	config.post.verbosity = verbosity;

	Mark mrk(img, imgorig, config.mk, corres_matrix, config.b.blockSize);

	cv::Mat_<uchar> post(markMatrix.rows, markMatrix.cols, (uchar) 0);
	if ( config.mk.computeCorrelMap )
	{
		std::vector<cv::Mat_<uchar> > correl_map = mrk.computeCorrelation(transformations,
																		  config.mk.writeCorrelation);
		// apply post processing to all individual correlation maps
		for (size_t i = 0; i < correl_map.size(); i++){
			PostProcess pp(config.post, correl_map[i]);
			if ( config.post.rmSmallAreas ){
				pp.applyAreaThreshold(markMatrix);
			}
			if ( ! config.post.morphOp.empty() ){
				pp.applyMorphFilter();
			}
			std::vector<std::vector<cv::Point> > contours;
			std::vector<cv::Vec4i> hierarchy;
			if ( config.post.fillHolesByContour ){
				pp.fillHolesByContour(contours, hierarchy);
			}
			cv::Mat tmp = pp.getMatrix();
			post += tmp;

			// dump post-processed correlation map
			if (verbosity >= 3){
				std::stringstream ss;
				ss << "_pp_" << i;
				mrk.dumpMatrix(ss.str(), tmp);
			}
		}

		// remove areas which don't have a correspondence
		PostProcess pp(config.post, post);
		pp.applyAreaThreshold(cv::Mat(), &corres_matrix);
		post = pp.getMatrix();
		// mark in the normal image the contour
		if ( config.mk.markContour ){
			mrk.markContour(post, NULL, NULL);
		}
	}
	else {
		PostProcess pp(config.post, markMatrix);
		// this doesn't make much sense for sparse point-clouds
		//		if ( config.post.rmSmallAreas ){
		//			pp.applyAreaThreshold(cv::Mat());
		//		}
		if ( ! config.post.morphOp.empty() ){
			pp.applyMorphFilter();
		}
		std::vector<std::vector<cv::Point> > contours;
		std::vector<cv::Vec4i> hierarchy;
		if ( config.post.fillHolesByContour ){
			pp.fillHolesByContour(contours, hierarchy);
		}

		post = pp.getMatrix();
		// mark in the normal image the contour
		// here we can reuse the already found contours
		if ( config.mk.markContour ){
			mrk.markContour(imgorig, &contours, &hierarchy);
		}
	}
	//mrk.markHull(point_groups);

	if (config.mk.writePost){
		//mrk.dumpMatrix("_post", post);
		// Adding line to label the detected regions with the same id
		//and write a labeled image
		mrk.labelRegions(point_groups,post);
	}
	if (veri == NULL && cps == NULL){
		throw std::runtime_error("Exe::markIt: veri = NULL -> call initSimilarity() before");
	}
	if (config.mk.writeChrom){
		mrk.writeChrom();
	}
	if (config.mk.writePairs){
		mrk.writePairs( config.groundTruthFile, false /*config.f.detectorType.empty()*/ );
	}
	if (config.f.useCps){
		mrk.markIt(cps->getSimilarityMatrix());
		return;
	}
	if (config.mk.markRegions){
		if (config.v.checkShift){
			mrk.markIt(veri->getShift());
		}
		else{
			mrk.markRegions(point_groups, config.mk.markShift);
		}
	}
	if ( config.v.fastsats
			|| ! config.v.estimateMethod.empty()
			|| config.v.cluster){
		mrk.markIt(markMatrix);
	}
	else{
		mrk.markIt(veri->getSimilarityMatrix());
	}

}

// conversion function
std::vector<cv::Mat> Execution :: corresToImg() const
{
    if (corres_matrix.empty())
        throw std::runtime_error("Execution::corresToImg: corres_matrix is empty!");

	if (verbosity >= 1)
		std::cout  << "compute correspondence matrix to img-matrix\n";

    // convert correspondence matrix in a format that can be saved as image
    std::vector<cv::Mat>  corres_in_img;
    for (size_t c = 0; c < corres_matrix.size(); c++){
        cv::Mat_<cv::Vec3w>	M(dimy, dimx, cv::Vec3w(0,0,0));
        for (int y = 0; y < dimy; y++){
            for (int x = 0; x < dimx; x++){
                if (corres_matrix[c](y,x).x < 0)
                    continue;
                cv::Point pb = corres_matrix[c](y,x);
                // one dimension won't be greater than 2^16
                // thus convert the ints to 16bit unsigned short
                M(y,x)[0] = static_cast<unsigned short>(pb.x);
                M(y,x)[1] = static_cast<unsigned short>(pb.y);
            }
        }
        // add to the vector:
        corres_in_img.push_back(M);
    }
    return corres_in_img;
}

// conversion function
void Execution :: imgToCorres(const std::vector<cv::Mat> & corres_in_img)
{
    if (corres_in_img.empty()){
        throw std::runtime_error("Execution::imgToCorres: corres_in_img is empty!");
    }
	if (verbosity >= 1)
		std::cout  << "compute correspondence matrix from img\n";

    for (size_t c = 0; c < corres_in_img.size(); c++){
        cv::Mat_<cv::Point>	M(dimy, dimx, cv::Point(-1,-1));
        for (int y = 0; y < dimy; y++){
            for (int x = 0; x < dimx; x++){
                if (corres_in_img[c].at<cv::Vec3w>(y,x)[0] == 0
                        && corres_in_img[c].at<cv::Vec3w>(y,x)[1] == 0)
                    continue;
                // cast the 16bit unsigned short to int
                M(y,x).x = static_cast<int>(corres_in_img[c].at<cv::Vec3w>(y,x)[0]);
				M(y,x).y = static_cast<int>(corres_in_img[c].at<cv::Vec3w>(y,x)[1]);
            }
        }
        // add to the vector:
        corres_matrix.push_back(M);
    }
}

int Execution :: getBlockCount() const
{
	if (!keypoints.empty()){
		return keypoints.size();
	}
	if (imgBlocks != NULL){
		return imgBlocks->size();
	}
	return 0;
}
int Execution :: getMatchCount() const
{
	int cnt = 0;
	for (size_t i = 0; i < corres_matrix.size(); i++){
		for (int y = 0; y < corres_matrix[i].rows; y++) {
			for (int x = 0; x < corres_matrix[i].cols; x++) {
				if (corres_matrix[i](y,x).x != -1){
					cnt++;
				}
			}
		}
	}
	return cnt;
}
