/*
 * Original work: Copyright (c) 2017, Thibaud Ehret <ehret.thibaud@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file main.cpp
 * @brief Main file
 *
 * @author Thibaud Ehret <ehret.thibaud@gmail.com>
 **/

#include "Utilities/LibImages.h"
#include "Utilities/cmd_option.h"
#include "Utilities/FeatManager/zMManager.h"
#include "Utilities/FeatManager/siftManager.h"
#include "Utilities/PatchMatch/patchmatch.h"
#include "Utilities/filters.h"
#include "Utilities/tools.h"



int main(int argc, char **argv) {

	// Read parameters
	clo_usage("Copy-move forgery detection based on Zernike moments or dense SIFT and PatchMatch");
	clo_help(" NOTE: Input (<) and output (>) sequences are specified by their paths in printf format.\n");

	//! Paths to input/output sequences
	using std::string;
	string input_path = clo_option("-i"    , "" , "< input image");
	string output_path = clo_option("-o"    , "" , "> output image");

	bool flip = (bool) clo_option("-flip"    , false, "< computes the matching on the flip image");

	// Parameters for the features
	int mtd = clo_option("-mtd"    , 0, "< method used to compute the features (0: Zernike moments, 1: sift), the default value is 0");
	int sp = clo_option("-sp"    , 8, "< radius of the zernike moments/bin size/patch size");
	int zm = clo_option("-zm"    , 1, "< if Zernike moments are used, choose which type of Zernike moments to use (0: basic Zernike moments, 1: resampled Zernike moments), the default value is 1");

	// Parameters for patchmatch
	int N = clo_option("-iter"    , 8, "< Nb of iterations");
	int th1 = clo_option("-th1"    , 64, "< minimum distance threshold for patchmatch");
	
	// Parameters post-processing
	int th_d = clo_option("-thd"    , 2500, "< minimum squared distance between clones");
	int th_e = clo_option("-the"    , 300, "< threshold on the error");
	int th_s = clo_option("-ths"    , 1200, "< minimum size for a clone");
	int radius_m = clo_option("-rdm"    , 4, "< radius of the median filter");
	int radius_e = clo_option("-rde"    , 6 , "< radius of the error filter");
	int radius_d = clo_option("-rdd"    , radius_e+radius_m, "< radius of the dilatation filter");

	if(input_path == "")
	{
		fprintf(stderr, "Options '-h', '-help' and '--help' list the different options available.\n",
				argv[0], argv[0]);
		return EXIT_FAILURE;
	}

	std::vector<float> image;
	ImageSize imSize;

	// Here we load the image
	loadImage(input_path.c_str(), image, imSize);
	transformColorSpace(image, imSize, true);

	// transform image to grayscale
	std::vector<float> imagegs(imSize.wh);

	for(int i = 0; i < imSize.wh; ++i)
		imagegs[i] = image[i];

	imSize.nChannels = 1;
	imSize.whc = imSize.wh;

	// Create the Zernike moments, or other SIFT descriptors
	FeatManager* fm;
	if(mtd == 0)
		fm = new ZMManager(imagegs, imSize, sp, 5, 26, 32, zm, flip);
	else if(mtd == 1)
		fm = new SiftManager(imagegs, imSize, sp, flip);

	ImageSize descrSize;
	descrSize = fm->sz();

	ImageSize visualSize;
	visualSize.width     = descrSize.width;
	visualSize.height    = descrSize.height;
	visualSize.nChannels = 1;
	visualSize.wh        = visualSize.width * visualSize.height;
	visualSize.whc       = visualSize.wh;

	ImageSize clvSize = visualSize;
	clvSize.nChannels = 3;
	clvSize.whc *= 3;
	std::vector<float> clvisual(clvSize.whc);

	/// Compute the matching using Patchmatch
	Patchmatch pm(*fm, th1);
	pm.initializeRandom();
	pm.propagateNtimes(N);

	// extract the displacement maps from the matching result of PatchMatch
	std::vector<int> dispX;
	std::vector<int> dispY;
	pm.extract(dispX,dispY);

	/// Median filtering of the displacement maps
	medianFilter(dispX, visualSize, radius_m, true);
	medianFilter(dispY, visualSize, radius_m, false); 
	
	// Save the filtered displacement map
	//view_displacement(clvisual, clvSize, dispX, dispY);
	// rescale(clvisual, clvSize);
	// string temp_path = "filteredPmDisp.png";
	// saveImage(temp_path.c_str(), clvisual, clvSize, 0., 255.); 

	std::vector<bool> detectionMask(visualSize.wh);
	std::vector<float> visual(dispX.size());

	/// compute error detection map
	errorDetectionFilter(detectionMask, dispX, dispY, visualSize, radius_e, th_e, sp, visual); 
	std::vector<int> matchIds(detectionMask.size()); // Label of each detected object
	pm.get_detected_matches(matchIds,detectionMask);

	// Save the error as well as the initial detection mask
//	rescale(visual, visualSize);
	//temp_path = "errorMap.png";
	//saveImage(temp_path.c_str(), visual, visualSize, 0., 255.); 
	//visual.assign(detectionMask.begin(), detectionMask.end());
	//rescale(visual, visualSize);
	//temp_path = "detectionMask.png";
	//saveImage(temp_path.c_str(), visual, visualSize, 0., 255.); 

	/// Compute connected components + remove the small ones
	sizeFilter(detectionMask, visualSize, th_s);

	/// Remove matches that are too close
	minDispFilter(detectionMask, dispX, dispY, th_d);

	/// Symmetrize the matches
	symmetrizationFilter(detectionMask, dispX, dispY, visualSize);

	/// Dilate the rest 
	dilationFilter(detectionMask, visualSize, radius_d);

	// Save the final detection mask
	visual.assign(detectionMask.begin(), detectionMask.end());
	rescale(visual, visualSize);
	//string temp_path;
	//if (output_path.compare("") == 0)
	//	temp_path = "filteredMask.png";
	//else
	//	temp_path = output_path;
	//saveImage(temp_path.c_str(), visual, visualSize, 0., 255.); 

	// Compute the final decision
	//bool forgery = false;
	//for(int i = 0; i < detectionMask.size(); ++i)
	//	if(detectionMask[i])
	//	{
	//		forgery = true;
	//		break;
	//	}

	//Save result image in a binary file in vector shape
	string temp_path;
	if (output_path.compare("") == 0){
		if (mtd)
		temp_path = "sift_result.bin";
		else
		temp_path = "zernike_result.bin";
	}
	else
	{
		temp_path = output_path;
	}
	

	labelDetectionMask(temp_path.c_str(),visual, matchIds, imSize, visualSize, radius_d);
	//function to save the objects detected ids
	

	/// Save the result in a separate txt file
	//ofstream file;
	//file.open("result.txt", ios::out);
	//if(forgery)
		//file << "This image IS a forgery" << endl;
	//else
		//file << "This image IS NOT a forgery" << endl;
	//file.close();
	delete fm;

	return 0;
}
