/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "mark.h"
#include "execution.h"
#include "matching.h"

#include <vector>
#include <queue>
#include <cstdio>
#include "highgui.h"


// TODO: this class should be revised!

Mark :: Mark(cv::Mat & _img,
			 cv::Mat & _origimg,
			 struct mark_cfg &_config,
			 const std::vector<cv::Mat_<cv::Point> > & _corres,
			 int _blockSize)
	: img(_img),
	  origimg(_origimg),
	  config(_config),
	  corres(_corres)
{
	block_size = _blockSize;
	params.push_back(CV_IMWRITE_PNG_COMPRESSION);
	params.push_back(9); // highest compression
}

void Mark :: writeMatrix(const cv::Mat& matrix) 
{
	if (matrix.empty() || !matrix.data){
		std::cerr << "WARNING: Mark::writeMatrix: matrix is empty -> return\n";
		return;
	}
	std::stringstream ss;
	ss << Execution::outputdir << "/"
	   << Execution::image_basename << "_matrix" << Execution::suffix << ".png";

	// some output
	if (Execution::verbosity >= 1) {
		int sum = 0;
		std::vector<cv::Mat> planes;
		cv::split(matrix, planes);
		for (int c = 0; c < matrix.channels(); c++){
			sum += cv::countNonZero(planes[c]);
		}
		std::cout << "mark " << sum << " pixels in "<< matrix.channels() << " plane(s) \n";
	}

	// just take first plane cause we marked only in that
	// and write out as gray image
	// TODO: change this behaviour perhaps
	/*
	if (chan != ALL && matrix.channels() > 1){
		std::vector<cv::Mat> matplanes;
		cv::split(matrix, matplanes);
		//TODO
		cv::imwrite(ss.str(), matrix);
	}
	else
	*/
	std::vector<int> parameters;
	parameters.push_back(CV_IMWRITE_PNG_COMPRESSION);
	parameters.push_back(9); // highest compression
	if	( ! cv::imwrite(ss.str(), matrix, parameters) ){
		std::cerr << "WARNING: couldn't write " << ss.str() << std::endl;
	}
	else if (Execution::verbosity >= 1){
		std::cout << "wrote " << ss.str() << std::endl;
	}
}

void Mark :: markIt(cv::Mat & matrix) 
{
	int mult = origimg.rows / img.rows; // get multiplier

	cv::Mat markImg;
	if (config.scaleUp) 
	{
		markImg = origimg;
		if (mult > 1){ // if decomposition just scale up our matrix...
			cv::Mat tmp;
			cv::resize(matrix, tmp, origimg.size(), 0, 0, cv::INTER_LINEAR);
			matrix = tmp;
		}
	}
	else {
		if (config.useOrig)
			markImg = origimg;
		else{
			markImg = img;

			if ( mult > 1 ) // decomposition -> normalize
			{
				std::vector<cv::Mat> planes;
				cv::split(markImg, planes);

				std::vector<cv::Mat > planesTo;
				planesTo.resize(markImg.channels());

				for (int i = 0; i < markImg.channels(); i++){
					//normalize
					double max;
					double min;
					cv::minMaxLoc(planes[i], &min, &max);
					planes[i] -= min; // move min to 0
					// TODO: use scale fct of opencv if this exists...
					if ( (max - min) > 0.0)
						planes[i] *= (255 / (max - min)); // scale  to 0-255
					
					// now convrt to uchar
					cv::Mat_<uchar> tmp = planes[i];
					planesTo[i] = tmp;
				}
				cv::merge(planesTo, markImg);
			}
		}
	}

	// apply a threshold minVal
	if (matrix.channels() == 1){
		if (Execution::verbosity >= 1)
			std::cout << "matrix.channels() == 1\n";
		cv::threshold(matrix, matrix, config.minVal, 255, cv::THRESH_BINARY); // 255 is just dummy
	}
	else {
		if (Execution::verbosity >= 1)
			std::cout << "matrix.channels() == 3\n";
		std::vector<cv::Mat> matplanes;
		cv::split(matrix, matplanes);
		for (int c = 0; c < matrix.channels(); c++){
			cv::threshold(matplanes[c], matplanes[c], config.minVal, 255, cv::THRESH_BINARY); // 255 is just dummy
		}
		cv::merge(matplanes, matrix);
	}

	// write the output-matrix = binary image
	if(config.writeMatrix)
	{
		writeMatrix(matrix); //, chan);
	}

	if (config.writeOutput){
		// the real marking step
		std::vector<cv::Mat> planes;
		cv::split(markImg, planes);
		uchar** rows = new uchar*[markImg.channels()];
		bool cl;
		for (int y = 0; y < markImg.rows; y++){
			for (int p = 0; p < markImg.channels(); ++p) {
				rows[p] = planes[p].ptr<uchar>(y);
			}
			for (int x = 0; x < markImg.cols; x++){

				if (matrix.channels() == 1)
					cl = (matrix.at<uchar>(y,x));
				else
					cl = (matrix.at<cv::Vec3b>(y,x)[0] || matrix.at<cv::Vec3b>(y,x)[1] || matrix.at<cv::Vec3b>(y,x)[2]);

				if (markImg.channels() >= 3) {
					// here we change *undetected* pixels to better expose
					// the detected ones
					if (!cl){ 
						float mix = (float)rows[0][x] + (float)rows[1][x]
							+ (float)rows[2][x];
						mix = round(0.333333f*0.75f * mix);
						rows[0][x] = rows[1][x] = rows[2][x] = mix;
					}
				} else {
					// here we can only mark in an inferior way
					if (cl) {
						rows[0][x] /= 2;
					}
				}
			}
		}
		cv::merge(planes, markImg);
		origimg = markImg;
		delete[] rows;
	}
}

void Mark ::labelRegions(
    std::vector<std::vector<std::pair<cv::Point, cv::Point> > > &clusters,
   const cv::Mat &post) {

  cv::Mat_<bool> visited(post.rows, post.cols, false);
  cv::Mat_<int> detectionMapIds(post.rows, post.cols, 0);
  int max_label = 0;
  int label;

  for (size_t g = 0; g < clusters.size(); g++) {
    // for every cluster in group
    cv::Scalar color = g;
    for (size_t c = 0; c < clusters[g].size(); c++) {
      cv::Point p1 = clusters[g][c].first;
      cv::Point p2 = clusters[g][c].second;
      // draw small circle
      if (p1.x > 0 && p2.x > 0 && post.at<uchar>(p1.y,p1.x) &&
          post.at<uchar>(p2.y,p2.x)) {
        // Perform labeling with neighborhood of p1
        if (!visited(p1.y, p1.x)) {
          if (detectionMapIds(p2.y, p2.x))
            label = detectionMapIds(p2.y, p2.x);
          else
            label = ++max_label;

		  // Perform label propagation
		  // Create a stack and insert non visited points
          std::queue<cv::Point> stack;
          stack.push(p1);
          if (!visited(p2.y, p2.x))
            stack.push(p2);

          while (!stack.empty()) {
            cv::Point p = stack.front();
            stack.pop();
			
            if(!visited(p.y,p.x)){
              visited(p.y, p.x) = true;
			        detectionMapIds(p.y,p.x) = label;

              std::vector<cv::Point> neighbors =
                  getNeighbors(visited, p );

              // Check if the neighboor was detected in post-processing
			  // If yes, insert it in the stack, otherwise ignore it
              for (size_t i = 0; i < neighbors.size(); i++) {
				if ( !visited(neighbors[i].y, neighbors[i].x)
			         && post.at<uchar>(neighbors[i].y, neighbors[i].x) != 0 ){
                  stack.push(neighbors[i]);
                }
              }
			}
          }
        }
      }
    }
  }
	//Write Final image with label
    dumpMatrix("_labeled", detectionMapIds);
}


// TODO put that in the markIt function and use markImg not origimg
void Mark :: markRegions(std::vector<std::vector
						  <std::pair<cv::Point,cv::Point> > > & clusters,
						  bool draw_shift_vecs)
{
	if (Execution::verbosity >= 1){
		std::cout << "mark regions with " << (draw_shift_vecs ? "" : "no ") << "shift vecs " << std::endl;
		std::cout << "num of clusters: " << clusters.size() << std::endl;
		// for every group
		if (Execution::verbosity >= 2){
			for  (size_t g = 0; g < clusters.size(); g++){
				std::cout << "\tcluster[" << g << "].size() = " << clusters[g].size() << std::endl;
			}
		}
	}
	// output image
	cv::Mat m;
	if (config.useOrig){
		m = origimg.clone();
	}
	else {
		m = img.clone();
	}

	// define colors
	std::vector<cv::Scalar> color_pool;
	color_pool.push_back(CV_RGB(255,0,0));
	color_pool.push_back(CV_RGB(0,0,255));
	color_pool.push_back(CV_RGB(255,0,255));
	color_pool.push_back(CV_RGB(0,255,255));
	color_pool.push_back(CV_RGB(255,255,0));

	cv::Scalar shift_color = CV_RGB(0,255,0); // green

//	cv::Point block_shift(0,0);
	cv::Point block_shift(block_size/2, block_size/2);

	// for every group
	for  (size_t g = 0; g < clusters.size(); g++){
		// for every cluster in group
		cv::Scalar color = color_pool[g % color_pool.size()];
//		cv::Mat tmp = origimg.clone();
		for (size_t c = 0; c < clusters[g].size(); c++){
			cv::Point p1 = clusters[g][c].first;
			// draw small circle
			if (p1.x > 0){
				cv::circle(m, p1 + block_shift, 1, color);
//				cv::circle(tmp, p1 + block_shift, 1, color);
			}
			cv::Point p2 = clusters[g][c].second;
			if (p2.x > 0){
				//cv::circle(m, p2 + block_shift, 1, color);
				if (cv::Scalar(m.at<cv::Vec3b>(p2.y, p2.x)) != color_pool[0] &&
						cv::Scalar(m.at<cv::Vec3b>(p2.y, p2.x)) != color_pool[1])
					cv::circle(m, p2 + block_shift, 1, cv::Scalar(0,127,255));
//				cv::circle(tmp, p2 + block_shift, 1, color);
			}
			if (draw_shift_vecs && p1.x > 0 && p2.x > 0){
				cv::line(m, p1 + block_shift, p2 + block_shift, shift_color);
			}
		}
//		std::stringstream ss;
//		ss << Execution::outputdir << "/" << Execution::image_basename << "_regions" << Execution::suffix << "_" << g << ".png";
//		cv::imwrite(ss.str(), tmp);
	}

	dumpMatrix("_regions", m);
}

void Mark :: markContour(const cv::Mat & cmfd_map,
						 std::vector<std::vector<cv::Point> > * contours,
						 std::vector<cv::Vec4i> * hierarchy)
{
	// output image
	cv::Mat m;
	if (config.useOrig){
		m = origimg.clone();
	}
	else {
		m = img.clone();
	}

	cv::Mat_<uchar> cmfd;
	if (cmfd_map.channels() == 3){
		std::vector<cv::Mat> planes;
		cv::split(cmfd_map, planes);
		cmfd = planes[0] + planes[1] + planes[2];
	}
	else{
		// we have to clone here as the image will be modified in
		// findContours!!
		cmfd = cmfd_map.clone();
	}
	//cv::imwrite("contour_img", orig_image);
	if (contours == NULL || contours->empty() ){
		contours = new std::vector<std::vector<cv::Point> >();
		hierarchy = new std::vector<cv::Vec4i>();
		cv::findContours(cmfd, *contours, *hierarchy,
						 CV_RETR_EXTERNAL, // no inner contours
						 CV_CHAIN_APPROX_NONE); // should the chain be approximated
	}

	cv::drawContours(m, *contours, -1,
//					 cv::Scalar(0,127,255), // color orange
					 cv::Scalar(0,255,0), // color green
					 7 // thickness, if CV_FILLED: filled out
					 );

	dumpMatrix("_contour", m);
}

void Mark :: markHull(const std::vector<std::vector<std::pair<cv::Point,cv::Point> > > & clusters)
{
	if (clusters.size() == 0){
		vout(2) << "no matches for a convex hull\n";
		return;
	}

	// output image
	cv::Mat m;
	if (config.useOrig){
		m = origimg.clone();
	}
	else {
		m = img.clone();
	}

	// overall hull
	std::vector<cv::Point> hull;
	for ( size_t i = 0; i < clusters.size(); i++ )
	{
		std::vector<cv::Point> cluster_single;
		std::vector<cv::Point> cluster_hull;
		for ( size_t k = 0; k < clusters[i].size(); k++ )
		{
			cluster_single.push_back( clusters[i][k].first );
		}
		cv::convexHull( cluster_single, cluster_hull );
		// insert the hull of one cluster
		hull.insert( hull.end(),
					 cluster_hull.begin(),
					 cluster_hull.end() );
	}

	// dump hull
	cv::Mat_<uchar> hull_m(m.rows, m.cols, (uchar) 0);
	for (size_t i = 0; i < hull.size(); i++)
	{
		hull_m( hull[i].y, hull[i].x ) = 255;
	}

	dumpMatrix("_hull", hull_m);

	return;
}

// marks it in a different way -> shift vecs connected by blue lines
void Mark :: markIt(std::vector< std::multimap<std::pair<int,int>,
					std::pair<cv::Point, cv::Point> > > shift)
{
	cv::Mat m;
	if (config.useOrig){
		m = origimg;
	}
	else {
		m = img;
	}
	for (size_t c = 0; c < shift.size(); c++){
		std::multimap<std::pair<int,int>, std::pair<cv::Point, cv::Point> >::iterator it;

		for (it = shift[c].begin(); it != shift[c].end(); ++it){
			cv::line(m, it->second.first, it->second.second, CV_RGB(0,0,255));
		}
	}
	writeMatrix(m);
}

void Mark :: writePairs(const std::string & groundTruthFile, bool shift)
{
	cv::Mat_<uchar> groundTruth = cv::imread(groundTruthFile, 0);
	if (!groundTruth.data){
		std::cerr << "WARNING: Can't write pairs: Mark::writePairs(): no ground Truth File given\n";
		return;
	}
	if (corres[0].size() != groundTruth.size()){
		std::cerr << "WARNING: can't generate pair-file as groundtruthfile-size differs with matrix size - perhaps decomposed?\n";
		return;
	}
	cv::Mat_<uchar> detection(origimg.rows, origimg.cols, (uchar)0);
	cv::Mat_<cv::Vec3b> detection_visual(origimg.rows, origimg.cols, cv::Vec3b(0,0,0));
	unsigned int count = 0;
	unsigned int countall = 0;
	unsigned int countallmatch= 0;
	bool hit1 = false;
	bool hit2 = false;
	for (int y = 0; y < origimg.rows - block_size/2; ++y) {
		for (int x = 0; x < origimg.cols - block_size/2; ++x) {
			cv::Point2f corr = Matching::getCorres(corres[0], x, y);
			if (corr.x < 0) {
				countallmatch++;	
				continue;
			}
			if (shift){ // shift by blockSize / 2
				if (groundTruth(y+block_size/2, x+block_size/2) == 255){
					hit1 = true;
				}
				if (groundTruth(corr.y+block_size/2, corr.x+block_size/2) == 255){
					hit2 = true;
				}
			}
			else {
				if (groundTruth(y,x) == 255){
					hit1 = true;
				}
				if (groundTruth(corr.y, corr.x) == 255){
					hit2 = true;
				}
			}
			if (hit1 && hit2)
			{
				// mark hits as red
				detection[y + block_size/2][x + block_size/2] = 255;
				detection[(int)corr.y + block_size/2][(int)corr.x + block_size/2] = 255;

				detection_visual[y + block_size/2][x + block_size/2] = cv::Vec3b(0,0,255); // red
				detection_visual[(int)corr.y + block_size/2][(int)corr.x + block_size/2] = cv::Vec3b(0,0,255);
				count++;
			}
			else {// otherwise mark them as gray
				detection_visual[y + block_size/2][x + block_size/2] = cv::Vec3b(128,128,128); // gray
			}
			if (hit1 || hit2)
				countall++;
			hit1 = false;
			hit2 = false;
		}
	}
	if (Execution::verbosity >= 1){
		std::cerr << "correspondences which would theoretically lie in correct region: " << count
				  << "\t(Note: this is not the number of correct marked pairs)"
				  << std::endl;
	}
	
	std::stringstream ss1, ss2;
	ss1 << Execution::outputdir << "/" << Execution::image_basename << "_pairs" << Execution::suffix << ".png";
	ss2 << Execution::outputdir << "/" << Execution::image_basename << "_pairs_visual" << Execution::suffix << ".png";
	if ( ! cv::imwrite(ss1.str(), detection) ){
		std::cerr << "couldn't write " << ss1.str() << std::endl;
	}
	else if (Execution::verbosity >= 1){
		std::cout << " wrote " << ss1.str() << std::endl;
	}
	if ( ! cv::imwrite(ss2.str(), detection_visual) ) {
		std::cerr << "couldn't write " << ss2.str() << std::endl;
	}
	else if (Execution::verbosity >= 1){
		std::cout << " wrote " << ss2.str() << std::endl;
	}
}

void Mark :: writeChrom(void)
{
	cv::Mat_<cv::Vec3b> bla(img.rows, img.cols, cv::Vec3b(0, 0, 0));
	int maxDim = std::max(img.rows, img.cols);
	for (int y = 0; y < img.rows; ++y) {
		for (int x = 0; x < img.cols; ++x) {
			cv::Point2f tmp = Matching::getCorres(corres[0], x, y);
			if (tmp.x < 0) {
				continue;
			}
			/* 
			//doesnt work yet
			int32_t x = static_cast<int32_t>(tmp.x);
			int32_t y = static_cast<int32_t>(tmp.y);
			uchar b,g,r;
			b = g = r = 0;
			b = static_cast<uchar>(y & 0x000000ff);
			int32_t lastByte = 15;
			g |= static_cast<uchar>(x & lastByte);
			g <<= 4;
			int32_t tmp2 = (y >> 8);
			g |= static_cast<uchar>(tmp2 & lastByte);
			tmp2 = x;
			tmp2 >>= 4;
			tmp2 &= 0x000000ff;
			r = static_cast<uchar>(tmp2);
			bla(y,x) = cv::Vec3b(b,g,r);
			*/
			
			// sei tmp.x/(dimx+dimy) = r/(r + g + b), tmp.y/(dimx + dimy) = b/(r+g+b)
			float chromr, chromg, chromb;
			chromr = tmp.x/maxDim;
			chromb = tmp.y/maxDim;
			if (chromr + chromb > 1) {
				chromr /= (chromr + chromb);
				chromb /= (chromr + chromb);
				chromg = 0;
			} else {
				chromg = 1.0 - chromr - chromb;
			}
			bla[y][x] = cv::Vec3b(chromb * 255, chromg * 255, chromr * 255);
		}
	}
	std::stringstream ss;
	ss << Execution::outputdir << "/" << Execution::image_basename << "_chrom" << Execution::suffix << ".png";
	if (! cv::imwrite(ss.str(), bla, params) ){
		std::cerr << "couldn't write " << ss.str() << std::endl;
	}
	else if (Execution::verbosity >= 1){
		std::cout << " wrote " << ss.str() << std::endl;
	}
}

// If we have computed a transformation with ransac or sats
// build up correleation maps: compare each
// pixel to its transformation, mark those pixels which have a high
// correlation of its neighborhood and the neighborhood around the
// transformed pixel location
std::vector<cv::Mat_<uchar> > Mark :: computeCorrelation(const std::vector<cv::Mat> & trans, bool write)
{
	vout(1)  << "Mark::writeCorrelation: have trans.size() " << trans.size() << std::endl;
	if (img.empty()){
		throw std::runtime_error("Mark :: writeCorrelation: if you want to build a correlation-map,"
					" you need to have an img, specify it with -I !\n");
	}

	// the output correlation map
	cv::Mat_<uchar> correlation_map_merged(img.rows, img.cols, (uchar)0);
	std::vector<cv::Mat_<uchar> > correlation_map;
	if (trans.empty()){
		return correlation_map;
	}
	// split image into planes
	std::vector<cv::Mat> planes;
	cv::split(img, planes);

	for (size_t i = 0; i < trans.size(); i++)
	{
		if (trans[i].empty()){
			continue;
		}
		cv::Mat_<double> trafo = trans[i];
//		if (trafo.rows == 3 && trafo(2,0) == 0.0 && trafo(2,1) == 0.0 && trafo(2,2) == 1.0){
//			trafo = trafo( cv::Range(0,2), cv::Range(0,3) );
//		}

//		int shift_x = static_cast<int>(trafo(0, 2));
//		int shift_y = static_cast<int>(trafo(1, 2));
//		double sx = std::sqrt( trafo(0,0)*trafo(0,0) + trafo(0,1)*trafo(0,1) );
//		if (trafo(0,0) < 0) sx *= -1;
//		double sy = std::sqrt( trafo(1,0)*trafo(1,0) + trafo(1,1)*trafo(1,1) );
//		if (trafo(1,1) < 0) sy *= -1;
//		double angle = std::atan(trafo(1,0) / trafo(1,1));
//		// radians -> degree, acos(-1.0) == pi
//		angle *= 180 / std::acos(-1.0);

//		// some output
//		if (Execution::verbosity >= 2){
//			std::cout  << "transformation matrix:\n" << trafo << std::endl;
//			std::cout << "measured shift: x: " << shift_x << " y: " << shift_y << std::endl;
//			std::cout << "measured scale: sx: " <<  sx << " sy: " << sy << std::endl;
//			std::cout << "measured angle: " << angle << std::endl;
//		}

// Note: we do that now in the verification class
//		// ignore abstruse scales
//		if ( ( (trafo.rows == 2)
//			 || ( (trafo.rows == 3) && (trafo(2,0) == 0.0
//									&& trafo(2,1) == 0.0 && trafo(2,2) == 1.0)) )
//			 && (std::abs(sx) > 1.6 || std::abs(sx) < 0.4 || std::abs(sy) > 1.6 || std::abs(sy) < 0.6) )
//		{
//			vout(2) << "skip as scale is too abstruse" << std::endl;
//			continue;
//		}

		for (int c = 0; c < img.channels(); c++){
			// compute correlation maps
			cv::Mat_<float> correl_f = computeCorrelMap(planes[c], trafo);
			if (correl_f.empty())
				std::cerr << "WARNING correl_f is empty\n";

			// compute 2nd correlation map which uses the inverse matrix
			cv::Mat_<float> correl_b = computeCorrelMap(planes[c], trafo, true);

			std::stringstream ss;
			ss.clear();
			ss << c << "_" << i;

			// dump the correlation maps
			if (Execution::verbosity >= 3){
				dumpMatrix(std::string("_correl_f_" + ss.str()), 255*correl_f);
				dumpMatrix(std::string("_correl_b_" + ss.str()), 255*correl_b);
			}

			// Apply Gaussian filter w. kernel size 7
			if (Execution::verbosity >= 2)
				std::cout  << "Apply Gaussian Filter on correlation map\n";
			cv::Mat_<float> correl_f_gauss;
			cv::Mat_<float> correl_b_gauss;
			cv::GaussianBlur(correl_f, correl_f_gauss, cv::Size(7,7), 0.0);
			cv::GaussianBlur(correl_b, correl_b_gauss, cv::Size(7,7), 0.0);
			// dump
			if (Execution::verbosity >= 3){
				dumpMatrix(std::string("_correl_gauss_f_" + ss.str()), 255*correl_f_gauss);
				dumpMatrix(std::string("_correl_gauss_b_" + ss.str()), 255*correl_b_gauss);
			}

			// Binarize image
			cv::Mat_<float> correl_f_bin;
			cv::Mat_<float> correl_b_bin;
			if (Execution::verbosity >= 2)
				std::cout  << "Binarize correlation map with th: " << config.binaryTh << std::endl;

			cv::threshold(correl_f_gauss, correl_f_bin, config.binaryTh, 1.0, cv::THRESH_BINARY);
			cv::threshold(correl_b_gauss, correl_b_bin, config.binaryTh, 1.0, cv::THRESH_BINARY);

			cv::Mat_<uchar> correl_f_bin_u;
			cv::Mat_<uchar> correl_b_bin_u;
			correl_f_bin.convertTo(correl_f_bin_u, CV_8UC1, 255);
			correl_b_bin.convertTo(correl_b_bin_u, CV_8UC1, 255);
			// dump
			if (Execution::verbosity >= 3){
				dumpMatrix(std::string("_correl_bin_f_" + ss.str()), correl_f_bin_u);
				dumpMatrix(std::string("_correl_bin_b_" + ss.str()), correl_b_bin_u);
			}

			// combined binarized correlation map
			correlation_map_merged += correl_f_bin_u;
			correlation_map_merged += correl_b_bin_u;
			correlation_map.push_back(correl_f_bin_u);
			correlation_map.push_back(correl_b_bin_u);
		}
	}

	// write final correlation matrix
	if (write){
		dumpMatrix("_correl", correlation_map_merged);
	}

	return correlation_map;
}

// the implementation differs from the formula of Pan in 2 points:
//	1. Pans formula lacks in a sqrt in the denominator
//	2. We have to use zero-mean data, otherwise the results are too large
//     this seems to be another bug of Pan!
cv::Mat_<float> Mark :: computeCorrelMap(const cv::Mat_<uchar> &plane,
										 const cv::Mat &trans,
										 bool inverse)
{
	if (Execution::verbosity >= 1)
		std::cout  << "--> Compute correlation map\n";
	if (plane.empty() || trans.empty()){
		std::cerr << "WARNING: plane is empty -> return empty map\n";
		return cv::Mat_<float>::zeros(plane.rows, plane.cols);
	}

	// warp the image
	cv::Mat_<uchar> warped;
	if (trans.rows == 2){
		if (inverse)
			cv::warpAffine(plane, warped, trans, img.size(), cv::WARP_INVERSE_MAP + cv::INTER_LINEAR);
		else
			cv::warpAffine(plane, warped, trans, img.size());
	}
	else if (trans.rows == 3){
		if (inverse) {
			cv::warpPerspective(plane, warped, trans, img.size(),
								cv::WARP_INVERSE_MAP + cv::INTER_LINEAR);
		} else {
			cv::Mat_<double> t2 = trans;
			cv::warpPerspective(plane, warped, t2, img.size());
		}
	}

	// the output correlation map, initialized w. 0
	cv::Mat_<float> correl_map(plane.rows, plane.cols, 0.0);


	// TODO: perhaps treat borders differently
	// TODO: this could heavily speeded up by using dft + integral images,
	//		 see: http://scribblethink.org/Work/nvisionInterface/nip.html
	// Let's do it partially:
	cv::Mat_<int> plane_int, warped_int;
	cv::integral(plane, plane_int, CV_32S);
	cv::integral(warped, warped_int, CV_32S);

	bool zero_mean = true; // this should be true, otherwise values too large!
	int area_dim = 2; // --> kernel-size = 5
	int area_size = (2*area_dim+1) * (2*area_dim+1);
	for (int y = area_dim; y < (plane.rows - area_dim); y++) {
		for (int x = area_dim; x < (plane.cols - area_dim); x++)
		{
			float sum_area = 0;
			float sum_area_corres = 0;
			float sum_area_both = 0;

			// This code is slower than mine and the result also differs
			// (in my opinion in a bad way, too)
			/*
			cv::Mat_<float> result;
			cv::matchTemplate(area, area_corres, result, CV_TM_CCOEFF_NORMED);
			//if (result(0,0) != 1) // need this as dividing by zero (see below) produces somehow 1
			//	correl_map(y,x) = result(0,0);
			*/

			float mean = 0.0;
			float mean_corres = 0.0;
			if (zero_mean){
				// compute the mean with the help of the integral images in a fast manner

				// top-left - top-right - bottom-left + top-right
				int sum = plane_int(y-area_dim,x-area_dim) // top-left
						- plane_int(y-area_dim,x+area_dim+1) // top-right
						- plane_int(y+area_dim+1,x-area_dim) // bottom left
						+ plane_int(y+area_dim+1,x+area_dim+1); // bottom-right
				mean = static_cast<float>(sum) / area_size;

				// top-left - top-right - bottom-left + top-right
				int sum_corres = warped_int(y-area_dim,x-area_dim) // top-left
						- warped_int(y-area_dim,x+area_dim+1) // top-right
						- warped_int(y+area_dim+1,x-area_dim) // bottom left
						+ warped_int(y+area_dim+1,x+area_dim+1); // bottom-right
				mean_corres = static_cast<float>(sum_corres) / area_size;

				//area_f -= cv::mean(area);
				//area_corres_f -= cv::mean(area_corres);
			}
			// this is nearly 4x slower?!
			/*
			sum_area_both = cv::sum(area_f.mul(area_corres_f))[0];
			sum_area = cv::sum(area_f.mul(area_f))[0];
			sum_area_corres = cv::sum(area_corres_f.mul(area_corres_f))[0];
			*/

			// compute sums in 5x5 neighborhood around point and corresponding point
			for (int t = -area_dim; t <= area_dim; t++){
				for (int s = -area_dim; s <= area_dim; s++){
			//for (int t = 0; t < area.rows; t++) {
				//for (int s = 0; s < area.cols; s++){
					//double pixel = area(t,s);
					//double pixel_corres = area_corres(t,s);
					float pixel = static_cast<float>( plane(y+t, x+s) )  - mean;
					float pixel_corres = static_cast<float>( warped(y+t, x+s) ) - mean_corres;

					sum_area_both += pixel * pixel_corres;
					sum_area += pixel * pixel;
					sum_area_corres += pixel_corres * pixel_corres;
				}
			}

			// avoid dividing by zero
			if ( sum_area * sum_area_corres > DBL_EPSILON ){
				// INFO: ERROR in the Pan paper! there must be a square root around the denominator
				float correlation = sum_area_both / std::sqrt(sum_area * sum_area_corres);
				correl_map(y, x) = correlation;
			}
		}
	}

	if (Execution::verbosity >= 1)
		std::cout  << "<-- compute correlation map\n\n";
	return correl_map;
}
void Mark :: dumpMatrix(std::string identifier, const cv::Mat & dump_matrix)
{
	// name
	std::stringstream ss;
	//ss << Execution::outputdir << "/" << Execution::image_basename
	ss << Execution::outputdir  << Execution::image_basename
	   << identifier << Execution::suffix << ".png";

	// png compression level
	std::vector<int> parameters;
	parameters.push_back(CV_IMWRITE_PNG_COMPRESSION);
	parameters.push_back(9); // highest compression

	if (! cv::imwrite(ss.str(), dump_matrix, parameters) ){
		std::cerr << "couldn't write " << ss.str() << std::endl;
	}
	else if (Execution::verbosity >= 1){
		std::cout << " wrote " << ss.str() << std::endl;
	}
}

