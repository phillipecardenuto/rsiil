/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include <algorithm>
#include <iostream> 
#include <queue>

#include "highgui.h"
#include "feature.h"
#include "matching.h"
#include "cps.h"
#include "execution.h"
#include "stdio.h"
#include "blockhandling.h"

Matching :: Matching(const std::vector<cv::KeyPoint> & _keypoints,
					 const cv::Mat & _features,
					 const cv::Size & _img_size,
					 struct mat_cfg config)
	: keypoints(_keypoints), features(_features), img_size(_img_size)
{
	this->config = config;
	// prepare criterion
	crit = 0;
	for (size_t i = 0; i < config.crit.size(); i++){
		int tmp = config.crit[i] - 1;
		// this pow we need cause boost::options dont save the values actually
		// specified in the enum but the position in the enum
		tmp = pow((double)2,(double)tmp);
		crit += tmp;
		critToId[tmp] = i;
		// test:
		//std::cout << EUCLIDIAN << " : " << critToId[EUCLIDIAN] << std::endl;
	}
	// check if user gave a threshold
	if (crit & EUCLIDIAN){
		if (config.th.empty()){
			throw std::runtime_error("no threshold set for euclidian criterion");
		}
		if (!config.euclidianTh){
			vout(2) << "adjust euclidian th\n";
			// this is how Mahdian et.al has defined his threshold
			// if the threshold is a similarity th (s. mahdian et al)
			// s(Bi, Bj) = 1 / (1 + ro(Bi, Bj) with ro = sqrt(squared euclid_dist)
			double t = config.th[critToId[EUCLIDIAN]];
			config.th[critToId[EUCLIDIAN]] = pow( (double)(1.0 / config.th[0]) - 1.0, (double)2);
			vout(2) << "similary threshold T = " << t
					<< " --> ro(Bi,Bj) <= " << config.th[critToId[EUCLIDIAN]] << std::endl;

		}
	}
	if (Execution::verbosity >= 2){
		if (crit & EUCLIDIAN_RATIO){
			std::cout << "EUCLIDIAN-RATIO-th: " << config.th[critToId[EUCLIDIAN_RATIO]] << std::endl;
		}
		if (crit & EUCLIDIAN){
			std::cout << "EUCLIDIAN-th: " << config.th[critToId[EUCLIDIAN]] << std::endl;
		}
	}
	// initial mask with true
	mask = cv::Mat_<uchar>(features.rows, 1, 1);
	if (Execution::verbosity >= 2){
		std::cout << "Matching::Matching(): size: " << img_size.width << "x" << img_size.height << std::endl;
	}
	// we need to have at least one corres-map
	corres.push_back(cv::Mat_<cv::Point>(img_size.height, img_size.width, cv::Point(-1,-1)));
}

Matching :: ~Matching(void)  
{
//	delete ind; // somehow this doesn work...
}

bool Matching ::euclidianDistanceCriterion(const cv::Mat_<float> & first, const cv::Mat_<float> & second) const
{
	double diff = cv::norm(first, second, cv::NORM_L2);
	if ( diff > config.th[critToId.find(EUCLIDIAN)->second] ){
		vout(4) << "eliminate one match" << std::endl;
		return false;
	}
	return true;
}

bool Matching :: LuoCriterion(const cv::Mat_<float> & first, const cv::Mat_<float> & second) const 
{
	std::vector<double> ths = config.th;
	if ((ths.size() - 2) < static_cast<unsigned int>(first.cols))
	{
		throw std::runtime_error("Matching::LuoCriterion: feature.dim > num of thresholds (-2) wrong");
	}
	CV_Assert(first.rows == 0);
	CV_Assert(first.rows == second.rows);

	double diff123 = 0.0;
	double diff4567 = 0.0;
	for (int k = 0; k < (int)ths.size() - 2; k++){ // for every threshold
		double diff = fabs(second(0,k) - first(0,k)); // need abs cause k > 0 could be unsorted
		if (diff > ths[k]){
			return false;
		}
		if (first.cols == 7){
			if (k <= 2)
				diff123 += diff;
			if (k == 2 && diff123 > ths[7])
				return false;
			else 
				diff4567 += diff;
		}
		else {
			if (k==0 && diff > ths[5])
				diff123 += diff;
			else {
				diff4567 += diff;
			}
		}
	}
	if (diff4567 > ths[8]){
		return false;
	}
	return true;
}

bool Matching :: cpsCriterion(const cv::Mat_<float> & first, const cv::Mat_<float> & second) const
{
	cv::Mat_<double> f1 = first;
	cv::Mat_<double> f2 = second;
	cv::Mat p;
	Cps::crossPowerSpectrum(f1, f2, p, true);
	double max;
    cv::minMaxLoc(p, NULL, &max, NULL, NULL);
	if (max < config.phaseCorrelTh)
		return false;
	return true;
}


bool Matching :: BravoCriterion(const cv::Mat_<float> & first, const cv::Mat_<float> & second) const
{
	CV_Assert(first.rows == 0);
	CV_Assert(first.cols == 4);
	// check avg colors
	for (int i = 1; i < 3; i++){
		double diff = fabs(first(0,i) - second(0,i));
		if (diff > config.th[0])
			return false;
	}
	// check entropy
	double diff2 = fabs(first(0,3) - second(0,3));
	if (diff2 > config.th[1])
		return false;

	// check correlation coefficient
	return correlCriterion(first, second);
}

// check the correleation coefficient of the fourier magnitudes of the log polar representation
bool Matching :: correlCriterion(const cv::Mat_<float> & first, const cv::Mat_<float> & second) const
{
	// compute log polar representation
	// true -> center, 180deg=rows, 8 width
	cv::Mat_<double> f1 = BlockHandling::logPolar(first, true, 180, 8);
	cv::Mat_<double> f2 = BlockHandling::logPolar(second, true, 180, 8);

	cv::Mat F1;
	cv::Mat F2;
	// compute fourier trafo
	cv::dft(f1, F1, cv::DFT_COMPLEX_OUTPUT);
	cv::dft(f2, F2, cv::DFT_COMPLEX_OUTPUT);

	// compute the complex magnitude
	std::vector<cv::Mat> realcomp1;
	std::vector<cv::Mat> realcomp2;
	// split in real and complex parts
	cv::split(F1, realcomp1); 	
	cv::split(F2, realcomp2); 
	// use only left part cause of DFT_COMPLEX_OUTPUT gains 1 array with 2 chans, 
	// every only half computed
	cv::Mat real1(realcomp1[0], cv::Range::all(), cv::Range(0, F1.cols/2 + 1));
	cv::Mat real2(realcomp2[0], cv::Range::all(), cv::Range(0, F2.cols/2 + 1));
	cv::Mat comp1(realcomp1[1], cv::Range::all(), cv::Range(0, F1.cols/2 + 1));
	cv::Mat comp2(realcomp2[1], cv::Range::all(), cv::Range(0, F2.cols/2 + 1));
	// compute complex magnitude
	cv::Mat mag1;
	cv::Mat mag2;
	cv::magnitude(real1, comp1, mag1); 
	cv::magnitude(real2, comp2, mag2); 

	if (mag1.cols > 1 && mag1.rows > 1)
		throw std::runtime_error("Matching::CorrelCriterion::mag1.cols > 1 && mag1.rows > 1");
	if (mag2.cols > 1 && mag2.rows > 1)
		throw std::runtime_error("Matching::CorrelCriterion::mag2.cols > 1 && mag2.rows > 1");

	
	double nom = mag1.dot(mag2);
	double denom1 = mag1.dot(mag1);
	double denom2 = mag2.dot(mag2);

	double correl = nom / sqrt( denom1 * denom2);

	if (Execution::verbosity >= 3) {
		// Note: getPosi gives left upper corner of the block
		cv::Point pa = BlockHandling::getPosi(first);
		cv::Point pb = BlockHandling::getPosi(second);
		std::cout << nom << " / sqrt( " << denom1 << " " << denom2 << " )\n";
		std::cout << "correlation coeff of block(" << pa.x << ", " << pa.y <<") and b("
				<< pb.x << ", " << pb.y << ") : "<< correl << std::endl;
	}
	
	if (correl < config.correlTh)
		return false;
	return true;
}

bool Matching :: chebyshevDistance(int dx, int dy) const 
{
	// already abs
	if (std::max(dx, dy) < config.minDistCheby) 
		return false;
	return true;
}

bool Matching :: euclidianDistance(int dx, int dy) const 
{
	if ((sqrt(dx*dx + dy*dy)) < config.minDistEuclidian){
		return false;
	}
	return true;
}

void Matching :: checkAndPut(const cv::Mat_<float> & first, const cv::Mat_<float> & second, 
							cv::Point pa, cv::Point pb)
{
	int dx, dy;
	if (checkDistAndCrit(first, second, pa, pb, &dx, &dy)){
		// insert in vector of correspondence matrices
		int i = 0;
		while(true){
			// corres.size() must > i if not add another plane
			if ((int)corres.size() == i){
				corres.push_back(cv::Mat_<cv::Point>(img_size.height,
													 img_size.width,
													 cv::Point(-1,-1)));
			}
			// add correspondence for that point
			if (corres[i](pa.y,pa.x).x == -1){
				corres[i](pa.y,pa.x) = pb;
				break;
			}
			// it may happen that we selected too many points in the first place
			// e.g. SIFT produces somehow more than 1 keypoint at the same position
			if (i == (config.numRows-1)){
				break;
			}

			// we haven't inserted it yet -> try with next plane
			i++;
		}
	}
}

bool Matching :: checkDistAndCrit(const cv::Mat_<float> & first, const cv::Mat_<float> & second,
		cv::Point & pa, cv::Point & pb,
		int *dxp, int *dyp)
{
	// abs for shift vectors!
	int dx = std::abs(pa.x - pb.x);
	int dy = std::abs(pa.y - pb.y); 

	if (dxp != NULL)
		*dxp = dx;
	if (dyp != NULL)
		*dyp = dy;

	// positions far enough away from each other?
	if (config.minDistEuclidian > 0) {
		if (!euclidianDistance(dx, dy))
			return false;
	}
	if (config.minDistCheby > 0) {
		if (!chebyshevDistance(dx, dy))
			return false;
	}

	// check similarity
	bool putin = true;
	
	if (crit & LUOCRIT){
		putin = LuoCriterion(first, second);
	}
	if (crit & EUCLIDIAN){
		putin = euclidianDistanceCriterion(first, second);
	}
	if (crit & CPS) {
		putin = cpsCriterion(first, second);
	}
	if (crit & CORREL) {
		putin = correlCriterion(first, second);
	}
	if (crit & BRAVOCRIT) {
		putin = BravoCriterion(first, second);
	}

	return putin;
}

void Matching :: computeMatchingKD(void)
{
	int N = features.rows; //size();
	int dim = features.cols; //getFeatureSize(c);
	if (Execution::verbosity >= 1){
		std::cout << " N = " << N << std::endl;
		std::cout << " dim = " << dim << std::endl;
	}
	if (Execution::verbosity >= 1){
		std::cout << " Initialize FLANN w. version " << FLANN_VERSION_ << std::endl;
	}
	double t = (double)cv::getTickCount();

	cv::setNumThreads(0); // does work if openmpi is linked
	//	ind = new cv::flann::Index(data, cv::flann::AutotunedIndexParams() );
	//	ind = cv::flann::Index(data, cv::flann::AutotunedIndexParams(0.99, 0.01, 0, 0.5) );

	ind = new cv::flann::Index( features, cv::flann::KDTreeIndexParams() );


	if (Execution::verbosity >= 1){
		std::cout << " get closest features" << std::endl;
	}

	// get numRows results around point i
	if (config.numRows == -1){
		if (Execution::verbosity >= 1){
			std::cout << "config.numRows == -1 -> set to 10";
		}
		config.numRows = 20; // just to avoid strange behaviour
	}
	
	// initialize with -1 to see if we get indices
	cv::Mat_<int> indices(N, config.numRows + 1, -1); // numRows + 1 cause first columns == block itself
	cv::Mat_<float> dists(N, config.numRows + 1);

	ind->knnSearch(features, indices, dists, config.numRows + 1, cv::flann::SearchParams(64) );

	// matching
	int	anz = config.numRows;
	for (int row = 0; row < features.rows; row++) {
		if (mask(row, 0) == 0){
			continue;
		}
		cv::Mat_<float> first = features.row(row);

		int count = 0;
		for (int col = 0; col <= anz; col++){
			// it may happen that the first index is not i,
			// but in general it will be, in that case -> continue
			if ( indices(row,col) == row || indices(row, col) == -1 || !mask(indices(row,col)) ){
				continue;
			}
			cv::Mat_<float> second = features.row( indices(row,col) );
			double curr_distance = cv::norm(first, second, cv::NORM_L2);

			// checks ratio of distances between i'th and (i+1)'th matched pair
			if ( (crit & EUCLIDIAN_RATIO) && config.numRows > 1)
			{
				if ( col + 1 > anz ){ // omit the last one
					break;
				}
				cv::Mat_<float> next = features.row( indices(row, col+1) );
				double next_distance = cv::norm( second, next, cv::NORM_L2 );
				// if th is greater than the ratio we stop adding
				if (curr_distance / next_distance > config.th[critToId[EUCLIDIAN_RATIO]]){
					break;
				}
			}

			checkAndPut(first, second, 
					keypoints[row].pt,
					keypoints[indices(row,col)].pt);

			// as our loop could go further we need to break it
			count++;
			if (count == config.numRows){
				break;
			}
		}
	}

	//measure time
	t = ((double)cv::getTickCount() - t)/cv::getTickFrequency();
	if (Execution::verbosity >=2){
		std::cout << "kdtree/flann finsished. time neaded: " << t << std::endl;
	}
}

// brute force nearest neighbor
void Matching :: computeNN()
{
	double t = (double)cv::getTickCount();

	cv::Ptr<cv::DescriptorMatcher> ind;
	if (features.type() == CV_8U
			|| features.type() == CV_8S) {
		ind = cv::DescriptorMatcher::create("BruteForce-Hamming");
	}
	else{
		ind = cv::DescriptorMatcher::create("BruteForce");
	}
	std::vector<std::vector<cv::DMatch> > matches;
	ind->knnMatch(features, features, matches, config.numRows+1);

	for (size_t i = 0; i < matches.size(); i++){
		for (size_t k = 0; k < matches[i].size(); k++){
			cv::DMatch dm = matches[i][k];
			if (dm.trainIdx == dm.queryIdx ) // normally that's true for k == 0
				continue;
			cv::Mat first = features.row(dm.queryIdx);
			cv::Mat second = features.row(dm.trainIdx);

			if ( crit & EUCLIDIAN_RATIO && config.numRows > 1 ){
				if ( k+1 == matches[i].size() ) // omit last one
					break;
				float curr_distance = dm.distance;
				float next_distance = matches[i][k+1].distance;
				if ( curr_distance / next_distance > config.th[critToId[EUCLIDIAN_RATIO]] ) {
					break;
				}
			}
			checkAndPut(first, second,
						keypoints[dm.queryIdx].pt,
						keypoints[dm.trainIdx].pt);
		}
	}

/*	// own version - slower than opencv version above
	for (int row = 0; row < features.rows; row++){
		cv::Mat_<float> first = features.row(row);

		// create a heap to store the distances
		// unfortunately I can't specify the size, which would be more efficient.
		// Maybe this fact makes the whole operation slower than just sorting the rows
		// of a one-dimensional distance matrix (can't take 2-dim cause of memory)
		// TODO: check this
		std::priority_queue<
				std::pair<double, int>,
				std::vector<std::pair<double, int> >,
				std::greater<std::pair<double, int> > > prio;
		for (int k = 0; k < features.rows; k++){
			if (k == row){
				continue;
			}
			cv::Mat second = features.row(k);
			double distance = cv::norm(first, second, cv::NORM_L2);
			prio.push(std::make_pair(distance,k));
		}

		for (int l = 0; l < config.numRows; l++){
			std::pair<double,int> top = prio.top();
			prio.pop();

			// checks ratio of distances between i'th and (i+1)'th matched pair
			if (crit & EUCLIDIAN_RATIO && config.numRows > 1){
				std::pair<double, int> next = prio.top();
				double curr_distance = top.first;
				double next_distance = next.first;
				// if th is greater than the ratio we stop adding
				if ( curr_distance / next_distance > config.th[critToId[EUCLIDIAN_RATIO]] ) {
					break;
				}
			}
			checkAndPut(first, features.row(top.second),
						keypoints[row].pt, keypoints[top.second].pt);
		}
	}
*/
	//measure time
	t = ((double)cv::getTickCount() - t)/cv::getTickFrequency();
	if (Execution::verbosity >=2){
		std::cout << "nearest neighbor brute force finsished. time neaded: " << t << std::endl;
	}
}

void Matching :: computeMatchingLexi(void)
{
	// create vector of feature rows
	std::vector<int> indices;
	indices.resize(features.rows);
	for (int i = 0; i < features.rows; i++){
		indices[i] = i;
	}
	// sort indices
	std::sort(indices.begin(), indices.end(), LessThanIdx(features));

	// the actual matching
	int anz1 = config.numRows;
	int anz2 = config.numRows;
	if (config.numRows == -1)
		anz1 = 10;
	for (int i = 0; i < (features.rows - anz1); i++) // for all sorted features
	{
		if (mask(i,0) == 0)
			continue;
		cv::Mat_<float> first = features.row(indices[i]);
		if (config.numRows == -1)
			anz2 = features.rows - i - 1;  // so check every block left

		double old_distance = DBL_MAX; // for EUCLIDIAN_RATIO-crit
		// check numRows neighbors
		for (int k = 1;  k <= anz2; k++){
			if (mask(i+k,0) == 0)
				continue;
			cv::Mat_<float> second = features.row(indices[i+k]);

			if (((crit & BRAVOCRIT) || (crit & LUOCRIT)) && config.th.size() > 0){ // just to shorten the things a little bit... (works cause sorted)
				if (second(0,0) > (first(0,0) + config.th[0]))
					break;
			}
			else {
				// checks ratio between i'th and (i+1)'th feature
				if (crit & EUCLIDIAN_RATIO && config.numRows > 1){
					double curr_distance = cv::norm(first, second, cv::NORM_L2);
					if (old_distance != DBL_MAX){ // then check
						if (old_distance / curr_distance > config.th[critToId[EUCLIDIAN_RATIO]]){
							break;
						}
					}
					old_distance = curr_distance;
				}

				checkAndPut(first, second, 
						keypoints[indices[i]].pt,
						keypoints[indices[i+k]].pt );
			}
		}
	}
	// some debug infos:
	if (Execution::verbosity >= 3){
		std::cout << "print first 3 sorted features" << std::endl;
		for (int i = 0; i < 3; i++){
			Feature::printFeature(features, indices[i]);
		}
		std::cout << "print last 3 sorted features" << std::endl;
		for (size_t i = indices.size()-4; i <= indices.size() - 1; i++){
			Feature::printFeature(features, indices[i]);
		}
	}
}

void Matching :: checkEntropy(int column)
{
	if (features.type() != CV_32F){
		throw std::runtime_error("Mathing::checkEntropy(): unsupported feature-type");
	}
	for (int y = 0; y < features.rows; y++){
		if (features.at<float>(y, column) < config.entropyTh){
			mask(y, 0) = 0;
		}
	}
}

void Matching :: computeMatching(void)
{
	if ( Execution::verbosity >=1 )
		std::cout << "Matching::computeMatching" << std::endl;

	// (for e.g.Bravo-Solorio et al.) we need sort out the ones with too low entropy
	if (config.entropyTh > 0 && (crit & BRAVOCRIT)){ 
		if ( Execution::verbosity >=1 )
			std::cout << "erase blocks with low entropy\n";
		checkEntropy(3);
	}

	if ( Execution::verbosity >=1 )
		std::cout << "\tsort # " << cv::sum(mask)[0] << " " << std::endl;

	// use a kd tree
	if (config.kdsort){
		computeMatchingKD();
	}
	else if (config.nnBrute){
		computeNN();
	}
	else {
		// lexicographic sorting
		computeMatchingLexi();
	}	
}

