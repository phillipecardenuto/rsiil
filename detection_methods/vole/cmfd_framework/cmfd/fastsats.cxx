/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

/*
   Copyright(c) 2010 Vincent Christlein <vincentchristlein@web.de>
   and Christian Riesse <christian.riess@cs.fau.de>

   This file may be licensed under the terms of of the GNU General Public
   License, version 3, as published by the Free Software Foundation. You can
   find it here: http://www.gnu.org/licenses/gpl.html

   If you use this code in your research, please cite:
   V.Christlein and C.Riess and E.Angelopoulou."On Rotation Invariance in Copy-Move Forgery Detection",
   Workshop on Information Forensics and Security, Seattle, December 2010

   Suggestions for improvements or other feedback are very welcome!
  */

#include "fastsats.h"
#include "verification.h"
#include <stdexcept>
#include <iostream>
#include <algorithm>

Fastsats :: Fastsats(std::vector<cv::Mat_<cv::Point> > &corres, 
					 int _block_size,
					 int _verbosity,
					 int _max_dist,
					 int _step,
					 bool _color,
					 bool _with_tree)
	: corres_map(corres)
{
	if (corres.empty() || corres[0].empty())
		throw std::runtime_error("Fastsats::Fastsats: correspondence map == empty");
	if (_block_size < 1)
		throw std::runtime_error("Fastsats::Fastsats: block size must be at least 1");

	block_size = _block_size;
	verbosity = _verbosity;
	max_dist = _max_dist;
	step = _step;
	color = _color;
	with_tree = _with_tree;

	dimx = corres[0].cols;
	dimy = corres[0].rows;
	output = cv::Mat_<uchar>::zeros(dimy, dimx);
	data_array = NULL;
}

Fastsats :: ~Fastsats(void)
{
	if (data_array != NULL){
		delete[] data_array;
	}
}

void Fastsats :: buildTree()
{
	// TODO maybe build up one tree for every channel of corres_map
	int N = 0;
	for (int y = 0; y < dimy; y++){
		for (int x = 0; x < dimx; x++){
			// just check first plane
			if (corres_map[0](y,x).x != -1){
				N++;
			}
		}
	}
	if (N == 0){
		std::cerr << "WARNING: no correspondences -> can't build a tree -> return\n";
		return;
	}

	// let's use the 'new'-mechanisme and hope that this gives us continuous memory
	try {
		data_array = new float[N*2];
	}
	// just give some more output if bad_alloc exception happens,
	// but throw it further
	catch(std::bad_alloc ex){
		std::cerr << "Fastsats: bad alloc error\n";
		throw std::bad_alloc(ex);
	}

	int n1 = 0;
	for (int y = 0; y < dimy; y++){
		for (int x = 0; x < dimx; x++){
			// insert only from the first plane of corres_map
			if (corres_map[0](y,x).x != -1){
				data_array[n1++] = static_cast<float>(x);
				data_array[n1++] = static_cast<float>(y);
			}
		}
	}

	data = cv::Mat_<float>(N, 2, data_array);
	if (data.empty()){
		throw std::runtime_error("Fastsats::sameAffineTrafo: data-matrix is empty");
	}

	if (verbosity >= 1){
		std::cout << "create tree with " << N << " points\n";
	}
	//tree = new cv::flann::Index(data, cv::flann::KDTreeIndexParams());
	tree = new cv::flann::Index(data, cv::flann::LinearIndexParams());

	// parameter "satsSearchRadius"
	int num_nn = 25; // TODO think over that parameter

	indices = cv::Mat_<int>(N, num_nn, -1);
	dists = cv::Mat_<float> (N, num_nn);
	if (verbosity >= 1){
		std::cout << "knnsearch for " << num_nn << " neighbors\n";
	}
	// TODO: maybe radius search is more intelligent when using opencv v2.3
	//       but this is still kinda broken - can't search for many features at the same time :/
	tree->knnSearch(data, indices, dists, num_nn, cv::flann::SearchParams(64) );
}


/*
void Fastsats :: localRansac(int minSameTransformPairsTh, int recomputationTh)
{
	if (corres_map.empty())
		return;

	std::vector<cv::Point2f> points;
	std::vector<cv::Point2f> points_c;
	for (int y = 0; y < corres_map[0].rows; y++) {
		for (int x = 0; x < corres_map[0].cols; x++) {
			for (size_t i = 0; i < corres_map.size(); i++){
				if (corres_matrix[i](y,x).x != -1){
					cv::Point c = corres_map[i](y,x);
					points.push_back( cv::Point2f(x,y) );
					points_c.push_back(c);
				}
				break;
			}
		}
	}

	cv::Ptr<cv::DescriptorMatcher> ind = cv::DescriptorMatcher::create("BruteForce");
	std::vector<std::vector<cv::DMatch> > matches;
	ind->knnMatch(cv::Mat(points), cv::Mat(points), matches, 21);

	cv::Mat_<uchar> output_corres(output.rows, output.cols, (uchar)0);
	for (size_t i = 0; i < matches.size(); i++){
		std::vector<cv::Point2f> from;
		for (size_t k = 0; k < matches[i].size(); k++){
			cv::DMatch dm = matches[i][k];
			from.push_back(points[dm.queryIdx]);
		}
		cv::Point2f c = points_c[i];
		// welchen index hat c in matches??
		for (size_t p = 0; p < points.size(); p++){
			if (c != points[p])
				continue;
			std::vector<cv::Point2f> to;
			for (size_t k = 0; k < matches[p].size(); k++){
				cv::DMatch dm = matches[p][k];
				to.push_back(points[dm.queryIdx]);
			}
			to.push_back(points_c[dm.queryIdx]);

			std::vector<uchar> inliers;
			int num_inliers = Verification::ransac("ownransac", from, to, normalize, trafo, 100, 0.05, inliers);
			if (num_inliers < 4)
				continue;
			std::vector<cv::Point2f> hyp, hyp_c
			for (size_t = k; k < inliers.size(); k++){
				hyp.push_back(from[k]);
				hyp_c.push_back(to[]);
			}


			// do we already have marked it due to the region growing
		// process? (need to correct it with the half of the blockSize
		// cause we want to mark the midpoint of the block)
		if (output(y + static_cast<int>(block_size/2),
				   x + static_cast<int>(block_size/2)) > 0){
			continue;
		}
		}
	}
}
*/


void Fastsats :: satsWithTree(int minSameTransformPairsTh, int recomputationTh)
{
	buildTree();
	cv::Mat_<uchar> output_corres(output.rows, output.cols, (uchar)0);
	for (int i = 0; i < data.rows; i++){
		int x = (int)data.at<float>(i, 0);
		int y = (int)data.at<float>(i, 1);
		// do we already have marked it due to the region growing
		// process? (need to correct it with the half of the blockSize
		// cause we want to mark the midpoint of the blocks)
		// This is also true for keypoints as then we assume that blockSize == 1
		if (output(y + static_cast<int>(block_size/2),
				   x + static_cast<int>(block_size/2) ) > 0){
			continue;
		}
		cv::Point2f src[3];
		cv::Point2f dst[3];
		int src_ind[3];
		src_ind[0] = i;

		bool found_three = false;
		// try to get three points
		for (size_t c = 0; c < corres_map.size(); c++){
			// Note: need Point2f (==float) instead of Point (==Point2i) cause of getAffineTrafo()
			cv::Point2f corres;
			if (!getCorres(x, y, c, corres)){
				break;
			}
			cv::Point2f base_point(static_cast<float>(x),
								   static_cast<float>(y));

			src[0] = base_point;
			dst[0] = corres;

			// get three nearest neighbors and its correspondences
			int found_count = 1;
			for (int k = 0; k < indices.cols; k++){
				int index = indices(i,k);
				// same index as actual point or no index at all
				if (index == i || index == -1) {
					continue;
				}
				// if distance to its nearest neighbor too high
				// -> break as they are distance sorted
				if (dists(i,k) > max_dist*max_dist){
					break;
				}

				cv::Point2f source( data(index, 0), data(index, 1) );
				// have already marked it
//				if (output(source.y + static_cast<int>(block_size/2),
//						   source.x + static_cast<int>(block_size/2) ) > 0){
//					continue;
//				}
				bool found = false;
				cv::Point2f test_corres;
				for (size_t c = 0; c < corres_map.size(); c++){
					// we can break here (see findThreePoints())
					if (!getCorres((int)source.x, (int)source.y, c, test_corres)){
						break;
					}
					// the correspondence point may lie too far away
					// allow max twice the distance between the points
					if (wrongDist(test_corres, dst[0], std::sqrt(dists(i,k))*2) ){
						continue;
					}
					// check that the points found are not colinear
					if (found_count == 2 && clockwise(src[0], src[1], source) == 0){
						continue;
					}
					found = true;
					break;
				}
				if (found == true){
					src[found_count] = source;
					dst[found_count] = test_corres;
					src_ind[found_count] = index;
					found_count++;
					// have found three points -> break
					if (found_count == 3){
						break;
					}
				}
			}
			// TODO : these continues/breaks are not really clear programmed
			//        --> use functions!
			// haven't found 3 points
			if (found_count != 3){
				continue;
			}
			else {
				found_three = true;
				break;
			}
		}
		if (!found_three){
			continue;
		}

		// TODO: maybe add here also a bounding box check

		// three points form a hypothesis on the rotation and scaling
		// of a copied area.
		// Now try to extend the matching area
		// solve system of equations and get transformation matrix
		cv::Mat trafo;
		// returns 2x3 matrix transformation matrix (see Page 653 in opencv.pdf)
		trafo = cv::getAffineTransform(src,dst);

		cv::Point2f boundary[3];
		cv::Point2f boundary_corres[3];
		for (int j = 0; j < 3; ++j){
			boundary[j] = src[0];
			boundary_corres[j] = dst[0];
		}
		if ((src[1].x > src[2].x) || (src[1].y < src[2].y)) {
			boundary[1] = src[2];
			boundary[2] = src[1];
			boundary_corres[1] = dst[2];
			boundary_corres[2] = dst[1];
		} else {
			boundary[1] = src[1];
			boundary[2] = src[2];
			boundary_corres[1] = dst[1];
			boundary_corres[2] = dst[2];
		}

		cv::Mat_<uchar> visited = cv::Mat_<uchar>(dimy, dimx, static_cast<uchar>(0));
		// mark the src points as visited
		for (int j = 0; j < 3; j++){
			visited[static_cast<int>(src[j].y)][static_cast<int>(src[j].x)] = 255;
		}

		int myRecomputationTh = recomputationTh;

		// put the src points in our todo queue and our hypothesis points vector
		std::vector<cv::Point2f> hypothesis_points;
		hypothesis_points.push_back(src[0]);
		hypothesis_points.push_back(src[1]);
		hypothesis_points.push_back(src[2]);
		// todo-queue -> breadth-first search
		std::queue<int> todo;
		addNeighbors(src_ind[0], todo);
		addNeighbors(src_ind[1], todo);
		addNeighbors(src_ind[2], todo);
		// add the correspondences in the correspondence vector
		std::vector<cv::Point2f> hypothesis_points_corres;
		hypothesis_points_corres.push_back(dst[0]);
		hypothesis_points_corres.push_back(dst[1]);
		hypothesis_points_corres.push_back(dst[2]);

		while(!todo.empty()) {
			// perform breadth-depth search
			// get the top element from the todo-queue
			int curr_index = todo.front();
			cv::Point2f hyp(data(curr_index,0), data(curr_index,1));
			// remove it
			todo.pop();
			int x2 = static_cast<int>(hyp.x);
			int y2 = static_cast<int>(hyp.y);

			// already visited? ||  have already marked it
			if (visited[y2][x2] == 255
					|| output(y2 + static_cast<int>(block_size/2),
							  x2 + static_cast<int>(block_size/2) ) > 0){
				continue;
			}
			// if not already visited, visit it
			visited[y2][x2] = 255;

			// get correspondence of the point
			bool found_hypo = false;
			cv::Point2f hyp_corr;
			for (size_t c = 0; c < corres_map.size(); c++){
				// we can break here (s.below)
				if (!getCorres(x2, y2, c, hyp_corr)){
					break;
				}

				// check if this point could fit the hypothesis
				// compute distance between real match location
				// and expected location
				cv::Point2f expected = transformPoint(hyp, trafo);
				double max_estimation_error = 3.0;
				// TODO: the more points in hypothesis the
				// stable the assumption -> adjust smaller distance
				if (wrongDist(hyp_corr, expected, max_estimation_error)){
					continue;
				}
				found_hypo = true;
				// it fits - save correspondence
				hypothesis_points_corres.push_back(hyp_corr);
				break;
			} // end for all correspondences of point x2,y2
			if (!found_hypo){
				continue;
			}

			// if it fits, add its (not-yet-visited)
			// neighbors to the todo-list
			hypothesis_points.push_back(hyp);
			updateExtremalPoints(boundary, boundary_corres, hyp, hyp_corr);
			// add its (not-yet-visited) neighbors to the todo-list
			addNeighbors(curr_index, todo);

			// TODO: solve that better

			if (hypothesis_points.size() % myRecomputationTh == 0) {
				// We assume: the more points are in the
				// hypothesis set the more stable it is ->
				// change the intermediate myRecomputationTh
				myRecomputationTh *= recomputationTh;
				// trafo is recomputed,
				// IDEA: points are removed from hypothesis
				// points if they suddenly don't match
				// anymore
				trafo = recomputeTransformation(boundary, boundary_corres, hypothesis_points);
				if (hypothesis_points.size() == 4){
					Verification::affinity( cv::Mat(hypothesis_points),
											cv::Mat(hypothesis_points_corres),
											trafo );
				}
				else{ // at least 4
//					std::vector<uchar> inliers;
//					cv::Mat tmp_trafo;
//					int num_inliers = Verification::ransac( "ownransac",
//									 cv::Mat(hypothesis_points),
//									 cv::Mat(hypothesis_points_corres),
//									 true, tmp_trafo, 200, 0.05, inliers );
//					std::vector<cv::Point2f> tmp(num_inliers);
//					std::vector<cv::Point2f> tmp_c(num_inliers);
//					if ( num_inliers > 3 ){
//						trafo = tmp_trafo;
//						int cnt = 0;
//						for (size_t in = 0; in < inliers.size(); in++){
//							if (inliers[in] > 0){
//								tmp[cnt] = hypothesis_points[in];
//								tmp_c[cnt] = hypothesis_points_corres[in];
//								cnt++;
//							}
//						}
//						hypothesis_points = tmp;
//						hypothesis_points_corres = tmp_c;
//					}
				}
			}
		} // end todo

		// copy hypothesis
/*		Don't allow double elements
		std::vector<cv::Point2f> hyp_copy = hypothesis_points;
		std::vector<cv::Point2f> hyp_corres_copy = hypothesis_points_corres;
		// sort them
		std::sort(hyp_copy.begin(), hyp_copy.end(), small);
		std::sort(hyp_corres_copy.begin(), hyp_corres_copy.end(), small);
		// remove unique elements
		hyp_copy.erase(std::unique(hyp_copy.begin(), hyp_copy.end(), same), hyp_copy.end());
		hyp_corres_copy.erase(std::unique(hyp_corres_copy.begin(), hyp_corres_copy.end(), same), hyp_corres_copy.end());
		std::cout << "old sizes: " << hypothesis_points.size() << " , " << hypothesis_points_corres.size() << std::endl;
		std::cout << "new sizes: " << hyp_copy.size() << " , " << hyp_corres_copy.size() << std::endl;

		// apply thresholding:
		// if amount of tranformed points found > threshold -> assume copied region
		if (static_cast<int>(hyp_copy.size()) >= minSameTransformPairsTh
				&& static_cast<int>(hyp_corres_copy.size()) >= minSameTransformPairsTh)
*/
		if (static_cast<int>(hypothesis_points.size()) >= minSameTransformPairsTh)
		{
			if (verbosity >= 3) {
				std::cout << "Fastsats: found " << hypothesis_points.size() << " points\n";
			}
			std::vector<std::pair<cv::Point,cv::Point> > cluster;
			// mark all hypothesis points in a matrix
			for (unsigned int i = 0; i < hypothesis_points.size(); i++){
				cv::Point p(static_cast<int>(hypothesis_points[i].x),
							static_cast<int>(hypothesis_points[i].y));
				cv::Point p_c(static_cast<int>(hypothesis_points_corres[i].x),
							static_cast<int>(hypothesis_points_corres[i].y));
				int shift = static_cast<int>(block_size/2);
				// with a shift of half the block size -> so that they
				// mark the midpoints of the blocks
				output.at<uchar>(p.y + shift, p.x + shift)= 255;
				//output.at<uchar>(p_c.y + shift, p_c.x + shift)= 255;
				// more accurate but minimal slower and generates 2 clusters
				// for every direction one
				output_corres(p_c.y + shift, p_c.x + shift)= 255;
				// save cluster
				cluster.push_back(std::make_pair(p, p_c));
			}
			// store transformation matrix
			trans.push_back(trafo);
			// store the hypothesis points and its correspondences
			all_hypo_points.push_back(cluster);
		} // end if
	}
	output += output_corres;
}

void Fastsats :: sats(int minSameTransformPairsTh, int recomputationTh)
{
	cv::Mat_<uchar> output_corres = cv::Mat_<uchar>::zeros(output.rows, output.cols);
	for (int y = 0; y < (dimy - static_cast<int>(block_size/2)); y+=step) {
		for (int x = 0; x < (dimx - static_cast<int>(block_size/2)); x+=step) {
			// do we already have marked it due to the region growing
			// process? (need to correct it with the half of the blockSize
			// cause we want to mark the midpoint of the block)
			if (output(y + static_cast<int>(block_size/2),
					   x + static_cast<int>(block_size/2)) > 0){
				continue;
			}

			cv::Point2f src[3];
			cv::Point2f dst[3];
			// try to find 3 suitable points
			// store them in src and their correspondences in dst
			if (!findThreePoints(x, y, src, dst)){
				continue;
			}

			// three points form a hypothesis on the rotation and scaling
			// of a copied area.
			// Now try to extend the matching area
			// solve system of equations and get transformation matrix
			cv::Mat trafo;
			// returns 2x3 matrix transformation matrix opencv-doku
			trafo = cv::getAffineTransform(src, dst);

			// set initial boundary
			cv::Point2f boundary[3];
			cv::Point2f boundary_corres[3];
			boundary[0] = src[0]; // left
			boundary_corres[0] = dst[0];
			boundary[1] = src[2]; // lower
			boundary_corres[1] = dst[2];
			boundary[2] = src[1]; // right
			boundary_corres[2] = dst[1];

			cv::Mat_<uchar> visited = cv::Mat_<uchar>::zeros(dimy, dimx);
			// mark the src points as visited
			for (int j = 0; j < 3; j++){
				visited[static_cast<int>(src[j].y)][static_cast<int>(src[j].x)] = 255;
			}

			int myRecomputationTh = recomputationTh;

			// put the src points in our todo vector and our hypothesis points vector
			std::vector<cv::Point2f> hypothesis_points;
			// add the correspondences in the correspondence vector
			std::vector<cv::Point2f> hypothesis_points_corres;
			// todo-queue -> breadth-first search
			std::queue<cv::Point2f> todo;
			for (int i = 0; i < 3; i++){
				hypothesis_points.push_back(src[i]);
				hypothesis_points_corres.push_back(dst[i]);
				addNeighbors(src[i], todo);
			}

			while(!todo.empty()) {
				// perform breadth-depth search
				// get the top element from the todo-queue
				cv::Point2f hyp = todo.front();
				// remove it
				todo.pop();
				int x2 = static_cast<int>(hyp.x);
				int y2 = static_cast<int>(hyp.y);

				// already visited? || have already marked it
				if (visited[y2][x2] == 255//){
						|| output(y2 + static_cast<int>(block_size/2),
								  x2 + static_cast<int>(block_size/2) ) > 0){
					continue;
				}

				// if not already visited, visit it
				visited[y2][x2] = 255;

				// get correspondence of the point
				bool found_hypo = false;
				cv::Point2f hyp_corr;
				for (size_t c = 0; c < corres_map.size(); c++){
					// we can break here (s.above)
					if (!getCorres(x2, y2, c, hyp_corr)){
						break;
					}

					// check if this point could fit the hypothesis
					// compute distance between real match location
					// and expected location
					cv::Point2f expected = transformPoint(hyp, trafo);
					double max_estimation_error = 3.0;
                    // TODO: the more points in hypothesis the more
					// stable the assumption -> adjust smaller distance
					if (wrongDist(hyp_corr, expected, max_estimation_error)){
						continue;
					}
					found_hypo = true;
					// it fits - save correspondence
					hypothesis_points_corres.push_back(hyp_corr);
					break;
				} // end for all correspondences of point x2,y2
				if (!found_hypo){
					continue;
				}
				// if it fits we have a point more in our hypothesis-set
				hypothesis_points.push_back(hyp);
				// TODO: instead of using 3 extremal points to
				// compute the trafo-matrix compute the convex
				// hull-points and use them for a rigid trafo
				// estimation
				updateExtremalPoints(boundary, boundary_corres, hyp, hyp_corr);
				// add its neighbors to the todo-list anyway if they are already
				// visited or not
				addNeighbors(hyp, todo);

				// TODO: solve that better
				if (hypothesis_points.size() % myRecomputationTh == 0) {
					// We assume: the more points are in the
					// hypothesis set the more stable it is ->
					// change the intermediate myRecomputationTh
					myRecomputationTh *= recomputationTh;
					// trafo is recomputed,
					// IDEA: points are removed from hypothesis
					// points if they suddenly don't match
					// anymore
					trafo = recomputeTransformation(boundary, boundary_corres, hypothesis_points);
				}
			} // end while (!todo.empty())

			// apply thresholding:
			// if amount of tranformed points found > threshold -> assume copied region
			// mark in an output matrix
			if (static_cast<int>(hypothesis_points.size()) > minSameTransformPairsTh) {
				std::vector<std::pair<cv::Point,cv::Point> > cluster;
				// mark all hypothesis points in a matrix
				for (unsigned int i = 0; i < hypothesis_points.size(); i++){
					cv::Point p(static_cast<int>(hypothesis_points[i].x),
								static_cast<int>(hypothesis_points[i].y));
					cv::Point p_c(static_cast<int>(hypothesis_points_corres[i].x),
								static_cast<int>(hypothesis_points_corres[i].y));
					int shift = static_cast<int>(block_size/2);
					// with a shift of half the block size -> so that they
					// mark the midpoints of the blocks
					output.at<uchar>(p.y + shift, p.x + shift)= 255;
					//output.at<uchar>(p_c.y + shift, p_c.x + shift)= 255;
					// This is minimal slower but more accurate (than the
					// out-commented line above)
					// Note: it generates a cluster for every direction
					output_corres(p_c.y + shift, p_c.x + shift)= 255;
					// save cluster
					cluster.push_back(std::make_pair(p, p_c));
				}
				// store transformation matrix
				trans.push_back(trafo);
				// store the hypothesis points and its correspondences
				all_hypo_points.push_back(cluster);
			} // end if
		} // end for x
	} // end for y
	output += output_corres;
}

cv::Mat_<uchar> Fastsats :: sameAffineTransformationSelection(int minSameTransformPairsTh,
															  int recomputationTh)
{
	// let's measure the time
	double t = (double)cv::getTickCount();
	
	if (with_tree){
		satsWithTree(minSameTransformPairsTh, recomputationTh);
	}
	else {
		sats(minSameTransformPairsTh, recomputationTh);
	}

	if ( verbosity >= 2 ){
		t = ((double)cv::getTickCount() - t)/cv::getTickFrequency();
		std::cout << "time needed for same affine trafo selection: " << t << std::endl;
	}
	return output;
}

// add neighbors in a circle fashion way
void Fastsats :: addNeighbors(const cv::Point2f &p,
							  std::queue<cv::Point2f> &todo) const
{
	int px = (int)p.x;
	int py = (int)p.y;
	/*
	int minX = std::max<int>(0, px-2);
	int maxX = std::min<int>(dimx, px+2);
	int minY = std::max<int>(0, py-2);
	int maxY = std::min<int>(dimy, py+2);
	for (int y = minY; y < maxY; ++y) {
		for (int x = minX; x < maxX; ++x) {
			todo.push(cv::Point2f(x, y));
		}
	}
	*/

	int num_circles = 2;
	int sign = 1;
	for (int dirs = 1; dirs <= 2*num_circles+1; dirs++){
		for (int nx = 1; nx <= dirs; nx++){
			// end of circle
			if (dirs == 2*num_circles+1 && nx == dirs){
				return;
			}
			px += sign*step;
			if (px >= dimx || px < 0 || py >= dimy || py <0 ) // avoid segfault
				continue;
			todo.push(cv::Point2f(px, py));
		}
		for (int ny = 1; ny <= dirs; ny++){
			py += sign*step;
			if (py >= dimy || py < 0 || px >= dimx || px < 0) // avoid segfault
				continue;
			todo.push(cv::Point2f(px, py));
		}
		sign *= -1;
	}
}

void Fastsats :: addNeighbors(int point_index,
							  std::queue<int> &todo) const
{
	for (int k = 0; k < indices.cols; k++){
		int index = indices(point_index, k);
		 // same index as actual point or no index at all
		if (index == point_index || index == -1){
			continue;
		}
		// if distance to its nearest neighbor too high break
		// (as the neighbors are distance-sorted)
		if (dists(point_index, k) > max_dist*max_dist){
			break;
		}
		// add to todo-list
		todo.push(index);
	}
}

// tries to find 3 points which fullfill certain criteria:
// 1. have correspondences
// 2. which are close to each other
// 3. bounding box of these three points is similar to the one
//    of the bounding box of the correspondence points
bool Fastsats :: findThreePoints(int x, int y, cv::Point2f *src,
								 cv::Point2f *dst)
{
	for (size_t c = 0; c < corres_map.size(); c++){
		cv::Point2f corres;
		// 1. check point -> does it have a correspondence
		//    can break here as we assume that further correspondence maps
		//    don't have a correspondence here any more
		if ( ! getCorres(x, y, c, corres)){
			break; // no match? leave...
		}
		cv::Point2f base_point(x, y);

		src[0] = base_point;
		dst[0] = corres;

		// 2. try to find neighboring points which have correspondences
		//    lying not too far away from the current correspondence

		// test some possible pts to get the 2nd point
		// TODO maybe test here some more points than just 4
		std::vector<cv::Point2f> test_points;
		test_points.reserve(4);
		test_points.push_back(cv::Point2f(x+3, y));
		test_points.push_back(cv::Point2f(x+2, y));
		test_points.push_back(cv::Point2f(x+3, y+1));
		test_points.push_back(cv::Point2f(x+2, y+1));

		if (!testPoints(base_point, corres, test_points,
						src[1], dst[1])){
			continue;
		}

		// try to get a third point
		test_points[0] = cv::Point2f(x, y+3);
		test_points[1] = cv::Point2f(x, y+2);
		test_points[2] = cv::Point2f(x+1, y+3);
		test_points[3] = cv::Point2f(x+1, y+2);

		if (!testPoints(base_point, corres, test_points,
						src[2], dst[2])){
			continue;
		}

		// 3. now check bounding box
		int bbx = (int)src[1].x - x;
		int bby = (int)src[2].y - y;
		// area bounding box
		int bb_area = bbx*bby;

		// get bounding box
		int minx = dimx;
		int miny = dimy;
		int maxx = 0;
		int maxy = 0;
		for (int i = 0; i < 3; i++){
			if (dst[i].x < minx)
				minx = (int)dst[i].x;
			if (dst[i].x > maxx)
				maxx = (int)dst[i].x;
			if (dst[i].y < miny)
				miny = (int)dst[i].y;
			if (dst[i].y > maxy)
				maxy = (int)dst[i].y;
		}
		// area of bounding box
		int bb_area_corres = (maxx-minx)*(maxy-miny);
		// check if area of bounding box of corresponding points is too big
		// allow double size -> 4
		double max_bb_area_factor = 4;
		double factor;
		if (bb_area_corres > bb_area){
			factor = static_cast<double>(bb_area_corres) / static_cast<double>(bb_area);
		} else {
			factor = static_cast<double>(bb_area) / static_cast<double>(bb_area_corres);
		}
		if (factor > max_bb_area_factor){
			continue;
		}
		// all tests passed
		return true;
	}
	// end searching for three suitable points
	// -> nothing found
	return false;
}

// take the first point which fits all criteria
// and return its position in the vector
bool Fastsats :: testPoints(const cv::Point2f & base,
							const cv::Point2f & corres,
							const std::vector<cv::Point2f> & test_points,
							cv::Point2f & result_point,
							cv::Point2f & result_corres) const
{
	for (size_t i = 0; i < test_points.size(); i++){
		// prevent adding already (globally) visited points
		if (output.at<uchar>((int)test_points[i].y + static_cast<int>(block_size/2),
							 (int)test_points[i].x + static_cast<int>(block_size/2)) > 0){
			continue;
		}
		// compute distance between base and the testing point
		cv::Point2f tmp = base - test_points[i];
		double dist = tmp.ddot(tmp); // no need for sqrt
        // go through the correspondence map and chose corres-point
		// which fits first		
		for (size_t c = 0; c < corres_map.size(); c++){
			cv::Point2f test_corres;
            // 1. has this a correspondence point at all
			//    (if not, than the lower corres_maps neither have one
			//    thus we can break here)
			if ( ! getCorres(static_cast<int>(test_points[i].x),
							 static_cast<int>(test_points[i].y),
							 c, test_corres) ){
				break;
			}
			// 2. test if it's correspondence point not too far away
			//    from the corres-point of the base-point
			if (wrongDist(test_corres, corres, dist*2)){
				continue;
			}

			// we have passed both tests
			result_point = test_points[i];
			result_corres = test_corres;
			return true;
		}
	}

	return false;
}

bool Fastsats :: getCorres(int x, int y, int c, cv::Point2f & result) const
{
	cv::Point2f tmp(static_cast<float>(corres_map[c](y,x).x),
					static_cast<float>(corres_map[c](y,x).y));
	if (tmp.x == -1){
		return false;
	}
	result = tmp;
	return true;
}

/*
   clockwise gives an indicator for the relative position of 3 points:
   - the points are colinear if clockwise == 0
   - the points are in clockwise order (= describe a turn right) if clockwise > 0
   - the points are in counterclockwise order if clockwise < 0
   this is an extremely cheap operation that can be used as a sorting criterion for the graham scan algorithm
   */
int Fastsats :: clockwise(const cv::Point2f idx1, const cv::Point2f idx2, const cv::Point2f idx3)
{
	double tmp = (idx2.x - idx1.x) * (idx3.y - idx1.y) - (idx2.y - idx1.y) * (idx3.x - idx1.x);
	if (fabs(tmp) < DBL_EPSILON) return 0;
	if (tmp > 0) return 1;
	return -1;
}

bool Fastsats :: wrongDist(const cv::Point2f &p1, const cv::Point2f &p2, double maxDist) const
{
	cv::Point2f diff = p2 - p1;
	/*
	double dist = cv::norm(diff.x * diff.x + diff.y * diff.y);
	return (dist > maxDist);
	*/
	// should be faster than that above
	double dist = diff.ddot(diff);
	return (dist > maxDist*maxDist);
}

void Fastsats :: updateExtremalPoints(cv::Point2f *extremalPoints,
									  cv::Point2f *extremalPoints_corres,
									  cv::Point2f newPoint,
									  cv::Point2f newPoint_corres) const
{
	// left
	if (extremalPoints[0].x > newPoint.x
			&& clockwise(extremalPoints[1], extremalPoints[2], newPoint) != 0) {
		extremalPoints[0] = newPoint;
		extremalPoints_corres[0] = newPoint_corres;
		return;
	}
	// lower
	if (extremalPoints[1].y < newPoint.y
			&& clockwise(extremalPoints[0], extremalPoints[2], newPoint) != 0) {
		extremalPoints[1] = newPoint;
		extremalPoints_corres[1] = newPoint_corres;
		return;
	}
	// right
	if (extremalPoints[2].x < newPoint.x
			&& clockwise(extremalPoints[0], extremalPoints[1], newPoint) != 0) {
		extremalPoints[2] = newPoint;
		extremalPoints_corres[2] = newPoint_corres;
		return;
	}
}

cv::Mat Fastsats :: recomputeTransformation( const cv::Point2f *extremalPoints,
											 const cv::Point2f *extremalPoints_corres,
											 std::vector<cv::Point2f> &hypothesis_points) const
{
	// solve system of equations and get transformation matrix
	cv::Mat trafo;
	// returns 2x3 matrix transformation matrix (see opencv-doku)
	// TODO maybe use more points + svd to compute a better approximation
	trafo = cv::getAffineTransform(extremalPoints, extremalPoints_corres);

	// TODO: double-check hypothesis points if they still fit the model

	return trafo;
}

cv::Point2f Fastsats :: transformPoint( const cv::Point2f &point, const cv::Mat &transformation )
{
	cv::Mat_<cv::Vec2f> oldPoint(1,1);
	oldPoint(0,0)[0] = point.x;
	oldPoint(0,0)[1] = point.y;
	cv::Mat newPoint;

	if (transformation.rows == 2){
		cv::transform(oldPoint, newPoint, transformation);
	}
	else {
		cv::perspectiveTransform(oldPoint, newPoint, transformation);
	}
	cv::Point2f dstPoint(newPoint.at<cv::Vec3f>(0,0)[0], newPoint.at<cv::Vec3f>(0,0)[1]);
	return dstPoint;
}
