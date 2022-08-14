/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <iostream> 
#include <iomanip>
#include <queue>
#include <limits.h>
#include <stdexcept> 

#include <cv.h>
#include <highgui.h>

#include "execution.h"

#include "verification.h"
#include "fastsats.h"
// Does not compile yet :(
//#include "homest.h"

Verification :: Verification(int blockSize, 
							 int width,
							 int height,
							 std::vector<cv::Mat_<cv::Point> > & _corres,
							 int decompLvl )
	:	corres(_corres)
{
	block_size = blockSize;
	dimx = width;
	dimy= height;
	decomposeLvl = decompLvl;
	if (corres.empty()){
		throw std::runtime_error("Verification::Verification(): corres.empty() !?");
	}
}

Verification :: ~Verification(void)
{}

double Verification :: rmse(const cv::Mat & a, const cv::Mat & b)
{
	CV_Assert(a.size == b.size);
//	cv::Mat diff = cv::abs(a) - cv::abs(b);
	cv::Mat diff = a - b;
	cv::Mat s = diff.mul(diff);
	double sum = cv::sum(s)[0];
	sum /= (a.rows * a.cols);
	return std::sqrt(sum);
}

// TODO: can we also skip if scales are negative or is that normal?
bool Verification :: skip(const cv::Mat & _trafo) const
{
	if ( _trafo.empty() ) {
		vout(2) << "skip (trafo = empty)\n";
		return true;
	}
	cv::Mat_<double> trafo = _trafo;
	vout(2) << "->skip check: " << printTrafo(trafo);
	double sx = std::sqrt( trafo(0,0)*trafo(0,0) + trafo(0,1)*trafo(0,1) );
	if (trafo(0,0) < 0) sx *= -1;
	double sy = std::sqrt( trafo(1,0)*trafo(1,0) + trafo(1,1)*trafo(1,1) );
	if (trafo(1,1) < 0) sy *= -1;
	// ignore abstruse scales
	if ( ( (trafo.rows == 2)
		 || ( (trafo.rows == 3) && (trafo(2,0) == 0.0
								&& trafo(2,1) == 0.0 && trafo(2,2) == 1.0)) )
		 && (std::abs(sx) > 2.2 || std::abs(sx) < 0.4 || std::abs(sy) > 2.2 || std::abs(sy) < 0.4) )
	{
		vout(2) << "<-skip (scale is too abstruse)\n";
		return true;
	}
	return false;
}

std::string Verification :: printTrafo(const cv::Mat & _trafo)
{
	std::stringstream ss;
	cv::Mat_<double> trafo = _trafo;
	double sx = std::sqrt( trafo(0,0)*trafo(0,0) + trafo(0,1)*trafo(0,1) );
	//				if (trafo(0,0) < 0) sx *= -1;
	double sy = std::sqrt( trafo(1,0)*trafo(1,0) + trafo(1,1)*trafo(1,1) );
	//				if (trafo(1,1) < 0) sy *= -1;
	double angle = std::atan(trafo(1,0) / trafo(1,1));
	// radians -> degree, acos(-1.0) == pi
	angle *= 180 / std::acos(-1.0);
	// doesn't work, opencv doesn't allow to pipe this in a stringstream
	ss << "trafo: " << trafo << std::endl;
	ss << "measured scale: sx: " <<  sx << " sy: " << sy << std::endl;
	ss << "measured angle: " << angle << std::endl;
	return ss.str();
}

// try to merge similar transformation matrix
void Verification :: mergeTrafos(std::vector<cv::Mat> & all_trafos,
								 std::vector<std::vector<cv::Vec2f> > & all_inliers,
								 /* ransac parameters */
								 const std::string & estimateMethod,
								 bool normalize,
								 int num_iterations,
								 double reprojTh)
{
	vout(2) << "->merge " << all_trafos.size() << " trafos\n";
	// Note: probablly I could merge more but it's already quite good
	std::vector<bool> merged(all_trafos.size(), false);
	for (size_t i = 0; i < all_trafos.size(); i++){

		for (size_t k = 0; k < all_trafos.size(); k++){
			// continue if we are looking at the same trafo or we've already merged the trafo
			if ( i == k || merged[i] || merged[k] ) {
				continue;
			}

			cv::Mat_<double> trafo1 = all_trafos[i];
			cv::Mat_<double> trafo2 = all_trafos[k];
			double e1 = rmse( trafo1(cv::Range(0,2), cv::Range(0,2)),
							  trafo2(cv::Range(0,2), cv::Range(0,2)) );
			double e2 = rmse( trafo1, trafo2);
			//				vout(2) << "--> rmse: i and k:" << e1 << " and " << e2 << std::endl << std::endl;

			// merge if rmse < 0.03 and if shifts are going in the same direction
			if ( e1 < 0.03 && e2 < 10)
				//						&& ( (trafo1(0,2) > 0 && trafo2(0,2) > 0 )
				//							 || (trafo1(0,2) < 0 && trafo2(0,2) < 0 ) )
				//						&& ( (trafo1(1,2) > 0 && trafo2(1,2) > 0 )
				//							 || (trafo1(1,2) < 0 && trafo2(1,2) < 0 ) ) )
			{
				//					std::cout << " merge clusters\n";
				// check signs of the angles
				//					double angle1 = std::atan(trafo1(1,0) / trafo1(1,1));
				//					double angle2 = std::atan(trafo2(1,0) / trafo2(1,1));
				//					bool same_sign = false;
				//					if ( (angle1 < 0 && angle2 < 0) || (angle1 > 0 && angle2 > 0) )
				//						same_sign = true;

				std::vector<cv::Vec2f> from = all_inliers[i*2];
				std::vector<cv::Vec2f> to = all_inliers[i*2+1];
				// take care about the direction points->correspondences vs correspondences->points
				//					if (same_sign){
				//						std::cout << "same sign\n";
				from.insert( from.end(), all_inliers[k*2].begin(), all_inliers[k*2].end() );
				to.insert( to.end(), all_inliers[k*2+1].begin(), all_inliers[k*2+1].end() );
				//					}
				//					else {
				//						from.insert( from.end(), all_inliers[k*2+1].begin(), all_inliers[k*2+1].end() );
				//						to.insert( to.end(), all_inliers[k*2].begin(), all_inliers[k*2].end() );
				//					}
				cv::Mat merged_trafo;
				std::vector<uchar> inliers;
				ransac(estimateMethod, cv::Mat(from), cv::Mat(to), normalize,
					   merged_trafo, num_iterations,
					   reprojTh, inliers);

				merged[i] = true;
				merged[k] = true;
				all_trafos.push_back(merged_trafo);
				all_inliers.push_back(from);
				all_inliers.push_back(to);
				// we may want to merge that again with another one
				merged.push_back(false);
			}
		}
	}

	//	std::cout << " -------- FINAL TRAFOS ------------\n";
	std::vector<cv::Mat> new_trafos;
	for (size_t i = 0; i < all_trafos.size(); i++){
		if ( ! merged[i] ){			
			new_trafos.push_back( all_trafos[i] );

			cv::Mat_<double> tmp = all_trafos[i].clone();
		}
	}
	all_trafos = new_trafos;
	vout(2) << "<-merge, have now " << all_trafos.size() << " trafos\n";
	if ( Execution::verbosity >= 2 && all_trafos.size() < 5 ){
		std::cout << "--- List all trafos ---\n";
		for (size_t i = 0; i < all_trafos.size(); i++){
			std::cout << printTrafo(all_trafos[i]);
		}
	}

	// TODO: do the same for all_inliers, but this is atm not used further
}

std::vector<cv::Mat> Verification :: ransac( const std::string & estimateMethod,
											 std::vector<std::vector<std::pair<cv::Point,
												 cv::Point> > > clusters,
											 int num_iterations,
											 double reprojTh,
											 bool normalize)
{
	std::vector<cv::Mat> all_trafos;
	std::vector<std::vector<cv::Vec2f> > all_inliers;
	// compute ransac for every cluster combination
	for (size_t first = 0; first < clusters.size(); first++){
		std::vector<cv::Vec2f> from;
		std::vector<cv::Vec2f> to;
		//cv::Mat_<uchar> from_m(dimy, dimx, (uchar)0);
		//cv::Mat_<uchar> to_m(dimy, dimx, (uchar)0);
		for (size_t i = 0; i < clusters[first].size(); i++){
			cv::Vec2f shift( block_size/2.f, block_size/2.f );
			from.push_back( cv::Vec2f(clusters[first][i].first.x, clusters[first][i].first.y) + shift);
			to.push_back( cv::Vec2f(clusters[first][i].second.x, clusters[first][i].second.y) + shift);

			//from_m( clusters[first][i].first.y ,clusters[first][i].first.x ) = 255;
			//to_m( clusters[first][i].second.y ,clusters[first][i].second.x ) = 255;
		}

		// get transformation matrix through ransac and add it to trafos
		cv::Mat trafo;
		std::vector<uchar> inliers;
		ransac(estimateMethod, cv::Mat(from), cv::Mat(to), normalize,
			   trafo, num_iterations,
			   reprojTh, inliers);
		//vout(2) << "Trafo: " << trafo << std::endl;
		if ( skip(trafo) ) {
			continue;
		}
		all_trafos.push_back(trafo);

		// inliers of from and to vectors
		std::vector<cv::Vec2f> f_i;
		std::vector<cv::Vec2f> t_i;
		if (!inliers.empty()){
			for (size_t i = 0; i < inliers.size(); i++){
				if (inliers[i] > 0){
					f_i.push_back(from[i]);
					t_i.push_back(to[i]);
				}
			}
			// save them
			all_inliers.push_back(f_i);
			all_inliers.push_back(t_i);
		}
		else {
			all_inliers.push_back(from);
			all_inliers.push_back(to);
		}

	}

	Verification::mergeTrafos(all_trafos, all_inliers, estimateMethod,
							  normalize, num_iterations, reprojTh);

	return all_trafos;
}

// This computes one global transformation from the correspondence maps, i.e.
// for images with multiple copies this procedure has to be applied multiple times
// TODO: actually also points from to have to be in from and vice versa
//       or put both points during the matching process to the correspondence maps
cv::Mat Verification :: ransac(const std::string & estimateMethod,
							   int num_iterations,
							   double reprojTh,
							   bool normalize)
{
	std::vector<cv::Vec2f> from;
	std::vector<cv::Vec2f> to;
	// get all relevant points from every channel of the correlation matrix
	// thus there could be more points if e.g. numRows > 1 (matching stage)
	for (size_t c = 0; c < corres.size(); c++) {
		for (int y = 0; y < corres[c].rows; y++){
			for (int x = 0; x < corres[c].cols; x++){
				if (corres[c](y,x).x != -1){
					from.push_back( cv::Vec2f(x,y) );
					to.push_back( cv::Vec2f( corres[c](y,x).x, corres[c](y,x).y) );
				}
			}
		}
	}
	cv::Mat trafo;
	std::vector<uchar> inliers;
	ransac(estimateMethod, cv::Mat(from), cv::Mat(to), normalize, trafo, num_iterations, reprojTh, inliers);
	if ( skip(trafo) ) {
		return cv::Mat();
	}
	return trafo;
}

// normalize points in src as suggested by Amerini (which refers to Hartley)
// i.e. the data gets zero-mean and the avg-distance to the center is sqrt(2)
void Verification :: normalizeCluster(cv::Mat & m,
									  cv::Mat & t)
{
	CV_Assert( m.cols == 1 );
	CV_Assert( m.type() == CV_32FC2 );

	// centroid:
	cv::Scalar mean = cv::mean(m);
	// translate
	//for (size_t i = 0; i < vec.size(); i++){
	for (int i = 0; i < m.rows; i++){
		//vec[i] -= cv::Vec2f(mean[0], mean[1]);
		m.at<cv::Vec2f>(i,0) -= cv::Vec2f(mean[0], mean[1]);
	}
	// compute avg-distance
	double avg_dist = 0.0;
	for (int i = 0; i < m.rows; i++){
	//for (size_t i = 0; i < vec.size(); i++){
		//avg_dist += cv::norm(vec[i]);
		avg_dist += cv::norm(m.at<cv::Vec2f>(i,0));
	}

	avg_dist /=  m.rows;
	// scale factor
	cv::Vec2f scale_factor( sqrt(2.0),sqrt(2.0) );
	scale_factor /= avg_dist;

	for (int i = 0; i < m.rows; i++){
		m.at<cv::Vec2f>(i,0) *= scale_factor[0];
	}

	t = (cv::Mat_<double>(3,3) <<
		 1. / avg_dist, 0., -mean[0] / avg_dist,
		 0., 1. / avg_dist, -mean[1] / avg_dist,
		 0., 0., 1.);
}

void Verification :: denormalizeH( const cv::Mat & t1, const cv::Mat & t2, cv::Mat & trafo )
{
	// Note: convertsion to double needed
	if (trafo.rows == 2){
		cv::Mat_<double> t = trafo;
		cv::Mat_<double> tmp(3,3);

		cv::Mat up(tmp(cv::Range(0,2),cv::Range(0,3)));
		t.copyTo(up);

		tmp(2,0) = tmp(2,1) = 0.0;
		tmp(2,2) = 1.0;
		trafo = t2.inv() * tmp * t1;
	}
	else {
		cv::Mat_<double> tmp = trafo;
		trafo = t2.inv() * tmp * t1;
	}

}

// computes the 3x3 affine homogenous transformation matrix
void Verification :: affinity(const cv::Mat & from,
							  const cv::Mat & to,
							  cv::Mat & trafo)
{
	CV_Assert(from.type() == CV_32FC2
			  && from.type() == to.type()
			  && from.cols == 1);
	// following Hartley, p.130
	cv::Mat f = from.clone();
	cv::Mat t = to.clone();

	// reshape to one channel: nr-points x 2 matrix (have to assign it to new matrices)
	cv::Mat_<float> f_m = f.reshape( 1, from.rows );
	cv::Mat_<float> t_m = t.reshape( 1, from.rows );

	// make data zero mean:
	// get centroid = mean
	cv::Scalar f_mean = cv::mean(f);
	cv::Scalar t_mean = cv::mean(t);
	// convert to matrices
	cv::Mat tA = (cv::Mat_<float>(1,2) << f_mean[0], f_mean[1]);
	cv::Mat tB = (cv::Mat_<float>(1,2) << t_mean[0], t_mean[1]);
	// subtract from data -> zero-mean
	cv::Mat_<float> one = cv::Mat_<float>::ones(f_m.rows, 1);
	f_m -= one * tA;
	t_m -= one * tB;

	cv::Mat_<float> data(from.rows, 4);
	cv::Mat tmp1(data(cv::Range(0,data.rows), cv::Range(0,2)));
	cv::Mat tmp2(data(cv::Range(0,data.rows), cv::Range(2,4)));
	f_m.copyTo(tmp1);
	t_m.copyTo(tmp2);
	// now data contains f_m and t_m

	// decompose data
	cv::SVD s(data, cv::SVD::FULL_UV);
	cv::Mat_<float> v = s.vt.t();

	// get the right singular-vectors of A corresponding to the two largest
	// singular values
	cv::Mat_<float> M = v(cv::Range(0,2),cv::Range(0,2));
	cv::Mat_<float> N = v(cv::Range(2,4), cv::Range(0,2));

	cv::Mat_<float> Hs = N * M.inv();
	cv::Mat_<float> t1 = Hs * tA.t();
	cv::Mat_<float> tmp = -(t1 - tB.t());

	cv::Mat H = (cv::Mat_<float>(3,3) << Hs(0,0), Hs(0,1), tmp(0,0),
				 Hs(1,0), Hs(1,1), tmp(1,0),
				 0.f, 0.f, 1.f);
	//	cv::Mat H = (cv::Mat_<float>(2,3) << Hs(0,0), Hs(0,1), tmp(0,0),
	//				 Hs(1,0), Hs(1,1), tmp(1,0));
	trafo = H;
}

int Verification::ransac(const std::string & estimateMethod,
						 const cv::Mat & _from,
						 const cv::Mat & _to,
						 bool normalize,
						 cv::Mat & trafo,
						 int num_iterations,
						 double reprojTh,
						 std::vector<uchar> & inliers)
{	
	if (estimateMethod != "ownransac"
			&& estimateMethod != "ransac"
			&& estimateMethod != "lmeds"
			&& estimateMethod != "pure"
			&& estimateMethod != "goldstandard"
			/*&& estimateMethod != "homest"*/)
	{
		throw std::runtime_error("Verification::ransac(): wrong estimateMethod given!");
	}
	CV_Assert( _from.type() == CV_32FC2
			   && _from.type() == _to.type()
			   && _from.size == _to.size
			   && _from.cols == 1);

	// make local copies which may be changed
	// e.g. due to normalization
	cv::Mat from = _from.clone();
	cv::Mat to = _to.clone();
//	cv::Mat from = from_tmp;
//	cv::Mat to = to_tmp;
//	// guarantee one col
//	if ( from_tmp.cols != 1 ){
//		from = from_tmp.reshape( 2, from_tmp.cols*from_tmp.rows );
//		to = to_tmp.reshape( 2, to_tmp.cols*to_tmp.rows );
//	}

	cv::Mat t1, t2;
	if ( normalize ){
		vout(2) << " normalize clusters\n";
		normalizeCluster(from,t1);
		normalizeCluster(to,t2);
	}

	// OpenCV-RANSAC
	// Gives homography	matrix
	int num_inliers = 0;
	if (estimateMethod == "ransac" || estimateMethod == "lmeds" || estimateMethod == "pure"){
		if (Execution::verbosity >= 1){
			std::cout << "-> Compute transformation matrix using RANSAC (OpenCV-version)"
					  << "-> gives perspective trafo\n";
		}

		// # of iterations of the Ransac algo can't be chosen here
		// thus it differs from the version of Pan/Amerini
		if ( estimateMethod == "pure" ){
			trafo = cv::findHomography(from, to);
		}
		if ( estimateMethod == "lmeds" ){
			trafo = cv::findHomography(from, to, CV_LMEDS);
		}
		else {
			trafo = cv::findHomography(from, to, CV_RANSAC, reprojTh, inliers);
		}
		for (size_t k = 0; k < inliers.size(); k++){
			if (inliers[k] > 0) {
				num_inliers++;
			}
		}
		if (Execution::verbosity >= 1){
			std::cout << "<- Compute transformation matrix using RANSAC, got "
					  << num_inliers << " inliers\n\n";
		}
	}

	/*
	else if (estimateMethod == "homest")
	{
		int n_points = from.size();
		// kinda uggly conversion
		double **from_a = new double*[n_points];
		double **to_a = new double*[n_points];
		for (size_t i = 0; i < from.size(); i++){
			from_a[i] = new double[2];
			from_a[i][0] = static_cast<double>(from[i][0]);
			from_a[i][1] = static_cast<double>(from[i][1]);
			to_a[i] = new double[2];
			to_a[i][0] = static_cast<double>(to[i][0]);
			to_a[i][1] = static_cast<double>(to[i][1]);
		}

		double H[9];
		homestaff((double(*)[2])from_a,(double(*)[2]) to_a, n_points, 0.9, H, true, NULL, NULL, Execution::verbosity);
		trafo = (cv::Mat_<float>(2,3) << H[0], H[1], H[2], H[3], H[4], H[5]);
		vout(2) << "test if last row is 0 0 1: " << H[6] << " " << H[7] << " " << H[8] << std::endl;
		// clenup:

		for (size_t i = 0; i < from.size(); i++){
			delete[] from_a[i];
			delete[] to_a[i];
		}

		delete[] from_a;
		delete[] to_a;
	}
	*/

	// --- own RANSAC ---
	else if (estimateMethod == "ownransac")
	{
		if (Execution::verbosity >= 1){
			std::cout << "-> Compute transformation matrix of " << from.rows
					  << " pairs using RANSAC (own version)"
					  << " with " << num_iterations << " iterations\n";
		}
		// initialization of parameters
		//srand(time(NULL));
		srand(0);
		int num_pairs = from.rows;
		// parameters which will be estimated/set
		cv::Mat final_trafo;
		int max_inliers = 0;
		// 1:hit ; 0:miss
		std::vector<uchar> tmp_inliers(num_pairs);

		// For the number of 'num_iterations' chose the transformation
		// which gives the most inliers
		for (int i = 0; i < num_iterations; i++)
		{
			// at least three points for affine trafo (at least four for perspective)
			// Note: 3 seems to return better results
			int num_points = 3;

			cv::Point2f *pts_from = new cv::Point2f[num_points];
			cv::Point2f *pts_to = new cv::Point2f[num_points];

			// chose 'num_points' random pairs which are not colinear
			// try this also only maximal num_iterations
			bool found = false;
			for (int t = 0; t < num_iterations; t++)
			{
				// chose 'numPoints' random pairs
				for (int k = 0; k < num_points; k++){
					int rnd_pos = rand() % from.rows;
					pts_from[k] = from.at<cv::Point2f>(rnd_pos,0);
					pts_to[k] = to.at<cv::Point2f>(rnd_pos, 0);
				}
				// if we use 3 pairs for the estimate of the affine transformation
				// they shouldn't be colinear
				if ( num_points == 3 ){
					// check that the points are not colinear
					if ( (Fastsats::clockwise(pts_from[0], pts_from[1], pts_from[2]) != 0 )
						 && (Fastsats::clockwise(pts_to[0], pts_to[1], pts_to[2]) != 0) ){
						found = true;
						break;
					}
				}
				else {
					found = true;
					break;
				}
			}
			// this may happen if too few correspondences were found
			// and if they are also colinear
			if ( !found ){
				vout(2) << "\tcouldn't find a transformation\n";
				trafo = cv::Mat();
				// cleanup
				delete[] pts_from;
				delete[] pts_to;
				return 0;
			}

			// get the transformation matrix
			// change that to getPerspectiveTransform if 4 points
			cv::Mat tmp_trafo;
			if (num_points == 3){
				// affine transformation -> 2x3 matrix
				tmp_trafo = cv::getAffineTransform(pts_from, pts_to);
			}
			else if (num_points == 4){
				tmp_trafo = cv::getPerspectiveTransform(pts_from, pts_to);
			}

			int tmp_num_inliers = 0;
			// check for every pair if it is an inlier or outlier
			for (int k = 0; k < num_pairs; k++)
			{
				if (tmp_trafo.empty()) {
					break;
				}
				// least square error
				cv::Point2f tmp = to.at<cv::Point2f>(k,0)
						- Fastsats::transformPoint( from.at<cv::Point2f>(k,0), tmp_trafo);
				double error = std::sqrt(tmp.dot(tmp)); // norm of the point
				// if backprojection-error < threshold -> inlier
				if (error <= reprojTh){
					tmp_inliers[k] = 1;
					tmp_num_inliers++;
				}
				else {// outlier
					tmp_inliers[k] = 0;
				}
			}

			// if we get a better result: update
			if (max_inliers < tmp_num_inliers){
				final_trafo = tmp_trafo;
				max_inliers = tmp_num_inliers;
				inliers = tmp_inliers;
			}

			// cleanup
			delete[] pts_from;
			delete[] pts_to;
			// robust enough -> break
			// TODO: maybe make the 0.97 also as input variable
			if ( max_inliers > 0.97*from.rows ){
				if (Execution::verbosity >= 2){
					std::cout << " break earlier cause : " << max_inliers << " > "
							  << 0.97*from.rows << std::endl;
				}
				break;
			}
		} // end of ransac-iterations

		if (max_inliers == 0 || final_trafo.empty()){
			trafo = cv::Mat();
			return 0;
		}

		if (Execution::verbosity >=2){
			std::cout << "size of from and to vectors: " << num_pairs << std::endl;
			std::cout << "estimated transformation matrix:\n"
					  << printTrafo(final_trafo);
			std::cout << "number of inliers " << max_inliers << std::endl;
		}
		if (Execution::verbosity >= 1)
			std::cout << "<- Compute transformation matrix using RANSAC: "
					  << max_inliers << " inliers\n\n";

		// assign only inliers as new from and to
		std::vector<cv::Vec2f> f(max_inliers);
		std::vector<cv::Vec2f> t(max_inliers);
		int cnt = 0;
		for (int i = 0; i < num_pairs; i++){
			if (inliers[i] > 0){
				f[cnt] = from.at<cv::Vec2f>(i,0);
				t[cnt] = to.at<cv::Vec2f>(i,0);
				cnt++;
				if (cnt == max_inliers){
					break;
				}
			}
		}
		// we don't want a refinement for all inliers if we have too many of them
		int max_inliers_for_gold = 1000;
		if ( max_inliers > max_inliers_for_gold ){
			std::vector<cv::Vec2f> f_rand(max_inliers_for_gold);
			std::vector<cv::Vec2f> t_rand(max_inliers_for_gold);
			for (int i = 0; i < max_inliers_for_gold; i++){
				int index = rand() % max_inliers;
				f_rand[i] = f[index];
				t_rand[i] = t[index];
			}
			f = f_rand;
			t = t_rand;
		}

		// compute a refined transformation matrix between f and t
		// skip that if we have enough inliers
		vout(2) << "\trefine trafo by gold-standard\n";
		Verification::affinity( cv::Mat(f), cv::Mat(t), trafo );
		num_inliers = max_inliers;
	} // end own ransac
	// pure refinement due to goldstandard
	else if ( estimateMethod == "goldstandard" ){
		vout(2) << " estimate trafo with the gold-standard algorithm, having "
				<< from.rows << " points, from which "
				<< std::min(from.rows,num_iterations) << " are taken\n";
		// let's adjust the maximum from and to with num_iterations
		// TODO: maybe this double-interpretation of num_iterations is not so good
		int max_inliers_for_gold = num_iterations;
		if ( from.rows > max_inliers_for_gold ){
			std::vector<cv::Vec2f> f_rand(max_inliers_for_gold);
			std::vector<cv::Vec2f> t_rand(max_inliers_for_gold);
			for (int i = 0; i < max_inliers_for_gold; i++){
				int index = rand() % from.rows;
				f_rand[i] = from.at<cv::Vec2f>(index,0);
				t_rand[i] = to.at<cv::Vec2f>(index,0);
			}
			Verification::affinity( cv::Mat(f_rand), cv::Mat(t_rand), trafo );
		}
		else {
			Verification::affinity( from, to, trafo );
		}
		vout(2) << printTrafo(trafo);
	}

	// de-normalize the homography matrix
	if ( !trafo.empty() && normalize ){
		denormalizeH(t1, t2, trafo);
	}

	return num_inliers;
}

void Verification :: convertCorresToShift(void)
{
	if (!shiftvec.empty()){
		return;
	}	
	// convert correspondence matrix to shift-vector multi-map
	shiftvec.resize(corres.size());
	for (size_t i = 0; i < corres.size(); i++){
		for (int y = 0; y < corres[i].rows; y++){
			for (int x = 0; x < corres[i].cols; x++){
				if (corres[i](y,x).x == -1){
					continue;
				}
				cv::Point pa(x,y);
				cv::Point pb(corres[i](y,x));
				int dx = std::abs(pa.x - pb.x);
				int dy = std::abs(pa.y - pb.y);
				shiftvec[i].insert(std::pair< std::pair<int, int>,
								   std::pair<cv::Point,cv::Point> >
								   (std::make_pair(dx,dy),
									std::make_pair(pa,pb)) );
			}
		}
	}
}

void Verification :: verificateShift(size_t minSameShift, int numMax, int shiftVariance)
{
	convertCorresToShift();

	if (Execution::verbosity >= 2){
		std::cout << "check shift\n";
		std::cout << "chans of shiftvecs: " << shiftvec.size() << std::endl;
	}

	double t = (double)cv::getTickCount();

	std::multimap<std::pair<int,int>, std::pair<cv::Point, cv::Point> >::iterator it;

	// delete in block set the ones with too less shift vectors
	std::pair<std::multimap<std::pair<int,int>, std::pair<cv::Point,cv::Point> >::iterator,
			std::multimap<std::pair<int,int>, std::pair<cv::Point,cv::Point> >::iterator> toDel;

	for (size_t i = 0; i < shiftvec.size(); i++)
	{
		if (shiftvec[i].begin() == shiftvec[i].end()){
            std::cerr << "Warning: Verification::verificate: there are no shiftvectors -> nothing matched\n";
			continue;
		}

		// vector of all maxima
		std::vector<size_t> max;
		max.push_back(0);
		size_t maxTh = 1;

		// check nr of max
		if (numMax > 0 && shiftVariance == 0){
			// get the shift vector with the highest frequency
			for (it = shiftvec[i].begin() ; it != shiftvec[i].end(); ) // ++it)
			{
				std::pair<int, int> tmp = it->first;
				size_t cnt = shiftvec[i].count(tmp);
				toDel = shiftvec[i].equal_range(tmp);
				if (cnt > minSameShift){
					max.push_back(cnt);
				}
				else // erase already these which have too small minSameShift
				{
					shiftvec[i].erase(toDel.first, toDel.second); // remove all pairs<int, int> in the shiftvec mmap
				}
				it = toDel.second;
			}

			// sort the maximum vector
			std::sort(max.begin(), max.end(), cmp);

			// delete duplicates
			std::vector<size_t>::iterator maxItEnd;
			maxItEnd = std::unique(max.begin(), max.end());
			max.erase(maxItEnd, max.end());

			maxTh = std::min(static_cast<size_t>(numMax), max.size());

			if (Execution::verbosity >= 1)
			{
				if (!max.empty()){
					if ( Execution::verbosity >= 1 ){
						std::cout << "Veri: accept the following shift vectors:\n";
					}
					for (size_t k = 0; k < maxTh; k++){
						if (max[k] >= minSameShift){
							std::cout << "\t" << max[k] << std::endl;
						}
					}
				}
				else{
					if ( Execution::verbosity >= 2 )
						std::cout << "Veri: no shift vector found " << std::endl;
					continue;
				}
			}

			for ( it = shiftvec[i].begin(); it != shiftvec[i].end(); )
			{
				std::pair<int, int> tmp = it->first;
				size_t cnt = shiftvec[i].count(tmp);
				toDel = shiftvec[i].equal_range(tmp);
				if (cnt < max[maxTh-1])
				{
					shiftvec[i].erase(toDel.first, toDel.second); // remove all pairs<int, int> in the shiftvec mmap
				}
				it = toDel.second;
			}
		}
		else if ( shiftVariance > 0 )
		{
			for (int var = 0; var < shiftVariance; var++)
			{
				// compute matrix which holds the counts of shiftvectors
				cv::Mat_<double> shiftCntMatrix(dimy, dimx);
				// distance allowed between two shiftvecs
				//int windowSize = shiftVariance; // actually windowsize = half of the window
				int windowSize = var; // actually windowsize = half of the window
				// contains sums of windows with 'windowSize'
				cv::Mat_<double> sum(dimy, dimx, 0.0);
				for ( it = shiftvec[i].begin(); it != shiftvec[i].end(); )
				{
					std::pair<int, int> tmp = it->first;
					size_t cnt = shiftvec[i].count(tmp);
					toDel = shiftvec[i].equal_range(tmp);
					it = toDel.second;
					shiftCntMatrix(tmp.second, tmp.first) = static_cast<double>(cnt);
				}
				// compute integral image out of shiftCntMatrix
				cv::Mat_<double> integral; // with size: (dimy+1) x (dimx+1)
				cv::integral(shiftCntMatrix, integral);
				// fast computation of sum of windows
				for (int y = windowSize; y < sum.rows-windowSize; y++) {
					for (int x = windowSize; x < sum.cols-windowSize; x++) {
						// fast computation of a sum of a window
						sum(y,x) = integral(y+windowSize, x+windowSize)
								- integral(y+windowSize, x-windowSize)
								- integral(y-windowSize, x+windowSize)
								+ integral(y-windowSize, x-windowSize);
					}
				}
				double max = 0.0;
				cv::minMaxLoc(sum, NULL, &max, NULL, NULL);
				if (max < minSameShift){
					continue;
				}
				if (Execution::verbosity >= 1){
					std::cout << "max: " << max << std::endl;
				}
				if (numMax > 0){
					cv::Point maxLoc;
					if (numMax > 1){
						cv::Mat_<double> tmp = sum.clone();
						for (int i = 0; i < numMax; i++){
							cv::minMaxLoc(tmp, NULL, &max, NULL, &maxLoc);
							tmp(maxLoc.y, maxLoc.x) = 0.0;
						}
					}
					else
						cv::minMaxLoc(sum, NULL, &max, NULL, &maxLoc);
				}
				// remove all pairs<int, int> in the shiftvec mmap
				for ( it = shiftvec[i].begin(); it != shiftvec[i].end(); ){
					std::pair<int, int> pos = it->first;
					toDel = shiftvec[i].equal_range(pos);
					if ( sum(pos.second, pos.first) < (numMax > 0 ? max : static_cast<double>(minSameShift)) )
					{
						shiftvec[i].erase(toDel.first, toDel.second);
					}
					it = toDel.second;
				}
				break;
			}
		}
		else
		{ // just delete the ones which have a too small number of same shift vecs
			for ( it = shiftvec[i].begin(); it != shiftvec[i].end(); )
			{
				std::pair<int, int> tmp = it->first;
				size_t cnt = shiftvec[i].count(tmp);
				toDel = shiftvec[i].equal_range(tmp);
				if (cnt < minSameShift)
				{
					// remove all pairs<int, int> in the shiftvec mmap
					shiftvec[i].erase(toDel.first, toDel.second);
				}
				it = toDel.second;
			}
		}
	} // end for chans

	t = ((double)cv::getTickCount() - t)/cv::getTickFrequency();
	if ( Execution::verbosity >= 2 ){
		std::cout << "time neaded for checking shift vecs: " << t << std::endl;
	}

	if (Execution::verbosity >= 1)
		std::cout << "Check shift done" << std::endl;
}

// TODO revise that!
cv::Mat& Verification :: getSimilarityMatrix(void)
{
	// is only empty if verificateShift() is used
	// TODO: really?
    if (matrix.empty())
	{
		if (Execution::verbosity >= 2)
			std::cout << "Get similarity matrix\n";

		//convertCorresToShift();
        if ( ! shiftvec.empty() ){

			std::vector<cv::Mat> planes;
			planes.resize(shiftvec.size());

			for (size_t c = 0; c < shiftvec.size(); c++){
				std::multimap<std::pair<int,int>, std::pair<cv::Point, cv::Point> >::iterator it;

				cv::Mat_<uchar> tmp(dimy, dimx, static_cast<uchar>(0));
				for(it = shiftvec[c].begin(); it != shiftvec[c].end(); ++it){
					int value = 1;

					//				size_t blockSize = 1;
					{ // TODO: make an option for that
						cv::Point p1 = it->second.first;
						cv::Point p2 = it->second.second;
                        // could also mark it that way to see which points got the most matches
//						if (static_cast<int>(p1.x + block_size/2) < dimx
//								&& static_cast<int>(p1.y + block_size/2) < dimy){
//							if ((tmp(p1.y + block_size/2, p1.x + block_size/2) + value) <= 255)
//								tmp(p1.y + block_size/2, p1.x + block_size/2) += value;
//						}
//						if (static_cast<int>(p2.x + block_size/2) < dimx
//								&& static_cast<int>(p2.y + block_size/2) < dimy){
//							if ((tmp(p2.y + block_size/2, p2.x + block_size/2) + value) <= 255)
//								tmp(p2.y + block_size/2, p2.x + block_size/2) += value;
//						}
                        tmp(p1.y + block_size/2, p1.x + block_size/2) = 255;
                        tmp(p2.y + block_size/2, p2.x + block_size/2) = 255;
					}

				}

				planes[c] = tmp;
			}
			if (planes.size() != 1){
				if (Execution::verbosity >=2){
					std::cout << "Verification::getSimilaritymatrix(): merge now planes";
				}
				cv::merge(planes, matrix);
			}
			else{
				matrix = planes[0];
			}

			// write a dumpfile
			if (Execution::verbosity >= 3) {
				std::cout << "imwrite detection matrix real: matrix2.png\n";
				std::stringstream ss;
				ss << Execution::outputdir << "matrix2.png";
				cv::imwrite(ss.str(), matrix);
			}
		}

		else {
			matrix = cv::Mat_<uchar>(dimy, dimx, (uchar) 0);
			for (int i = 0; i < (int)corres.size(); i++){
				for (int y = 0; y < corres[i].rows; y++){
					for (int x = 0; x < corres[i].cols; x++){
						if (corres[i](y,x).x != -1){
							matrix.at<uchar>(y,x) = 255;
						}
					}
				}
			}
		}
	} 
	return matrix;
}

void Verification :: printSet(void) const {
	std::cout << "----PRINT similar block pairs -----\n";
	for (unsigned int i = 0; i < shiftvec.size(); i++)
		std::cout << "num of blockpairs: " << shiftvec[i].size() << std::endl;
	/*	std::set<blockpair>::iterator it;
  for (it = bset.begin(); it != bset.end(); ++it){
  std::cerr << "( " << it->first->getX() << " , " << it->first->getY() << " )";
  std::cerr << "  ~ ";
  std::cerr << " ( " << it->second->getX() << " , " << it->second->getY() << " )";
  std::cerr <<"\t[";
  unsigned int fsize = it->first->getFeatureSize();
  for (unsigned int k = 0; k < fsize; k++){
  std::cerr << it->first->getFeatureAt(k) << " ";
  }
  std::cerr << "] ~ [";
  for (unsigned int k = 0; k < fsize; k++){
  std::cerr << it->first->getFeatureAt(k) << " ";
  }
  std::cerr << "]\n";

  }
  */
	std::cout <<"----PRINT of similar block pairs done ----\n";
}

