/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "buildcluster.h"
#include "execution.h"
#include "matching.h"
#include "mark.h"

#include <iostream>
#include <queue>
#include <opencv2/highgui/highgui.hpp>

namespace cluster{

BuildCluster :: BuildCluster( const std::vector<cv::Mat_<cv::Point> > & _corres,
							  float _cut_threshold,
							  const std::string _linkage_method)
	: corres(_corres),
	  linkage_method(_linkage_method)
{
	// TODO make perhaps an option out of these
	max_clusters = 10;
	min_links = 4;
	if (_cut_threshold < 1.0){
		cut_threshold = _cut_threshold * cluster_matrix.rows * cluster_matrix.cols;
	}
	else {
		cut_threshold = _cut_threshold;
	}
	cluster_matrix = cv::Mat_<cv::Vec3b>::zeros(corres[0].rows, corres[0].cols);
}

// Note we only add one direction here (in opposite to the other compute) as this is normally
// used with keypoint-based methods which should have been run with bruteforce-matcher, thus
// is already very accurate (i.e. both pairs( (x,y) and its correspondence) are with a very high
// probablity matched)
void BuildCluster :: compute()
{
	if (corres.empty())
		return;

	//cv::Mat_<uchar> merged_corres(cluster_matrix.rows, cluster_matrix.cols, (uchar) 0);
	// assign each correspondence point to a cluster
	for (int y = 0; y < corres[0].rows; y++){
		for (int x = 0; x < corres[0].cols; x++){
			for (size_t c = 0; c < corres.size(); c++){
				if (corres[c](y,x).x > 0){
		//			merged_corres(y,x) = 255;
					cv::Point corr = corres[c](y,x);
					Cluster * c = new Cluster( std::make_pair( cv::Point(x,y), corr ) );
					clusters.push_back( cv::Ptr<Cluster>(c) );
				}
				else {
					break;
				}
			}
		}
	}

//	cv::imwrite("merged_corres.png", merged_corres);

	_compute();
}

void BuildCluster :: compute( const std::vector<std::vector<std::pair<cv::Point,
								cv::Point> > > & point_groups )
{
	if (point_groups.empty())
		return;

	/*
	vout(2) << "num point_groups: " << point_groups.size() << std::endl;
	for (size_t i = 0; i < point_groups.size(); i++){
		cv::Mat_<uchar> tmp = cv::Mat_<uchar>::zeros(corres[0].rows, corres[0].cols);
		cv::Mat_<uchar> tmp_c = cv::Mat_<uchar>::zeros(corres[0].rows, corres[0].cols);
		for (size_t k = 0; k < point_groups[i].size(); k++){
			tmp( point_groups[i][k].first.y, point_groups[i][k].first.x ) = 255;
			tmp_c( point_groups[i][k].second.y, point_groups[i][k].second.x ) = 255;
		}
		std::stringstream ss;
		ss << "tmp_point_group_" << i;
		cv::imwrite(std::string(ss.str() + "_0.png"), tmp);
		cv::imwrite(std::string(ss.str() + "_1.png"), tmp_c);
	}
	*/
	clusters.clear();

	// assigns point-groups to clusters, also assign its correspondences to clusters
	for (size_t i = 0; i < point_groups.size(); i++){
		Cluster *c = new Cluster();
		c->objects = point_groups[i];
		clusters.push_back( cv::Ptr<Cluster>(c) );

		// also add the other direction as an independent cluster
		Cluster *c_corres = new Cluster();
		std::vector< std::pair<cv::Point,cv::Point> > other_direction;
		other_direction.resize(point_groups[i].size());
		for (size_t k = 0; k < point_groups[i].size(); k++){
			other_direction[k] = std::make_pair( point_groups[i][k].second,
												 point_groups[i][k].first );
		}
		c_corres->objects = other_direction;
		clusters.push_back( cv::Ptr<Cluster>(c_corres) );
	}

	_compute();
}

void BuildCluster :: _compute()
{
	vout(2) << "num clusters: " << clusters.size() << std::endl;
/*
	for (size_t i = 0; i < clusters.size(); i++){
		cv::Mat_<uchar> tmp = cv::Mat_<uchar>::zeros(corres[0].rows, corres[0].cols);
		for (size_t k = 0; k < clusters[i]->objects.size(); k++){
			tmp( clusters[i]->objects[k].first.y, clusters[i]->objects[k].first.x ) = 255;
		}
		std::stringstream ss;
		ss << "tmp_cluster_" << clusters[i]->id << ".png";
		cv::imwrite(ss.str(), tmp);
	}
*/
	double t = (double)cv::getTickCount();
// TODO: unfortunately this takes too much memory if we have many clusters
/*
	// stores the distances
	std::vector<std::vector<Similarity> > dists(clusters.size());

	// keeps track of active clusters
	std::vector<bool>active(clusters.size(), true);
	// priority queue of distances between clusters
	std::vector< std::multiset<Similarity,LessThanIdx> > priority;
	priority.resize(clusters.size());
	for (size_t n = 0; n < clusters.size(); n++){
		std::vector<Similarity> dists_tmp;
		dists_tmp.reserve(clusters.size()-1);
		for (size_t i = 0; i < clusters.size(); i++){
			if (n == i) continue;
			float distance = getDistance(*(clusters[n]), *(clusters[i]), linkage_method);
			Similarity s(distance,n,i);
			priority[n].insert(s);
			dists_tmp.push_back(s);
		}
		dists[n] = dists_tmp;
//		if (n % 100 == 0){
//			std::cout << "n = " << n << std::endl;
//		}
	}

	if (Execution::verbosity >= 4){
		std::cout << "start now loop\n";
	}
//	while(!priority.empty()){
	while( true ) {
		// get smallest distance from all priority-queues of clusters
		float dist_min = FLT_MAX;
		Similarity min_sim;
		for (size_t i = 0; i < priority.size(); i++){
			if (active[i]){
//				Similarity & s = *priority[i].begin();
//				// maybe we already have de-activated it -> remove
//				if ( !active[s.ind_b] ) {
//					priority[i].erase( priority[i].begin() );
//					i--;
//					continue;
//				}
//				std::cerr << "dist " << s.dist << " <? " << dist_min << std::endl;
//				if (s.dist < dist_min){
				if (priority[i].begin()->dist < dist_min){
					dist_min = priority[i].begin()->dist;
					min_sim = *(priority[i].begin());
				}
			}
		}
		if (Execution::verbosity >= 4){
			std::cout << "min_dist = " << dist_min << std::endl;
//			std::cout << " a: " << min_sim.ind_a << " b: " << min_sim.ind_b << std::endl;
		}

		// can we stop merging?
		if ( dist_min == FLT_MAX ){
			break;
		}
		if ( dist_min > cut_threshold ){
			break;
		}
		// merge the clusters with lowest distance
		clusters[min_sim.ind_a]->merge(*(clusters[min_sim.ind_b]));

		// deactivate the cluster which has been merged
		active[min_sim.ind_b] = false;
		priority[min_sim.ind_a].clear();

		// remove elements from that cluster
		//clusters.erase(clusters.begin() + min_sim.ind_b);
		for (size_t i = 0; i < clusters.size(); i++){
			if ( !active[i]  || i == min_sim.ind_a )
				continue;
			priority[i].erase(dists[i][min_sim.ind_a]);
			priority[i].erase(dists[i][min_sim.ind_b]);

			float distance = getDistance(*(clusters[i]), *(clusters[min_sim.ind_a]), linkage_method);
			Similarity s1(distance,i,min_sim.ind_a);
			Similarity s2(distance,min_sim.ind_a,i);
			// update
			dists[i][min_sim.ind_a] = s1;
			dists[min_sim.ind_a][i] = s2;

			priority[i].insert(s1);
			priority[min_sim.ind_a].insert(s2);
		}
	}

	std::vector<cv::Ptr<Cluster> > clusters_tmp;
	for (size_t i = 0; i < clusters.size(); i++){
		if (active[i])
			clusters_tmp.push_back(clusters[i]);
	}
	clusters = clusters_tmp;
*/
	// the naive way, slow but not memory-intensive
	// merge clusters while the distance is smaller than a treshold
	// TODO: could be speeded-up by a more proper data structure, e.g. kd-tree/p-tree, etc.
	// TODO: probably the brute-force nn-search of opencv is faster -> try out
	while(true){
		float min_distance = FLT_MAX;
		int ind_a;
		int ind_b;
		// search for the two clusters with the minimum distance
		// TODO: this could probably solved better with an appropriate
		// data-structure (heap, priority-queue) so that one doesn't
		// need to search all the time the whole clusters
		for (size_t i = 0; i < clusters.size(); i++){
			for (size_t k = i+1; k < clusters.size(); k++){
				float distance = FLT_MAX;
				distance = getDistance(*(clusters[i]), *(clusters[k]), linkage_method);

				if (distance < min_distance){
					bool allow_merge = true;

					// check if a correspondences of cluster i is in cluster k
					// if so, we do not allow that these clusters can be merged
					// Note: we have a problem here when we have more correspondence maps
//					for (size_t s = 0; s < clusters[i]->objects.size(); s++){
//						for (size_t c = 0; c < corres.size(); c++){
//							cv::Point corres_p = Matching::getCorres(corres,
//																	 clusters[i]->objects[s],
//																	 c);
//							if (corres_p.x == -1)
//								break;
//							for (size_t t = 0; t < clusters[k]->objects.size(); t++){
//								if (corres_p == clusters[k]->objects[t]){
//									allow_merge = false;
//								}
//							}
//						}
//					}

//					for (size_t s = 0; s < clusters[i]->objects.size(); s++){
//						cv::Point corr = clusters[i]->objects[s].second;
//						for (size_t t = 0; t < clusters[k]->objects.size(); t++){
//							if (corr == clusters[k]->objects[t].first){
//								allow_merge = false;
//							}
//						}
//					}

					if (allow_merge){
						min_distance = distance;
						ind_a = i;
						ind_b = k;
					}
				}
			}
		}		

		// no merge possible
		if (min_distance == FLT_MAX){
			break;
		}
		if (min_distance <= cut_threshold){
			vout(3) << "min_distance " << min_distance << " < " << cut_threshold
					<< " --> merge: " << clusters[ind_a]->id << " and "
					<< clusters[ind_b]->id << std::endl;
			// merge clusters a and b
			clusters[ind_a]->merge(*(clusters[ind_b]));
			// erase cluster b
			clusters.erase(clusters.begin() + ind_b);
		}
//		// adapt cut-threshold
//		else if (clusters.size() > max_clusters){
//			cut_threshold += 5;
//		}
		else {
			break;
		}
	} // end while(true)

	//measure time
	t = ((double)cv::getTickCount() - t)/cv::getTickFrequency();
	if (Execution::verbosity >=2){
		std::cout << "clustering finsished. time neaded: " << t << std::endl;
		std::cout << "whole number of clusters: " << clusters.size() << std::endl;
		std::cout << "cut_threshold: " << cut_threshold << std::endl;
	}
}

// assign every cluster to a unique number
cv::Mat_<uchar> BuildCluster :: getMatrix()
{
	cv::Mat_<uchar> cluster_matrix(corres[0].rows, corres[0].cols, (uchar)0);
	for (size_t i = 0; i < clusters.size(); i++){
		for (size_t k = 0; k < clusters[i]->objects.size(); k++){
			cv::Point p = clusters[i]->objects[k].first;
			cluster_matrix(p.y, p.x) = i + 1;
		}
	}
	/*
	// visualize our clusters and correspondences
	cv::Mat_<cv::Vec3b> visual(cluster_matrix.rows, cluster_matrix.cols, cv::Vec3b::all(0));
	for (int y = 0; y < visual.rows; y++){
		for (int x = 0; x < visual.cols; x++){
			if ( cluster_matrix(y,x) != 0 ){
				visual(y,x) = cv::Vec3b(0,0,255);
				for (size_t i = 0; i < corres.size(); i++){
					cv::Point2f corres_point = corres[i](y,x);
					if (corres_point.x == -1){
						break;
					}
					if (cluster_matrix(corres_point.y, corres_point.x) != 0){
						visual(y,x) = cv::Vec3b(0,255,0);
					}
				}
			}
		}
	}
	cv::imwrite("visual_cluster.png", visual);
	*/
	return cluster_matrix;
}

// this will return individual clusters
std::vector<std::vector<std::pair<cv::Point,cv::Point> > > BuildCluster :: getSplitPointGroups()
{
	vout(1) << "--> BuildCluster::getSplitPointGroups()\n";

	cv::Mat_<uchar> cluster_matrix = getMatrix();

	// remove clusters which have less then min_links connections
	// to another cluster
	std::vector<std::vector<std::pair<cv::Point,cv::Point> > > split_clusters;
	for (size_t i = 0; i < clusters.size(); i++)
	{
		std::vector<std::vector<std::pair<cv::Point,cv::Point> > > links = getLinks( cluster_matrix,
																					 *(clusters[i]) );
		for (size_t k = 1; k < links.size(); k++)
		{
			if (links[k].size() < static_cast<size_t>(min_links)){
				continue;
			}
			split_clusters.push_back(links[k]);
		}
	}

	vout(1) << "<-- BuildCluster::getSplitPointGroups() finished and returned "
			<< split_clusters.size() << " split-clusters" << std::endl;

	if (Execution::verbosity >= 3){
		// visualize our clusters and correspondences
		cv::Mat_<cv::Vec3b> visual(cluster_matrix.rows, cluster_matrix.cols, cv::Vec3b::all(0));
		for (size_t y = 0; y < split_clusters.size(); y++){
			// chose a color for the cluster
			cv::Vec3b color(rand()*200, rand()*200, rand()*200);
			// guarantee that the color is not too dark
			color += cv::Vec3b(55,55,55);
			for (size_t x = 0; x < split_clusters[y].size(); x++){
				visual( split_clusters[y][x].first.y, split_clusters[y][x].first.x ) = color;
			}
		}
		Mark::dumpMatrix("_split_point_groups", visual);
	}
	return split_clusters;
}

std::vector<std::vector<std::pair<cv::Point, cv::Point> > > BuildCluster :: getPointGroups()
{
	removeSmallClusters();

	std::vector<std::vector<std::pair<cv::Point, cv::Point> > > ret;
	for (size_t i = 0; i < clusters.size(); i++){
		/*
		std::vector<cv::Point> tmp;
		for (size_t k = 0; k < clusters[i]->objects.size(); k++){
			tmp.push_back(clusters[i]->objects[k].first);
		}
		ret.push_back(tmp);
		*/
		ret.push_back(clusters[i]->objects);
	}
	return ret;
}

// TODO as every cluster-point now also have his correspondence saved,
// maybe this can be changed?
std::vector<std::vector<std::pair<cv::Point,cv::Point> > > BuildCluster ::
	getLinks(const cv::Mat_<uchar> & cluster_matrix, const Cluster & c)
{
	std::vector<std::vector<std::pair<cv::Point,cv::Point> > > links;
	if ( c.objects.empty() ){
		std::cerr << "WARNING:BuildCluster::getLinks(): cluster c has no points -> no links";
		return links;
	}

	links = std::vector<std::vector<std::pair<cv::Point,cv::Point> > >(clusters.size()+1);
	// compute links to other clusters
	// note link_count(0) = no correspondence for a point found
	// (should actually never happen)
	// but this may happen as in Matching we only add one direction (not both)
	// of a match
	for (size_t k = 0; k < c.objects.size(); k++)
	{

		/*
		cv::Point p = c.objects[k].first;
		// add all possible correspondences
		for (size_t i = 0; i < corres.size(); i++)
		{
			cv::Point corres_point = corres[i](p.y, p.x);
			if (corres_point.x == -1){
				break;
			}

			uchar cluster_nr = cluster_matrix( corres_point.y, corres_point.x );
			// add point to the links of specific cluster_nr;
			links[cluster_nr].push_back( std::make_pair<cv::Point,cv::Point>( p, corres_point ) );
		}
		*/
		uchar cluster_nr = cluster_matrix( c.objects[k].second );
		links[cluster_nr].push_back( c.objects[k] );

	}
	if ( Execution::verbosity >= 3 ){
		//std::cout << "BuildCluster::getLinks():\n";
		for (size_t k = 0; k < links.size(); k++){
			if (links[k].size() > 0)
			//if (links[k].size() > min_links)
				std::cout << "\tlinks[" << k << "]: " << links[k].size() << std::endl;
		}
		//std::cout << std::endl;
	}
	return links;
}
void BuildCluster :: removeSmallClusters()
{
	// assign every cluster to a unique number
	cv::Mat_<uchar> cluster_matrix = getMatrix();

	// remove clusters which have less then min_links connections
	// to another cluster
	std::vector<cv::Ptr<Cluster> > new_clusters;
	// go through all clusters
	for (size_t i = 0; i < clusters.size(); i++)
	{
		// compute links to other clusters
		// note link_count(0) = no correspondence for a point found
		// (should actually never happen)
		std::vector<std::vector<std::pair<cv::Point, cv::Point> > > links = getLinks(cluster_matrix,
																					 *clusters[i]);
		if (links.empty()){
			continue;
		}
		// now remove the links which are neglible, i.e. assign the ones with more than
		// min_links to new ones
		std::vector< std::pair<cv::Point,cv::Point> > objects;
		for (size_t k = 1; k < links.size(); k++){
			if ( links[k].size() < static_cast<size_t>(min_links) ){
				continue;
			}
			for (size_t l = 0; l < links[k].size(); l++){
				objects.push_back(links[k][l]);
			}
		}
		// assign points to new clusters
		if ( ! objects.empty() )
		{
			cv::Ptr<Cluster> c(new Cluster());
			c->objects = objects;
			vout(2) << "assign new cluster with " << objects.size() << " points" << std::endl;
			new_clusters.push_back(c);
		}
	}
	// replace clusters with the reduced clusters
	clusters.clear();
	clusters.insert( clusters.end(), new_clusters.begin(), new_clusters.end() );
	vout(1) << "<-- BuildCluster::removeSmallClusters() returned "
			<< clusters.size() << " clusters" << std::endl;
}

float BuildCluster :: getDistance(const Cluster & a, const Cluster & b, const std::string & linkage_method)
{
	float distance;
	if (linkage_method == "single"){
		distance = singleLink(a, b);
	}
	else if (linkage_method == "ward"){
		distance = wardsLink(a, b);
	}
	else{
		throw std::runtime_error("BuildCluster::getDistance: wrong linkage method given");
	}
	return distance;
}

void BuildCluster :: markClustersInMatrix(int shift)
{
	removeSmallClusters();

	// mark the clusters in the matrix
	srand(time(NULL));
	int count = 0;
	for (size_t i = 0; i < clusters.size(); i++){
		count++;
		// chose a color for the cluster
		cv::Vec3b color(rand()*200, rand()*200, rand()*200);
		// guarantee that the color is not too dark
		color += cv::Vec3b(55,55,55);
		for (size_t k = 0; k < clusters[i]->objects.size(); k++){
			cv::Point & p = clusters[i]->objects[k].first;
			//std::cout << "mark: " << p.x << " " << p.y << std::endl;
			cluster_matrix(p.y + shift, p.x + shift) = /*cv::Vec3b(255,255,255); */color;
		}
	}
	if (Execution::verbosity >=1){
		std::cout << "total number of marked clusters: " << count << std::endl;
	}
}

// Single Linkage Method = smallest distance between objects in the two clusters
float BuildCluster :: singleLink(const Cluster & a, const Cluster & b) const
{
	float min_dist = DBL_MAX;
	for (size_t i = 0; i < a.objects.size(); i++){
		for (size_t k = 0; k < b.objects.size(); k++){
			cv::Point c = a.objects[i].first - b.objects[k].first;
			float dist = cv::norm(c);
			if (dist < min_dist){
				min_dist = dist;
			}
		}
	}
	return min_dist;
}

// Ward's Linkage method
float BuildCluster :: wardsLink(const Cluster & a, const Cluster & b) const
{
	// cluster of a and b
	Cluster ab;
	ab.objects = a.objects;
	ab.objects.insert(ab.objects.end(), b.objects.begin(), b.objects.end());
	// compute sum of sqaures:

	float ward_dist = getSumOfSquares(ab) - (getSumOfSquares(a) + getSumOfSquares(b));
//	std::cerr << ward_dist << " = "
//			  << getSumOfSquares(ab) << " - " << getSumOfSquares(a) << " + "
//			  << getSumOfSquares(b) << std::endl;
	return ward_dist;
}

float BuildCluster :: getSumOfSquares(const Cluster &c) const
{
	cv::Point2d center = c.getClusterCenter();
	float whole_dist = 0.0;
	for (size_t i = 0; i < c.objects.size(); i++){
		cv::Point2d tmp(c.objects[i].first.x-center.x, c.objects[i].first.y-center.y);
		float dist = cv::norm(tmp);
		whole_dist += dist;
	}
	return whole_dist;
}

} // namespace cluster
