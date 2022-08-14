/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef BUILDCLUSTER_H
#define BUILDCLUSTER_H

#include <opencv2/core/core.hpp>
#include <string>
#include <vector>

namespace cluster{

struct Similarity{
	Similarity(){}
	Similarity(float _dist, int _ind_a, int _ind_b){
		dist = _dist;
		ind_a = _ind_a;
		ind_b = _ind_b;
	}
	float dist;
	// cluster indices
	int ind_a;
	int ind_b;
};

struct LessThanIdx {
	bool operator()(const Similarity & a, const Similarity & b) const{
		return a.dist < b.dist;
	}
};

class Cluster{
public:
	explicit Cluster(){
		id = getid();
	}
	explicit Cluster(std::pair<cv::Point,cv::Point> point){
		id = getid();
		add(point);
	}
	void add(std::pair<cv::Point,cv::Point> c){
		objects.push_back(c);
	}

	cv::Point2d getClusterCenter() const{
		cv::Point2d center(0.0, 0.0);
		for (size_t i = 0; i < objects.size(); i++){
			center.x += objects[i].first.x;
			center.y += objects[i].first.y;
		}
		center *= 1.0/objects.size();
		return center;
	}
	void merge(Cluster & b){
		objects.insert(objects.end(),
					   b.objects.begin(),
					   b.objects.end());
	}	
	std::vector< std::pair<cv::Point,cv::Point> > objects;
	int id;
private:
	int getid(){
		static int count = 0;
		return count++;
	}
};

class BuildCluster
{
public:
	BuildCluster(const std::vector<cv::Mat_<cv::Point> > & corres,
				 float _cut_threshold,
				 const std::string _linkage_method = "single");
	inline void setCutThreshold(float _cut_threshold) {
		cut_threshold = _cut_threshold;
	}
	/// executes clustering
	void compute( const std::vector<std::vector<std::pair<cv::Point,
				  cv::Point> > > & point_groups );
	void compute();
	/// return the matrix with the clusters marked in different
	/// colors
	cv::Mat getMarkedMatrix(int shift){
		markClustersInMatrix(shift);
		return cluster_matrix;
	}
	/// matrix where each cluster has its own id
	cv::Mat_<uchar> getMatrix();
	/// returns clusters of points
	std::vector<std::vector<std::pair<cv::Point,cv::Point> > > getPointGroups();
	/*! from the computed clusters it returns the clusters of points
	  * grouped by correspondences (like sats also returns the clusters)
	  */
	std::vector<std::vector<std::pair<cv::Point,cv::Point> > > getSplitPointGroups();

private:
	/// the actual computation method
	void _compute();
	/// chose the correct distance
	float getDistance(const Cluster & a, const Cluster & b, const std::string & linkage_method);
	/// remove too small clusters
	void removeSmallClusters();
	/// returns the links of a specific cluster to all other clusters
	/// i.e. vector[x] contains all pairs (point & correspondence) of cluster c to cluster x
	std::vector<std::vector<std::pair<cv::Point,cv::Point> > >
		getLinks(const cv::Mat_<uchar> & cluster_matrix, const Cluster & c);

	// linkage methods:
	/*! Single Linkage
	  \return dist(a,b) = min(||x_a_i,x_a_j)||)
	*/
	float singleLink(const Cluster & a, const Cluster & b) const;
	/*! Ward's Linkage
	 * ESS: Error sum of squares
	 *\return dist(a,b) = ESS(a+b)-ESS(a)-ESS(b)
	 */
	float wardsLink(const Cluster & a, const Cluster &) const;
	/// helper method for wardsLink()
	float getSumOfSquares(const Cluster &c) const;
	/// marks the cluster
	void markClustersInMatrix(int shift);

	/// correspondence map
	const std::vector<cv::Mat_<cv::Point> > & corres;
	/// output matrix
	cv::Mat_<cv::Vec3b> cluster_matrix;
	/// the clusters
	std::vector<cv::Ptr<Cluster> > clusters;
	/// which cluster is active
	//std::vector<bool> active;
	/// where to cut the tree (although we don't
	/// set up a real tree here)
	float cut_threshold;
	/// linkage method decides how to merge clusters
	const std::string linkage_method;
	/// maximum number of clusters
	int max_clusters;
	/// minimum links to another cluster
	int min_links;
};

} // end namespace cluster
#endif // BUILDCLUSTER_H
