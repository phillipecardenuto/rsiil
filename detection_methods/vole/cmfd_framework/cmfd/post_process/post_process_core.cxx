#include "post_process_core.h"
#include <queue>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <ctime>
PostProcess :: PostProcess(const PPConfig & c, const cv::Mat & _cmfd_img)
	: config(c)
{
	cmfd_img = _cmfd_img.clone();
}

PostProcess :: ~PostProcess(void){}

// check in a breadth-depth search manner the size of a
// detected region, if this is too small, then unmark it
void PostProcess :: applyAreaThreshold(const cv::Mat & matrix,
									   std::vector<cv::Mat_<cv::Point> > * _corres)
{
	if ( !matrix.empty() ){
		CV_Assert(matrix.rows == cmfd_img.rows && matrix.cols == cmfd_img.cols);
	}
	int min_area_size;
	// if we have an area threshold smaller than 1 treat it as a percentage value
	if (config.areaThreshold < 1.0)
		min_area_size = config.areaThreshold * cmfd_img.cols * cmfd_img.rows;
	else
		min_area_size = static_cast<int>(config.areaThreshold);
	if (config.verbosity >= 2){
		std::cout << "\tmin_area_size = " << min_area_size
				  << " cmfd_img.channels(): " << cmfd_img.channels() << std::endl;
	}

	std::vector<cv::Mat> planes;
	cv::split(cmfd_img, planes);
//    bool * tmp = new bool[cmfd_img.rows*cmfd_img.cols];
    cv::Mat_<bool> visited(cmfd_img.rows, cmfd_img.cols, false);
	for (int c = 0; c < cmfd_img.channels(); c++) 
	{
		cv::Mat_<uchar> cmfd_img_plane = planes[c];
//        memset(tmp, false, cmfd_img.rows * cmfd_img.cols * sizeof(bool) );
//        cv::Mat_<bool> visited(cmfd_img.rows, cmfd_img.cols, tmp);

        visited.setTo(false);
        for ( int y = 0; y < cmfd_img.rows; ++y ) {
            for ( int x = 0; x < cmfd_img.cols; ++x )
			{
				if ( visited(y,x) || cmfd_img_plane(y,x) == 0 )
					continue;
			
				visited(y,x) = true;
				std::queue<cv::Point> queue;
				std::vector<cv::Point> path;
				queue.push(cv::Point(x,y));
				path.push_back(cv::Point(x,y));
				
				while( !queue.empty() )
				{
					cv::Point node = queue.front();
					queue.pop();
					std::vector<cv::Point> neighbors = getNeighbors(cmfd_img_plane, node);

					for (size_t i = 0; i < neighbors.size(); i++){
						if ( !visited(neighbors[i].y, neighbors[i].x)
							 && cmfd_img_plane(neighbors[i].y, neighbors[i].x) != 0 ){
							queue.push(neighbors[i]);
							path.push_back(neighbors[i]);
							visited(neighbors[i].y, neighbors[i].x) = true;
						}
					}
				}
				// if the area is too small unmark it		
				if ( path.size() < static_cast<size_t>(min_area_size) ){
					unmark(path, c);
					continue;
				}
				if ( !matrix.empty() ){
					// go through the path and check if we have a feature on that path
					bool found_feature = false;
					for (size_t i = 0; i < path.size(); i++){
						if ( (matrix.channels() == 1)
							 && (matrix.at<uchar>(path[i].y,path[i].x) > 0) ){
							found_feature = true;
							break;
						}
						else if ( (matrix.channels() == 3)
								  // maybe just the channel c is too harsh here...
								  && (matrix.at<cv::Vec3b>(path[i].y,path[i].x)[c] > 0) ){
							found_feature = true;
							break;
						}
					}
					if ( !found_feature ){
						unmark(path, c);
						continue;
					}
				}
				if ( _corres != NULL) {
					bool found_marked_corres = false;
					for (size_t i = 0; i < path.size(); i++){
						cv::Point p = path[i];
						// have we a correspondence?
						for (size_t t = 0; t < _corres->size(); t++){
							cv::Point p_c = (*_corres)[t](p.y, p.x);
							if (p_c.x == -1) // no correspondence
								break;
							// ... and is it marked?
							if ( cmfd_img_plane(p_c.y, p_c.x) > 0 ){
								found_marked_corres = true;
								break;
							}
						}
						if ( found_marked_corres )
							break;
					}
					// no correspondence found -> unmark that area
					if ( !found_marked_corres ){
						if (config.verbosity > 2)
							std::cout << "area has no correspondence area -> unmark\n";
						unmark(path, c);
					}
				} // end if
			} // end for img.cols
		} // end for img.rows  
	} // end for img.channels()
//    delete[] tmp;
}

void PostProcess :: unmark(const std::vector<cv::Point> & path, int c)
{
	for (size_t i = 0; i < path.size(); i++){
		if (cmfd_img.channels() == 1){
			cmfd_img.at<uchar>(path[i].y, path[i].x) = 0;
		}
		if (cmfd_img.channels() == 3){
			cmfd_img.at<cv::Vec3b>(path[i].y, path[i].x)[c] = 0;
		}
	}
}

void PostProcess :: applyMorphFilter()
{
	// determine iterations
	int iterations = 0;
	if (config.morphIter.size() == 0)
		iterations = 1;
	else if (config.morphIter.size() == 1)
		iterations = config.morphIter[0];
	else if (config.morphIter.size() != config.morphOp.size()){
		throw std::runtime_error("morphIter.size() != morphOp.size()");
	}

	for (size_t i = 0; i < config.morphOp.size(); i++){
		// # of iterations
		int iter = (iterations == 0 ? config.morphIter[i] : iterations);
		std::string op = config.morphOp[i];
		if (op == "ERODE"){
			cv::Mat tmp;
			cv::erode(cmfd_img,		// src
					tmp,			// dst
					cv::Mat(),	// element (cv::Mat() == 3x3)
					cv::Point(-1,-1), // anchor-point (-1,-1):center
					iter);
			cmfd_img = tmp;
		}
		else if (op == "DILATE"){
			cv::Mat tmp;
			cv::dilate(cmfd_img,	// src
					tmp,			// dst
					cv::Mat(),		// element (cv::Mat() == 3x3)
					cv::Point(-1,-1), // anchor-point (-1,-1):center
					iter);
			cmfd_img = tmp;
		}
		// more advanced operations
		else if (op == "OPEN"		// dilate(erode(src,element))
				|| op == "CLOSE"	// erode(dilate(src,element))
				|| op == "ADIEND" // dilate(src,element) - erode(src,element)
				|| op == "TOPHAT"	// src - open(src,element)
				|| op == "BLACKHAT") // close(src,element) - src
		{
			// determine operation and set the cv-type accordingly (see comments above)
			int op_cv;
			if (op == "OPEN")		op_cv = cv::MORPH_OPEN;
			if (op == "CLOSE")		op_cv = cv::MORPH_CLOSE;
			if (op == "GRADIEND")	op_cv = cv::MORPH_GRADIENT;
			if (op == "TOPHAT")		op_cv = cv::MORPH_TOPHAT;
			if (op == "BLACKHAT")	op_cv = cv::MORPH_BLACKHAT;
			
			cv::Mat tmp;
			cv::morphologyEx(cmfd_img, // src
					tmp,	// dst
					op_cv,		// the operation
					cv::Mat(), // element let's hope that works here too
					cv::Point(-1,-1), // anchor point (see above)
					iter);
			cmfd_img = tmp;
		}
	} // end for
}

void PostProcess :: fillHolesByContour(std::vector<std::vector<cv::Point> > & contours,
									   std::vector<cv::Vec4i> & hierarchy )
{
	CV_Assert(!cmfd_img.empty());
	cv::Mat_<uchar> img;
	/// merge channels if we have a 3-channel cmfd map
	if (cmfd_img.channels() == 3){
		std::vector<cv::Mat> planes;
		cv::split(cmfd_img, planes);
		img = planes[0] + planes[1] + planes[2];
	}
	else{
		// we have to clone here as the image will be modified in
		// findContours!!
		img = cmfd_img.clone();
	}

	if (contours.empty()){
		cv::findContours(img, contours, hierarchy,
						 CV_RETR_EXTERNAL, // no inner contours
						 CV_CHAIN_APPROX_NONE); // should the chain be approximated
	}
	cv::drawContours(cmfd_img, contours, -1,
					 (cmfd_img.channels() == 1 ? cv::Scalar(255) : cv::Scalar(255,255,255)), // color
					 CV_FILLED // filled out
					 );
}

std::vector<cv::Point> PostProcess :: getNeighbors(const cv::Mat &img, const cv::Point &p) const
{
	std::vector<cv::Point> neighbors;
	// add the neighbors
	if ( p.y+1 < img.rows )
		neighbors.push_back(cv::Point(p.x,	 p.y+1));	// mid below
	if ( p.y-1 >= 0)
		neighbors.push_back(cv::Point(p.x,	 p.y-1));	// mid above
	if ( p.x-1 >= 0)
		neighbors.push_back(cv::Point(p.x-1, p.y));		// left mid
	if ( p.x+1 < img.cols)
		neighbors.push_back(cv::Point(p.x+1, p.y));		// right mid

	// for eight-neighborhood:
	if ( p.x-1 >= 0 && p.y+1 < img.rows )
		neighbors.push_back(cv::Point(p.x-1, p.y+1));	// left below
	if ( p.x+1 < img.cols && p.y+1 < img.rows )
		neighbors.push_back(cv::Point(p.x+1, p.y+1));	// right below
	if ( p.x-1 >= 0 && p.y-1 >= 0 )
		neighbors.push_back(cv::Point(p.x-1, p.y-1));	// left above
    if ( p.x+1 < img.cols && p.y-1 >= 0 )
        neighbors.push_back(cv::Point(p.x+1, p.y-1));	// right above

	return neighbors;
}
