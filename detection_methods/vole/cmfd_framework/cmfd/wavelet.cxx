/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "wavelet.h"

#include "highgui.h"
#include <iostream>

void Wavelet :: subband(const cv::Mat_<uchar> &img, std::vector<cv::Mat_<double> > &dest, 
						const cv::Mat_<double>  &maska, const cv::Mat_<double>  &maskb, bool normalize)
{
	// convolve image
	cv::Mat_<double> out;
	cv::Mat_<double> imgdouble = img;
	cv::sepFilter2D( imgdouble, out, -1, maska, maskb);

	// TODO do I need this? /= (2*sqrt(0.5);
	double sum = (cv::sum(maska)[0] + cv::sum(maskb)[0]);
	if (sum > 0)
		out /= sum;

	// scale image
	if ( out.rows > 1 && out.cols > 1 ) { // only resize if posible
		cv::Mat_<double> scaled;
		cv::resize(out, scaled, cv::Size(), 0.5, 0.5, cv::INTER_NEAREST);

		if (normalize){
			double max;
			double min;
			cv::minMaxLoc(scaled, &min, &max);			
			scaled -= min; // move min to 0
			scaled *= (255 / (max - min)); // scale  to 0-255
		}

		dest.push_back(scaled);
	}
}

void Wavelet :: decompose(const cv::Mat &img, std::vector<cv::Mat_<double> > &dest,
						  const cv::Mat_<double>  &high, const cv::Mat_<double>  &low, int n, bool normalize)
{
	if (img.empty()) 
		return;

	cv::Mat_<uchar>  band = img;
	for (int i = 0; i < n; i++) // for every depth compute the four subbands
	{
		subband(band, dest, high, high, normalize);
		subband(band, dest, low,  high, normalize);
		subband(band, dest, high, low, normalize);
		subband(band, dest, low,  low, normalize );

		band = dest.back();

		/* forget the last subband of every level (the low-low filtered one),
		 * only keep the the one of the last level (rest)
		 */
		if (i < n - 1) {
			dest.pop_back();
		}
	}

}

void Wavelet :: decompose(const cv::Mat &img, cv::Mat &destimg,
						  const cv::Mat_<double>  &high, const cv::Mat_<double>  &low, int n, bool normalize)
{
	std::vector<cv::Mat_<double> > vec;
	decompose(img, vec, high, low, n, normalize);

	cv::Mat_<double> tmp = vec.back();
	vec.pop_back();
	// set the levels together again
	rek(vec, tmp);

	cv::Mat_<uchar> uu = tmp; // cast to uchar
	destimg = uu;
}

void Wavelet :: rek(std::vector<cv::Mat_<double> > &vec, cv::Mat_<double> &m)
{
	if (vec.empty()) // recursion end
		return;

	cv::Mat_<double> a(m.rows * 2, m.cols*2);

	for (int j = 0; j < m.rows; j++){
		for (int i = 0; i < m.cols; i++){
			a(j, i) = m(j,i);
			a(j, m.cols + i ) = vec[vec.size()-1](j,i);
			a(m.rows + j, i) = vec[vec.size()-2](j,i);
			a(m.rows + j, m.cols + i) =vec[vec.size()-3](j,i);
		}
	}

	vec.erase(vec.end() - 3, vec.end());
	m = a;
	rek(vec, m);
}

void Wavelet :: decomposeHaar(const cv::Mat &img, cv::Mat &destimg, int n)
{
	/* setup high- and low-pass filters */
	cv::Mat_<double> high(1,2, sqrt(0.5));
	high(0,1) *= -1;
	cv::Mat_<double> low(1,2, sqrt(0.5));

	decompose(img, destimg, high, low, n, true);
}

void Wavelet :: decomposeHaar(const cv::Mat &img, std::vector<cv::Mat_<double> >  &dest, int n, bool normalize)
{
	/* setup high- and low-pass filters */
	cv::Mat_<double> high(1,2, sqrt(0.5));
	high(0,1) *= -1;
	cv::Mat_<double> low(1,2, sqrt(0.5));

	decompose(img, dest, high, low, n, normalize);
}


