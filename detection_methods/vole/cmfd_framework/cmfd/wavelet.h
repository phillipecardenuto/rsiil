/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef WAVELET_H
#define WAVELET_H

#include "cv.h"
#include <vector>

class Wavelet
{
	public:
		Wavelet(void){}

		/*! Computes a wavelet-decomposition on an image.
		 *
		 * \param img  the image to decompose
		 * \param dest return-location for the subbands of img
		 * \param n    number of decomposition levels
		 *
		 * The subbands are stored in the dest parameter per channel.
		 * per channel order: HH, LH, HL, LL with ascending order,
		 * e.g. n = 2: HH_1,LH_1,HL_1,HH_2,LH_2,HL_2,LL2
		 */
		static void decompose(const cv::Mat &img, std::vector<cv::Mat_<double> > &dest,
				const cv::Mat_<double> &maska, const cv::Mat_<double> &maskb, int n, bool normalize = false);
		
		/*! Computes a wavelet-decomposition on an image.
		 *
		 * \param img  the image to decompose
		 * \param dest subbands build up in an image
		 * \param n    number of decomposition levels
		 *
		 * The subbands are stored in the dest parameter.
		 */
		static void decompose(const cv::Mat &img, cv::Mat &destimg, 
				const cv::Mat_<double> &maska, const cv::Mat_<double> &maskb, int n, bool normalize = false);

		/// shortcut for haar (= daubechies D2) wavelet decomposition
		static void decomposeHaar(const cv::Mat &img, cv::Mat &destimg, int n);
	
		/// shortcut for haar (= daubechies D2) wavelet decomposition
		static void decomposeHaar(const cv::Mat &img, std::vector<cv::Mat_<double> > &dest, int n, bool normalize = false);

		/*! Computes a subband of an image using the given filters.
		 *
		 * \param img         image to compute subband from
		 * \param dest        return-location (list) for the subband
		 * \param maska       horizontal filter
		 * \param maskb       vertical filter
		 *
		 * The computed subband is appended to the list of subbands dest.
		 */
		static void subband(const cv::Mat_<uchar> &img, std::vector<cv::Mat_<double> > &dest, 
				const cv::Mat_<double> &maska, const cv::Mat_<double> &maskb, bool normalize = false); 
	private:
		/// helper to recursivly build up the image from the vector of subbands
		static void rek(std::vector<cv::Mat_<double> > &vec, cv::Mat_<double> &m);
};
#endif
