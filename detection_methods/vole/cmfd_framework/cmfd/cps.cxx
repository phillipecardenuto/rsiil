/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "cps.h"

#include "highgui.h"
#include "execution.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>

Cps :: Cps(cv::Mat &_img, cv::Mat &_origimg, struct ft_cfg cfg, int _waveletLvl)
: img(_img),
	origimg(_origimg)
{
	divx = cfg.partsInX;
	divy = cfg.partsInY;
	waveletMult = _waveletLvl;
	correlTh = cfg.correlTh;
	pixelDiff = cfg.pixelDiff;
	waveletMult = _waveletLvl * 2;
	if (waveletMult < 1)
		waveletMult = 1;
	matrix = cv::Mat_<cv::Vec3b>::zeros(img.rows * waveletMult, img.cols * waveletMult);
}

// divide in non overlipping sub images
void Cps :: computeSubImg(void)
{

	int width = img.cols / divx;
	//	int restwidth = img.cols % divx;
	int height = img.rows / divy;
	//	int restheight = img.rows % divy;

	// cut subimages out of image
	for (int y = 0; y < divy; y++){
		for (int x = 0; x < divx; x++){
			int newheight = (y+1) * height;
			int newwidth = (x+1) * width;
			// THIS doesnt work, cause mulsprectrums have to be equal dimensional
			//	if (y == divy - 1)
			//		newheight += restheight;
			//	if (x == divx - 1)
			//		newwidth += restwidth;

			subimages.push_back(cv::Mat(img, cv::Range(y * height, newheight),
						cv::Range(x * width, newwidth) ) );
		}
	}
}


void Cps :: crossPowerSpectrum(cv::Mat& src1, cv::Mat& src2, cv::Mat& dest, bool normalize)
{
	// convert to double
	cv::Mat_<double> f1 = src1;
	cv::Mat_<double> f2 = src2;

	cv::Mat F1;
	cv::Mat F2;
	// compute fourier trafo
	cv::dft(f1, F1, cv::DFT_COMPLEX_OUTPUT);
	cv::dft(f2, F2, cv::DFT_COMPLEX_OUTPUT); 

	// compute cross power spectrum of f1 & f2
	cv::Mat P;
	cv::mulSpectrums(F1, F2, P, cv::DFT_COMPLEX_OUTPUT, true); // true: conjugate 2.param

	if (normalize){
		// compute the complex magnitude
		std::vector<cv::Mat> realcomp;
		cv::split(P, realcomp); // split P in real and complex part
		//	cv::Mat real = realcomp[0];
		//	cv::Mat comp = realcomp[1];
		// use only left part cause of DFT_COMPLEX_OUTPUT gains 1 array with 2 chans, 
		// every only half computed
		cv::Mat real(realcomp[0], cv::Range::all(), cv::Range(0, P.cols/2 + 1));
		cv::Mat comp(realcomp[1], cv::Range::all(), cv::Range(0, P.cols/2 + 1));

		// compute complex magnitude
		cv::Mat mag;
		cv::magnitude(real, comp, mag);  

		// normalize cross power spectrum
		real /= mag;
		comp /= mag;
		cv::merge(realcomp, P);
	}

	// inverse fourier transform of P -> write in dest
	cv::idft(P, dest, cv::DFT_REAL_OUTPUT | cv::DFT_SCALE);
}

// FIXME: does not work atm
void Cps :: findRegion(void)
{
	// compute sub images
	computeSubImg();

	// split img into channel planes
	std::vector<cv::Mat> imgplanes;
	cv::split(origimg, imgplanes);

	if (Execution::verbosity >=3)
		std::cerr << "waveletMult: " << waveletMult << std::endl; 
	for (size_t i = 0; i < subimages.size(); i++)
	{
//		for (size_t k = 0; k < subimages.size(); k++)
		for (size_t k = i; k < subimages.size(); k++)
		{
			if (i == k ) continue; // TODO remove this if you know how autocorrel works...
			std::vector<cv::Mat> planes1;
			cv::split(subimages[i], planes1);
			std::vector<cv::Mat> planes2;
			cv::split(subimages[k], planes2);

			std::vector<cv::Mat_<cv::Vec3b> > matrices;
			matrices.resize(4);
			for (int l = 0; l < 4; l++)
			{
				matrices[l] = cv::Mat_<cv::Vec3b>(img.rows*waveletMult, img.cols*waveletMult, cv::Vec3b((uchar)0,(uchar)0,(uchar)0));
			}

			for (size_t c = 0; c < planes1.size(); c++) 
			{
				cv::Mat_<double> theImg = imgplanes[c];

				cv::Mat p;
				//				if (i != k)
				crossPowerSpectrum(planes1[c], planes2[c], p, true); 
				//				else
				//					crossPowerSpectrum(planes1[c], planes2[c], p, false); 
				
				if (Execution::verbosity >=3)
					std::cerr << "p.cols p.rows " << p.cols << " " << p.rows << " ,p.channels() " << p.channels() << " ,p.type() " << p.type() << std::endl;

				// write a dump
				if (Execution::verbosity >= 3){
					std::ostringstream ss;
					ss << "P" << i << k << c;
					std::cerr << "write dump " << ss.str() << std::endl;
					FILE * outfile = fopen( (ss.str()).c_str(), "w" ); // fopen + frwite is faster than ofstream
					//				std::ofstream of((ss.str()).c_str());
					for (int y = 0; y < p.rows; y++){
						for (int x = 0; x < p.cols; x++){
							cv::Mat_<double> tmp = p;
							fprintf(outfile, "%i %i %f.10\n",x,y,tmp(y,x));
							//						of << x << " " << y << " " << std::setprecision(10) << tmp(y,x) << "\n";
						}
						//					of << "\n";
						fprintf(outfile, "\n");
					}
				}

				//--- get maximum peak and spatial location (= position in cross power spectrum)
				cv::Mat_<uchar> mask(p.rows, p.cols, 1);
				// the autocorrelation case
				if (i == k){
					// prepare mask so, that the peak is not at or near the origin
					for (int y = 0; y < 8; y++) { // 8 is an arbitrary value
						for (int x = 0; x < 8; x++) {
							mask(y,x) = 0;
						}
					}
				}
				double max;
				cv::Point maxLoc;
				cv::minMaxLoc(p, NULL, &max, NULL, &maxLoc, cv::Mat()) ; //mask);

				if (max < correlTh) 
				{
					continue;
				}

				if (Execution::verbosity >= 3){
					std::cerr << "real max " << max << " between subimg " << i << " " << " and " << k << " with maxLoc: " 
						<< maxLoc.x << " " << maxLoc.y << std::endl;
				}

				// get location in subimages
				cv::Size s;
				cv::Point f1Loc; 
				cv::Point f2Loc; 
				subimages[i].locateROI(s, f1Loc);
				subimages[k].locateROI(s, f2Loc);

				// correct offset -> (dx',dy')
				/*
				   if (maxLoc.x > p.cols / 2)
				   {
				   maxLoc.x -= p.cols;
				   }
				   if (maxLoc.y > p.rows / 2)
				   {
				   maxLoc.y -= p.rows;
				   }
				   */

				if (Execution::verbosity >= 3){
					std::cerr << "corrected max " << max << " between subimg " << i << " " << " and " << k << " with maxLoc: " 
						<< maxLoc.x << " " << maxLoc.y << std::endl;
				}
				f1Loc *= waveletMult;
				f2Loc *= waveletMult;
				maxLoc *= waveletMult;

				// dx,dy
				int dx, dy;
				// TODO: is that correct??
				//maxLoc.x = maxLoc.x % p.cols;
				//maxLoc.y = maxLoc.y % p.rows;

				// try 4 directions
				for (int l = 0; l < 4; l++) {
					//				if (f1Loc.x <= f2Loc.x)
					//					dx = f2Loc.x - f1Loc.x + maxLoc.x;
					//				else 
					//					dx = f2Loc.x + maxLoc.x;
					//				if (f1Loc.y <= f2Loc.y)
					//					dy = f2Loc.y - f1Loc.y + maxLoc.y;
					//				else 
					//					dy = f2Loc.y + maxLoc.y;

					if (l == 0){
						dx = p.cols - maxLoc.x + (f2Loc.x - f1Loc.x);
						dy = p.rows - maxLoc.y + (f2Loc.y - f1Loc.y);
					}
					else if (l == 1){
						dx = p.cols - maxLoc.x + (f2Loc.x - f1Loc.x);
						dy = maxLoc.y + (f2Loc.y - f1Loc.y);
					}
					else if (l == 2){
						dx = maxLoc.x + (f2Loc.x - f1Loc.x);
						dy = p.rows - maxLoc.y + (f2Loc.y - f1Loc.y);
					}
					else if (l == 3){
						dx = maxLoc.x + (f2Loc.x - f1Loc.x);
						dy = maxLoc.y + (f2Loc.y - f1Loc.y);
					}
				
					if (Execution::verbosity >= 3){
						std::cerr << "f1Loc: " << f1Loc.x << " " << f1Loc.y << std::endl;
						std::cerr << "f2Loc: " << f2Loc.x << " " << f2Loc.y << std::endl;
						std::cerr << "dx dy " << dx << " " << dy << std::endl;

						std::cerr << "borders: " << subimages[i].cols * waveletMult + f1Loc.x << " " << subimages[i].rows * waveletMult + f1Loc.y << std::endl;
					}

					int count = 0;
					for (int y = f1Loc.y; y < f1Loc.y + subimages[i].rows * waveletMult; y++)
					{
						for (int x = f1Loc.x; x < f1Loc.y + subimages[i].cols * waveletMult; x++)
						{
							if ( ((x + dx) < 0) || ((x + dx) >= origimg.cols) 
									|| ((y + dy) < 0) || ((y + dy) >= origimg.rows) )
								continue;
							if (count < 5 && Execution::verbosity >=3){
								std::cerr << "mark px\n";
								count++;
							}	
							double pxl = theImg(y,x);
							double pxlCopy = theImg(y + dy, x + dx);

							double sub = fabs(pxl - pxlCopy);

							if (sub < pixelDiff) // check pixel difference threshold
							{
								(matrices[l])(y,x)[c] = 255;
								(matrices[l])(y + dy, x + dx)[c] = 255;
							//	(matrix)(y,x)[c] = 255;
						//		(matrix)(y + dy, x + dx)[c] = 255;
							}
						}
					}

					if (Execution::verbosity >= 3){
						std::stringstream ss;
						ss << "matrix_" << i << k << c<< l << ".png";
						cv::imwrite(ss.str(), matrices[l]);
					}
				} // end 4 cases	
			
				int sum[4];
				int maxi = 0;
				int maxInd = 0;
				for (int l = 0; l < 4; l++){
					cv::Scalar sc = cv::sum(matrices[l]);
					sum[l] = sc[0] + sc[1] + sc[2];
					if (sum[l] > maxi){
						maxi = sum[l];
						maxInd = l;
					}
				}
				if (sum[maxInd] > 0){
					for (int y=0; y < img.rows * waveletMult; y++){
						for (int x=0; x < img.cols * waveletMult; x++){			
							for (int c = 0; c < 3; c++){
								if ( (matrix(y,x)[c] + (matrices[maxInd](y,x)[c])) < 256)
									matrix(y,x)[c] += (matrices[maxInd])(y,x)[c];
							}
						}
					}
				}
			
			} // end for all channels 
		}
	}

	if (Execution::verbosity >= 3)
		cv::imwrite("matrix.png", matrix);
}
