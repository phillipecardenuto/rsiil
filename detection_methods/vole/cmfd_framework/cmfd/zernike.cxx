/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "zernike.h"
#include "execution.h"

#include <cmath>

Zernike :: Zernike(struct ft_cfg _config, const BlockHandling & b, int numThreads)
	:	Feature(_config, b, numThreads)
{
	feature_size = 0;
	for(int i=0; i <= config.degZernike; i++){
		feature_size += (i/2+1);
	}
	allocateFeature(feature_size);
	pi = std::acos(-1.0);

	// prepare coords for polar coords
	cv::Mat_<double> X(block_size, block_size);
	cv::Mat_<double> Y(block_size, block_size);
	for (int y = 0; y < block_size; y++){
		for (int x = 0; x < block_size; x++){
			X(y,x) = x - (block_size / 2);
			Y(y,x) = y - (block_size / 2);
		}
	}
	// computes polar coord in Radiens around the center
	cv::Mat_<double> rho; // magnitude
	cv::Mat_<double> theta; // angle
	cv::Mat r,t; // tmp-matrices as cartToPolar can't deal with templated version?!
	cv::cartToPolar(X, Y, r, t);
	rho = r;
	theta = t;

	// create mask with circle
	cv::Mat_<uchar> mask(block_size, block_size, static_cast<uchar>(0));
	cv::circle(mask, cv::Point(block_size/2, block_size/2),
				static_cast<int>(block_size/2), cv::Scalar(1,0), -1);

	zernikeMoments.resize(block_size*block_size);
	for (int y = 0; y <block_size; y++){
		for (int x = 0; x < block_size; x++){
			zernikeMoments[y*block_size + x] = cv::Mat_<cv::Vec2d>(config.degZernike+1, config.degZernike+1, 0.0);
			if (mask(y,x) == 0)
				continue;
			for (int n = 0; n <= config.degZernike; n++){
				for (int m = 0; m <= n; m++) {
					if ( ((n-m) % 2) != 0)
						continue;
					double poly = zernikePoly(n, m, rho(y,x) / (block_size/2)); // scale to unit-radius
					if (cvIsNaN(poly)){
						std::cerr << " n " << n << " m " << m << " rho: " << rho(y,x) << std::endl;
						throw std::runtime_error("Zernike::Zernike(): maeh is nan");
					}
					
					// express e^(j*m*theta) w. sin & cos terms
					// real part
					double real = poly * std::cos(m*theta(y,x));
					// complex part
					double comp = - poly * std::sin(m*theta(y,x));
					// magnitude
					zernikeMoments[y*block_size + x](n,m)[0] = real;
					zernikeMoments[y*block_size + x](n,m)[1] = comp;
				}
			}
		}
	}
}

unsigned long long  Zernike :: fak(int n)
{
	unsigned long long f = 1;
	for (int i = 2; i <= n; i++){
		f *= i;
	}
	return f;
}

double Zernike :: zernikePoly(int n, int m, double rho)
{
	double r = 0.0;
	for (int s = 0; s <= static_cast<int>((n - m) / 2); s++){
		r += (std::pow((double)-1,(double)s) * fak(n-s) * std::pow((double)rho, (double)(n-2*s)) ) /
				(fak(s) * fak((n+m)/2-s) * fak((n-m)/2-s));
	}
	return r;
}

void Zernike :: computeOne(const cv::Mat &curr, int row)
{
	if (Execution::verbosity >= 4){
		std::cerr << row << " ";
	}

	// create circular mask
	cv::Mat_<uchar> mask(curr.rows, curr.cols, static_cast<uchar>(0));
	cv::circle(mask, cv::Point(curr.cols/2, curr.rows/2),
			   static_cast<int>(curr.rows/2), cv::Scalar(1,0), -1);

	int x = 0;
	for (int n = 0; n <= config.degZernike; n++){
		for (int m = 0; m <= n; m++) {
			if ( ((n-m) % 2) != 0)
				continue;
			double moment_r = 0.0;
			double moment_i = 0.0;
			for (int y = 0; y < curr.rows; y++){
				for (int x = 0; x < curr.cols; x++){
					if (mask(y,x) == 0)
						continue;
					moment_r += static_cast<double>(curr.at<uchar>(y, x))
							* zernikeMoments.at(y*block_size + x)(n, m)[0];
					moment_i += static_cast<double>(curr.at<uchar>(y, x))
							* zernikeMoments[y*block_size + x](n, m)[1];
				}
			}
			// normalize
			moment_r *= (n+1) / pi;
			moment_i *= (n+1) / pi;
			// magnitude
			feature(row, x++) = static_cast<float>(std::sqrt(moment_r*moment_r + moment_i*moment_i));
		}
	}
}
