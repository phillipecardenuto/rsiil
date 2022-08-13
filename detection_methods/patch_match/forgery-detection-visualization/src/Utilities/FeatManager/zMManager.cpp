/*
 * Original work Copyright (c) 2017, Thibaud Ehret <ehret.thibaud@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file zMManager.cpp
 * @brief Feature manager for Zernike moments
 *
 * @author Thibaud Ehret <ehret.thibaud@gmail.com>
 **/

#include "zMManager.h"

#define V 4

// Compute a bicubic interpolation
inline float bicubicKernel(float x)
{
	x = std::abs(x);
	return ((x <= 1) ? 1.5*x*x*x - 2.5*x*x + 1 : ((x <= 2) ? -0.5*x*x*x + 2.5*x*x - 4*x + 2 : 0));
}

ZMManager::ZMManager(const std::vector<float>& image, ImageSize dsz, int sp, int order, int numRho, int numTheta, int zm, bool flip)
{
	descSize = dsz;
	computeZM(image, dsz, sp, order, numRho, numTheta, zm, features1);
	
	if(flip)
	{
		// If we want to compute the matching between the image and the flipped version

		// First the flipped image is computed
		std::vector<float> flipped_image(image.size());
		for(int x = 0; x < dsz.width; ++x)
		for(int y = 0; y < dsz.height; ++y)
			flipped_image[(dsz.width-x-1) + dsz.width*y] = image[x + dsz.width*y];
		
		// Then the features are computed on the flipped image
		std::vector<std::vector<float> > temp_features;
		computeZM(flipped_image, dsz, sp, order, numRho, numTheta, zm, temp_features);

		// Finally flip back these features so they match the position of the ones of the original image
		features2.resize(temp_features.size());
		for(int x = 0; x < dsz.width; ++x)
		for(int y = 0; y < dsz.height; ++y)
			features2[(dsz.width-x-1) + dsz.width*y] = temp_features[x + dsz.width*y];
	}
	else
	{
		// Otherwise we match the image and the features have already been computed
		features2 = features1;
	}
}

void ZMManager::computeZM(const std::vector<float>& image, ImageSize dsz, int sp, int order, int numRho, int numTheta, int zm, std::vector<std::vector<float> >& features)
{
	int widthPatch = 2*sp;

	std::vector<std::pair<int, int> > indexes;

	// Compute the list of possible pairs of coefficients for Zernike moments up to the order "order"
	for(int n = 0; n <= order; ++n)
	{
		for(int l = 0; l <= n; ++l)
		{
			if((n - l) % 2 == 0)
				indexes.push_back(std::make_pair(n,l));
		}
	}

	if(zm == 1)
	{
		printf("Using resampled Zernike Moments\n");
		// Start by initializing factorial
		std::vector<double> factorial(order+1);
		factorial[0] = 1.;
		for(int i = 1; i <= order; ++i)
			factorial[i] = factorial[i-1]*i;

		// compute the Zernike coefficients that does not depend on the position
		std::vector<std::vector<float> > coefs(indexes.size(), std::vector<float>(order));
		std::vector<float> weights(indexes.size());
		for(int i = 0; i < indexes.size(); ++i)
		{
			for(int s = 0; s <= (indexes[i].first-indexes[i].second)/2; ++s)
				coefs[i][s] = ((s % 2 == 0) ? 1 : -1) * factorial[indexes[i].first-s] / ((indexes[i].first - 2*s + 2) * factorial[s] * factorial[(indexes[i].first + indexes[i].second)/2 - s] * factorial[(indexes[i].first - indexes[i].second) / 2 - s]);
			weights[i] = (indexes[i].first + 1) / PI;
		}

		std::vector<std::vector<std::complex<float> > > convFilter(indexes.size(), std::vector<std::complex<float> >(widthPatch*widthPatch, std::complex<float>(0.,0.)));

		// sampling along rho + sampling along theta to create a convolution filter  for the current Zernike moment
		for(int k = 0; k < sp; ++k)
			for(int l = 0; l < V*(2*k+1); ++l)
			{
				// Compute the corresponding cartesian position
				float rx = k*std::cos(l * 2*PI/((2*k+1)*V));
				float ry = k*std::sin(l * 2*PI/((2*k+1)*V));

				for(int i = 0; i < indexes.size(); ++i)
				{
					// Compute the position specific coefficients
					std::complex<float> w(0,0);
					for(int s = 0; s <= (indexes[i].first-indexes[i].second)/2; ++s)
						w += coefs[i][s] * (pow((float)(k+1)/sp, indexes[i].first - 2*s + 2) - pow((float)k/sp, indexes[i].first - 2*s + 2));

					w *= (indexes[i].second == 0) ? 2*PI/((2*k+1)*V) : std::complex<float>(0,1) * (std::polar(1.f, -(float)indexes[i].second * (l+1) * 2*(float)PI/((2*l+1)*V)) - std::polar(1.f, -(float)indexes[i].second * l * 2*(float)PI/((2*k+1)*V))) / (float)indexes[i].second;

					// Do a bicubic interpolation for the "equivalent pixel" value
					float rx = k*std::cos(l * 2*(float)PI/((2*k+1)*V));
					float ry = k*std::sin(l * 2*(float)PI/((2*k+1)*V));
					int initIdx1 = sp + std::floor(rx);
					int initIdx2 = sp + std::floor(ry);

					for(int idx1 = std::max(initIdx1-1, 0); idx1 < std::min(initIdx1+2, widthPatch); ++idx1)
						for(int idx2 = std::max(initIdx2-1, 0); idx2 < std::min(initIdx2+2, widthPatch); ++idx2)
						{
							convFilter[i][idx1 + widthPatch*idx2] += bicubicKernel(rx - idx1 + sp) * bicubicKernel(ry - idx2 + sp) * w;
						}

				}
			}

		features.clear();
		features.resize(descSize.wh, std::vector<float>(indexes.size()));
		// Compute the convolution of the image and the Zernike moments to create the final feature map
		for(int n = 0; n < indexes.size(); ++n)
		{
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(NTHREAD)
#endif
			for(int id = 0; id < descSize.wh; ++id)
			{
				std::complex<float> temp(0.,0.);
				int px = id % descSize.width;
				int py = id / descSize.width;
				for(int x = std::max(px-sp,0), dx=x-px+sp; x < std::min(px+sp, (int)descSize.width); ++x, ++dx) 
					for(int y = std::max(py-sp,0), dy=y-py+sp; y < std::min(py+sp, (int)descSize.height); ++y, ++dy) 
						temp += (image[x + descSize.width*y] * convFilter[n][dx + dy*widthPatch]);

				features[id][n] = std::abs(weights[n] * temp);
			}
		}
	}
	else
	{
		printf("Using basic Zernike Moments\n");

		// Start by initializing factorial
		std::vector<double> factorial(order+1);

		factorial[0] = 1.;
		for(int i = 1; i <= order; ++i)
			factorial[i] = factorial[i-1]*i;

		// compute the Zernike coefficients that does not depend on the position
		std::vector<std::vector<float> > coefs(indexes.size(), std::vector<float>(order));
		std::vector<float> weights(indexes.size());
		for(int i = 0; i < indexes.size(); ++i)
		{
			for(int j = 0; j <= (indexes[i].first-indexes[i].second)/2; ++j)
				coefs[i][j] = ((j % 2 == 0) ? 1 : -1) * factorial[indexes[i].first-j] / (factorial[indexes[i].first] * factorial[(indexes[i].first + indexes[i].second)/2 - j] * factorial[(indexes[i].first - indexes[i].second) / 2 - j]);
			weights[i] = (indexes[i].first + 1) / PI;
		}

		std::vector<std::vector<std::complex<float> > > convFilter(indexes.size(), std::vector<std::complex<float> >(widthPatch*widthPatch, std::complex<float>(0.,0.)));

		for(int x = 0; x < 2*sp; ++x)
			for(int y = 0; y < 2*sp; ++y)
			{
				// Compute the equivalent polar coordinates
				float r = sqrt((x-sp)*(x-sp) + (y-sp)*(y-sp));
				float theta = atan2(y-sp, x-sp);
				for(int i = 0; i < indexes.size(); ++i)
				{
					// Compute the position specific coefficients
					std::complex<float> w(0,0);
					for(int j = 0; j <= (indexes[i].first-indexes[i].second)/2; ++j)
						w += coefs[i][j] * pow(r, indexes[i].first - 2*j);

					w *= std::polar(1.f, -indexes[i].second*theta);

					convFilter[i][x + widthPatch*y] += w;
				}
			}

		features.resize(descSize.wh, std::vector<float>(indexes.size()));
		// Compute the convolution of the image and the Zernike moments to create the final feature map
		for(int n = 0; n < indexes.size(); ++n)
		{
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(NTHREAD)
#endif
			for(int id = 0; id < descSize.wh; ++id)
			{
				std::complex<float> temp(0.,0.);
				int px = id % descSize.width;
				int py = id / descSize.width;
				for(int x = std::max(px-sp,0), dx=x-px+sp; x < std::min(px+sp, (int)descSize.width); ++x, ++dx) 
					for(int y = std::max(py-sp,0), dy=y-py+sp; y < std::min(py+sp, (int)descSize.height); ++y, ++dy) 
						temp += (image[x + descSize.width*y] * convFilter[n][dx + dy*widthPatch]);

				features[id][n] = std::abs(weights[n] * temp);
			}
		}
	}
}


float ZMManager::distance(unsigned id1, unsigned id2)
{
	float dist = 0.f, dif;
	for (unsigned c = 0; c < features1[0].size(); ++c)
		dist += (dif = features1[id1][c] - features2[id2][c]) * dif;
	return dist;
}

unsigned ZMManager::getNeighboringFeat(unsigned id, int direction)
{
	int width = descSize.width;
	int height = descSize.height;

	int x = id % width;
	int y = id / width;

	// Right
	if(direction == 0)
	{
		if(x+1 < width)
			return (x+1)+y*width;
		return id;
	}
	// Up
	else if(direction == 1)
	{
		if(y+1 < height)
			return x+(y+1)*width;
		return id;
	}
	// Up-Right
	else if(direction == 2)
	{
		if(y+1 < height && x+1 < width)
			return (x+1)+(y+1)*width;
		return id;
	}
	// Up-Left
	else if(direction == 3)
	{
		if(x-1 >= 0 && y+1 < height)
			return (x-1)+(y+1)*width;
		return id;
	}
	// Left
	else if(direction == 4)
	{
		if(x-1 >= 0)
			return (x-1)+y*width;
		return id;
	}
	// Down
	else if(direction == 5)
	{
		if(y-1 >= 0)
			return x+(y-1)*width;
		return id;
	}
	// Down-Left
	else if(direction == 6)
	{
		if(x-1 >= 0 && y-1 >= 0)
			return (x-1)+(y-1)*width;
		return id;
	}
	// Down-Right
	else if(direction == 7)
	{
		if(y-1 >= 0 && x+1 < width)
			return (x+1)+(y-1)*width;
		return id;
	}
	else
	{
		// This direction of propagation doesn't exist ...
		return 0;
	}
}


unsigned ZMManager::getInverseNeighboringFeat(unsigned id, int direction)
{

	int width = descSize.width;
	int height = descSize.height;

	int x = id % width;
	int y = id / width;

	// Right
	if(direction == 4)
	{
		if(x+1 < width)
			return (x+1)+y*width;
		return id;
	}
	// Up
	else if(direction == 5)
	{
		if(y+1 < height)
			return x+(y+1)*width;
		return id;
	}
	// Up-Right
	else if(direction == 6)
	{
		if(y+1 < height && x+1 < width)
			return (x+1)+(y+1)*width;
		return id;
	}
	// Up-Left
	else if(direction == 7)
	{
		if(x-1 >= 0 && y+1 < height)
			return (x-1)+(y+1)*width;
		return id;
	}
	// Left
	else if(direction == 0)
	{
		if(x-1 >= 0)
			return (x-1)+y*width;
		return id;
	}
	// Down
	else if(direction == 1)
	{
		if(y-1 >= 0)
			return x+(y-1)*width;
		return id;
	}
	// Down-Left
	else if(direction == 2)
	{
		if(x-1 >= 0 && y-1 >= 0)
			return (x-1)+(y-1)*width;
		return id;
	}
	// Down-Right
	else if(direction == 3)
	{
		if(y-1 >= 0 && x+1 < width)
			return (x+1)+(y-1)*width;
		return id;
	}
	else
	{
		// This direction of propagation doesn't exist ...
		return 0;
	}
}

ZMManager::~ZMManager()
{

}
