/*
 * Original work: Copyright (c) 2017, Thibaud Ehret <ehret.thibaud@gmail.com>
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
 * @file filters.cpp
 * @brief File containing the different filters used by the algorithm
 *
 * @author Thibaud Ehret <ehret.thibaud@gmail.com>
 **/


#include "filters.h"
#include <fstream>
#include <iostream>

int computeCircularDomain(std::vector<bool>& domain, int radius)
{
	int width = (2*radius+1);
	domain.resize(width*width, false);

	int radius2 = radius*radius;
	int sizeDomain = 0;

	for(int x = 0; x < width; ++x)
	for(int y = 0; y < width; ++y)
	{
		// Test if the position is in the disk define by the radius
		if((x-radius)*(x-radius)+(y-radius)*(y-radius) <= radius2)
		{
			domain[x+width*y] = true;
			sizeDomain++;
		}
	}
	return sizeDomain;
}

void medianFilter(std::vector<int>& dispMap, ImageSize imSize, int radius, bool xAxis)
{
	std::vector<int> backup(dispMap);

	std::vector<bool> domain;
	computeCircularDomain(domain, radius);

	int widthDomain = 2*radius+1;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(NTHREAD)
#endif
	for(unsigned id = 0; id < dispMap.size(); ++id)
	{
		int px = id % imSize.width;
		int py = id / imSize.width;

		std::vector<int> values;
		// Add all elements inside the disk of radius "radius" to a list
		for(int x = std::max(px-radius,0), dx=x-px+radius; x < std::min(px+radius+1, (int)imSize.width); ++x, ++dx) 
		for(int y = std::max(py-radius,0), dy=y-py+radius; y < std::min(py+radius+1, (int)imSize.height); ++y, ++dy) 
		{
			if(domain[dx+widthDomain*dy])
			{
				values.push_back(backup[x+imSize.width*y]);
			}
		}
		// Compute the median of this list
		std::nth_element(values.begin(), values.begin() + values.size()/2, values.end());
		dispMap[id] =  values[values.size()/2];

		// check if the new displacement is not illegal, if it is replace it by the closest legal one
		if(xAxis)
		{
			if((px + dispMap[id]) < 0)
				dispMap[id] = -px;
			if((px + dispMap[id]) >= imSize.width)
				dispMap[id] = imSize.width-1-px;
		}
		else // yAxis
		{
			if((py + dispMap[id]) < 0)
				dispMap[id] = -py;
			if((py + dispMap[id]) >= imSize.height)
				dispMap[id] = imSize.height-1-py;
		}
	}
}

int computeConnectedComponent(std::vector<bool>& binary, ImageSize imSize, std::vector<std::vector<unsigned> >& coords, std::vector<unsigned>& size)
{
	std::vector<bool> visited(binary.size(), false);

	int nb = 0;
	std::stack<unsigned> stack;

	for(unsigned id = 0; id < visited.size(); ++id)
	{
		if(!visited[id] && binary[id])
		{
			nb++;
			int count = 1;
			visited[id] = true;

			int idx = coords.size();
			coords.push_back(std::vector<unsigned>());
			coords[idx].push_back(id);

			stack.push(id);
			while(!stack.empty())
			{
				unsigned index = stack.top();
				stack.pop();

				unsigned a,b;
				a = index % imSize.width;
				b = index / imSize.width;

				// Define the list of the neighbors used to compute the component
				std::vector<std::pair<unsigned,unsigned>> neighbors{{a+1,b},{a-1,b},{a,b+1},{a,b-1}};

				for(auto n : neighbors) 
				{
					// Add the possible neighbors to the current component
					if((n.first >= 0) && (n.second >= 0) && (n.first < imSize.width) && (n.second < imSize.height) && !visited[n.first+imSize.width * n.second] && binary[n.first+imSize.width * n.second]) 
					{
						stack.push(n.first + imSize.width*n.second);
						coords[idx].push_back(n.first + imSize.width*n.second);
						count++;
						visited[n.first + imSize.width*n.second] = true;
					}
				}
			}
			size.push_back(count);
		}
	}
	return nb;
}

void sizeFilter(std::vector<bool>& detectionMask, ImageSize imSize, int minCompSize)
{
	std::vector<std::vector<unsigned> > compCoords;
	std::vector<unsigned> compSize;
	// Compute the list of all connected component with their size
	int nbComp = computeConnectedComponent(detectionMask, imSize, compCoords, compSize);

	for(int comp = 0; comp < nbComp; ++comp)
	{
		// If the component is too small, remove it
		if(compSize[comp] < minCompSize)
		{
			for(int pixel = 0; pixel < compCoords[comp].size(); ++pixel)
				detectionMask[compCoords[comp][pixel]] = false;
		}
	}
}

void errorDetectionFilter(std::vector<bool>& detectionMask, const std::vector<int>& dispX, const std::vector<int>& dispY, 
						ImageSize imSize, int radius, float th_e, int sp, std::vector<float>& visual)
{
	std::vector<bool> domain;
	int sizeDomain = computeCircularDomain(domain, radius);

	int widthDomain = 2*radius+1;

	// Load P (homogeneous coordinates of the points inside the disk of radius "radius")
	std::vector<float> P(3*sizeDomain, 0);
	int k = 0;
	for(int x = -radius, dx=0; x < radius+1; ++x, ++dx) 
	for(int y = -radius, dy=0; y < radius+1; ++y, ++dy) 
	{
		if(domain[dx+widthDomain*dy])
		{
			P[3*k] = x;
			P[1 + 3*k] = y;
			P[2 + 3*k] = 1;
			++k;
		}
	}

	// Compute H
	std::vector<float> H(sizeDomain*sizeDomain, 0);
	std::vector<float> temp(3*3, 0);
	productMatrix(temp, P, P, 3, 3, sizeDomain, false, true, true);
	inverseMatrix(temp, 3);
	std::vector<float> temp2(3*sizeDomain, 0);
	productMatrix(temp2, P, temp, sizeDomain, 3, 3, true, false, true);
	productMatrix(H, temp2, P, sizeDomain, sizeDomain, 3, false, false, true);

	// Compute I-H in H
	for(int i = 0; i < H.size(); ++i)
	{
		if((i%sizeDomain) == (i/sizeDomain))
			H[i] = 1 - H[i];
		else
			H[i] = -H[i];
	}
	
	std::vector<float> errors(dispX.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(NTHREAD)
#endif
	for(unsigned id = 0; id < dispX.size(); ++id)
	{
		int px = id % imSize.width;
		int py = id / imSize.width;

		errors[id] = 0.;

		std::vector<float> deltaX(sizeDomain);
		std::vector<float> deltaY(sizeDomain);
		int k2 = 0;

		// Extract the set of displacement on which the error will be computed
		for(int x2 = px-radius, dx2=0, k2=0; x2 < px+radius+1; ++x2, ++dx2) 
		for(int y2 = py-radius, dy2=0; y2 < py+radius+1; ++y2, ++dy2) 
		{
			if((x2 >= 0) && (y2 >= 0) && (x2 < (int)imSize.width) && (y2 < (int)imSize.height) && domain[dx2+widthDomain*dy2])
			{
				deltaX[k2] = dispX[x2+imSize.width*y2];
				deltaY[k2] = dispY[x2+imSize.width*y2];
				++k2;
			}
			else if(domain[dx2+widthDomain*dy2])
			{
				deltaX[k2] = 0.;
				deltaY[k2] = 0.;
				++k2;
			}
		}

		// Compute the error
		std::vector<float> outputX(sizeDomain);
		std::vector<float> outputY(sizeDomain);
		productMatrix(outputX, H, deltaX, sizeDomain, 1, sizeDomain, false, false, true);
		productMatrix(outputY, H, deltaY, sizeDomain, 1, sizeDomain, false, false, true);
		for(int p = 0; p < sizeDomain; ++p)
		{
			errors[id] += (outputX[p] * outputX[p]);
			errors[id] += (outputY[p] * outputY[p]);
		}
	}

	// Save the errors only for visualization purpose
	for(int i = 0; i < imSize.wh; ++i)
		visual[i] = std::max(std::log(errors[i]), -50.f);

	detectionMask.resize(errors.size());
	// Threshold the error map
	for(int i = 0; i < errors.size(); ++i)
		detectionMask[i] = (errors[i] < th_e);


	// Remove the boundary
	// Bottom
	for(int x = 0; x < imSize.width; ++x)
	for(int y = 0; y < sp; ++y)
		detectionMask[x+y*imSize.width] = false;
	// Top
	for(int x = 0; x < imSize.width; ++x)
	for(int y = imSize.height - sp; y < imSize.height; ++y)
		detectionMask[x+y*imSize.width] = false;
	// Left
	for(int x = 0; x < sp; ++x)
	for(int y = 0; y < imSize.height; ++y)
		detectionMask[x+y*imSize.width] = false;
	// Right
	for(int x = imSize.width - sp; x < imSize.width; ++x)
	for(int y = 0; y < imSize.height; ++y)
		detectionMask[x+y*imSize.width] = false;
}

void dilationFilter(std::vector<bool>& detectionMask, ImageSize imSize, int radius)
{
	std::vector<bool> backup(detectionMask);
	std::vector<bool> domain;
	computeCircularDomain(domain, radius);

	int widthDomain = 2*radius+1;

	for(unsigned id = 0; id < detectionMask.size(); ++id)
	{
		int px = id % imSize.width;
		int py = id / imSize.width;

		// if a pixel contains a detection at a distance smaller than radius then set the current pixel to a detection
		for(int x = std::max(px-radius,0), dx = x-px+radius; x < std::min(px+radius+1, (int)imSize.width); ++x, ++dx) 
		{
			if(detectionMask[id])
				break;

			for(int y = std::max(py-radius,0), dy = y-py+radius; y < std::min(py+radius+1, (int)imSize.height); ++y, ++dy) 
			{
				if(domain[dx+widthDomain*dy] && backup[x + imSize.width*y])
				{
					detectionMask[id] = true;
					break;
				}
			}
		}
	}
}

void minDispFilter(std::vector<bool>& detectionMask, const std::vector<int>& dispX, const std::vector<int>& dispY, float th_d)
{
	// if a forgery is detected without a large enough displacement, remove it
	for(int i = 0; i < dispX.size(); ++i)
		detectionMask[i] = detectionMask[i] && ((dispX[i] * dispX[i] + dispY[i] * dispY[i]) >= th_d);
}

void symmetrizationFilter(std::vector<bool>& detectionMask, const std::vector<int>& dispX, const std::vector<int>& dispY, ImageSize imSize)
{
	std::vector<bool> backup(detectionMask);

	for(int id = 0; id < dispX.size(); ++id)
	{
		// Symmetrize the detections
		detectionMask[id + dispX[id] + dispY[id]*imSize.width] = (backup[id] || backup[id + dispX[id] + dispY[id]*imSize.width]);
	}
}

void labelDetectionMask( const char* outname, std::vector<float>& detectionMask, \
						std::vector<int>& matchIds, ImageSize imSize, ImageSize visualSize, int radius){



	std::vector<int> detectionMapIds(matchIds.size(),0);
	std::vector<bool> domain;
	computeCircularDomain(domain, radius-1);
	int matchidI;
	int widthDomain = 2*radius+1;
	int max_label = 0;
	int label;
	int imgx, imgy;
	int vx, vy;

	for(int id=0; id<matchIds.size();id++){
		matchidI=matchIds[id];
		if(detectionMask[id] && matchidI && !detectionMapIds[id])
		{
			// If this the object was copied more than once, it will have a lot of others copied regions
			//
			if(detectionMapIds[matchidI]){
				label = detectionMapIds[matchidI]; 
			}
			else{
			label =++max_label;
			}

			fillObjectWithLabel(detectionMask,detectionMapIds, visualSize , id, label, 2);
			fillObjectWithLabel(detectionMask,detectionMapIds, visualSize , matchidI, label, 2);

		}
	}

	std::vector<int> writeMap(imSize.wh ,0);
	for(int i=0; i<visualSize.wh; i++){
			vx = i % visualSize.width;
			vy = i / visualSize.width;
			writeMap[vx + vy*imSize.width] = detectionMapIds[vx + vy*visualSize.width];
	}


	std::ofstream outFile;
	outFile.open(outname, std::ios::binary);
	for (int i = 0; i < writeMap.size(); i++)
		outFile.write(reinterpret_cast<const char *>(&writeMap[i]), sizeof(int));
	outFile.close();

}

/*
	Insert the same label inside a conncected componnet of an image
*/
void fillObjectWithLabel(std::vector<float>& detectionMask, std::vector<int>& detectionMapIds, \
 ImageSize imSize, int position, int label, int radius )
{
	std::vector<bool> visited(detectionMask.size(), false);

	int id ;
	int px;
	int py;
	std::stack<int> stack;
	stack.push(position);

	std::vector<bool> domain;
	computeCircularDomain(domain, radius-1);
	int widthDomain = 2*radius+1;

	while (!stack.empty())
	{
		id = stack.top();
		stack.pop();
		if(!visited[id] && detectionMask[id])
		{
			visited[id] = true;
			detectionMapIds[id] = label;

			px = id % imSize.width;
			py = id / imSize.width;

			for(int x = std::max(px-radius,0), dx = x-px+radius; x < std::min(px+radius+1, (int)imSize.width); ++x, ++dx) 
			for(int y = std::max(py-radius,0), dy = y-py+radius; y < std::min(py+radius+1, (int)imSize.height); ++y, ++dy) 
				if(domain[dx+widthDomain*dy] )
					if(detectionMask[x+imSize.width*y] && !visited[x+imSize.width*y])
						stack.push(x+imSize.width*y);
		}
	}
}