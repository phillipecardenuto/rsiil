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
 * @file patchmatch.cpp
 * @brief Apply PatchMatch algorithm on a feature manager
 *
 * @author Thibaud Ehret <ehret.thibaud@gmail.com>
 **/

#include "patchmatch.h"

static int rand_gen(int i)
{
	return std::rand() % i;
}

Patchmatch::Patchmatch(FeatManager& fm_, int th_)
{
	fm = &fm_;
	imSize = fm->sz();
	th = th_;
}

int Patchmatch::initializeRandom()
{
	matches.resize(imSize.wh);
	distances.resize(imSize.wh);

	std::default_random_engine generator;

	std::uniform_int_distribution<unsigned> match_distribution(0,imSize.wh-1);

	int matchId;
	int offset;
	int dif1, dif2;

	// For each patch compute a valid random match as well as the distance the patch and its match
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(NTHREAD)
#endif
	for (unsigned id = 0; id < imSize.wh; ++id)
	{
		do {
			matchId = match_distribution(generator);
			offset = (dif1 = ((id%imSize.width) - (matchId%imSize.width)))*dif1 + (dif2 = (id/imSize.width - matchId/imSize.width))*dif2;
		} while(offset < th);
		matches[id] = matchId;
		distances[id] = fm->distance(id, matchId);
	}

	return 0;
}


void Patchmatch::propagate()
{
	for(int id = 0; id < imSize.wh; ++id)
	{
		// Propagate from the left
		propToId(id, 0);
		// Propagate from the top
		propToId(id, 1);
		// Propagate from the top-left
		propToId(id, 2);
		// Propagate from the top-right
		propToId(id, 3);


		// Propagate from the left
		propToId1(id, 0);
		// Propagate from the top
		propToId1(id, 1);
		// Propagate from the top-left
		propToId1(id, 2);
		// Propagate from the top-right
		propToId1(id, 3);
	}

	// Symmetrize all the matches if necessary
	for(int id = 0; id < imSize.wh; ++id)
	{
		symmetrize(id);
	}
}

void Patchmatch::backpropagate()
{
	for(int id = 0; id < imSize.wh; ++id)
	{
		// Propagate from the right
		propToId(id, 4);
		// Propagate from the bottom
		propToId(id, 5);
		// Propagate from the bottom-right
		propToId(id, 6);
		// Propagate from the bottom-left
		propToId(id, 7);

		// Propagate from the right
		propToId1(id, 4);
		// Propagate from the bottom
		propToId1(id, 5);
		// Propagate from the bottom-right
		propToId1(id, 6);
		// Propagate from the bottom-left
		propToId1(id, 7);
	}

	// Symmetrize all the matches if necessary
	for(int id = 0; id < imSize.wh; ++id)
	{
		symmetrize(id);
	}
}

void Patchmatch::propagateNtimes(int N)
{
	for(int i = 0; i < N; ++i)
	{
		if(i % 2 == 0)
			propagate();
		else
			backpropagate();

		randomSearch();
	}	
}

void Patchmatch::extract(std::vector<int>& dispX, std::vector<int>& dispY)
{
	dispX.resize(imSize.wh);
	dispY.resize(imSize.wh);

	for(int ind = 0; ind < imSize.wh; ++ind)
	{
		int x = ind % imSize.width;
		int y = ind / imSize.width;

		int nx = matches[ind] % imSize.width;
		int ny = matches[ind] / imSize.width;

		dispX[ind] = nx - x;	
		dispY[ind] = ny - y;	
	}
}

// Zero order propagation (basic propagation from patchmatch)
inline void Patchmatch::propToId(unsigned id, int direction)
{
	unsigned propId = fm->getInverseNeighboringFeat(id, direction);
	unsigned newId = fm->getNeighboringFeat(matches[propId], direction);
	int offset, dif1, dif2;
	if(newId != id)
	{
		float dist = fm->distance(id, newId);
		offset = (dif1 = ((id%imSize.width) - (newId%imSize.width)))*dif1 + (dif2 = (id/imSize.width - newId/imSize.width))*dif2;
		if(dist < distances[id] && offset >= th)
		{
			matches[id] = newId;
			distances[id] = dist;
		}
	}
}

// First order propagation (specific from this paper)
inline void Patchmatch::propToId1(unsigned id, int direction)
{
	int idX = id % imSize.width;
	int idY= id / imSize.width;

	unsigned propId  = fm->getInverseNeighboringFeat(id, direction);
	int x = propId % imSize.width;
	int y = propId / imSize.width;
	unsigned propId1 = fm->getInverseNeighboringFeat(propId, direction);
	int x1 = propId1 % imSize.width;
	int y1 = propId1 / imSize.width;

	int px  = matches[propId] % imSize.width;
	int py  = matches[propId] / imSize.width;
	int px1 = matches[propId1] % imSize.width;
	int py1 = matches[propId1] / imSize.width;

	int newDeltaX = 2*(px-x) - (px1-x1);
	int newDeltaY = 2*(py-y) - (py1-y1);

	int offset, dif1, dif2;
	if(idX+newDeltaX >= 0 && idX+newDeltaX < imSize.width &&  idY+newDeltaY >= 0 && idY+newDeltaY < imSize.height)
	{
		int newId = (idX+newDeltaX) + (idY+newDeltaY)*imSize.width;
		if(newId != id)
		{
			float dist = fm->distance(id, newId);
			offset = newDeltaX*newDeltaX + newDeltaY*newDeltaY;
			if(dist < distances[id] && offset >= th)
			{
				matches[id] = newId;
				distances[id] = dist;
			}
		}
	}
}

void Patchmatch::randomSearch()
{
	std::default_random_engine generator;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(NTHREAD)
#endif
	for(int id = 0; id < imSize.wh; ++id)
	{
		int init_size_search_region = std::max(imSize.width, imSize.height); 

		int px = matches[id] % imSize.width;
		int py = matches[id] / imSize.width;

		int dif1, dif2;
		for (int imSize_sr = init_size_search_region; imSize_sr >= 1; imSize_sr /= 2) 
		{
			// Sampling window 
			unsigned x_min = std::max((int)px-imSize_sr, 0), x_max = std::min(px+imSize_sr+1, (int)imSize.width-1);
			unsigned y_min = std::max((int)py-imSize_sr, 0), y_max = std::min(py+imSize_sr+1, (int)imSize.height-1);
			std::uniform_int_distribution<unsigned> x_distribution(x_min,x_max);
			std::uniform_int_distribution<unsigned> y_distribution(y_min,y_max);

			// Compute the new candidate
			unsigned x_feat = x_distribution(generator);
			unsigned y_feat = y_distribution(generator);

			unsigned newId = y_feat*imSize.width + x_feat;

			if(newId != id)
			{
				// Verify if the new candidate is better than the previous match and if so update the displacement map
				float dist = fm->distance(id, newId);
				int offset = (dif1 = (x_feat - (id % imSize.width)))*dif1 + (dif2 = (y_feat - (id / imSize.width)))*dif2;
				if(dist < distances[id] && offset >= th)
				{
					matches[id] = newId;
					distances[id] = dist;
				}
			}
		}
	}
}

inline void Patchmatch::symmetrize(unsigned id)
{
	if(distances[id] < distances[matches[id]])
		matches[matches[id]] = id;
}
