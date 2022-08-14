/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "dwtfeat.h"
#include "wavelet.h"
Dwtfeat :: Dwtfeat(struct ft_cfg config,
				   const BlockHandling &blocks,
				   int numThreads)
	:	Feature(config, blocks, numThreads)
{
	allocateFeature(block_size*block_size);
	if (config.decomposeLvl <= 0){
		std::cout << "WARNING: chose DWTFEAT with decomposelevel <= 0\n";
	}
}

void Dwtfeat :: computeOne(const cv::Mat & cur, int row)
{
	std::vector<cv::Mat_<double> > dest;
	if (config.decomposeLvl <= 0){
		Feature::computeOne(cur, row);
		return;
	}

	if (config.modified)
		Wavelet::decomposeHaar(cur, dest, config.decomposeLvl, true);
	else 
		Wavelet::decomposeHaar(cur, dest, config.decomposeLvl, false);
	
	int f = 0;
	// need reverse order LL2 --> HH1, so reverse iterator
	std::vector<cv::Mat_<double> >::reverse_iterator rit;
	for (rit = dest.rbegin(); rit != dest.rend(); ++rit)
	{
		cv::MatIterator_<double> it = rit->begin();
		for (; it != rit->end() ; ++it)
		{
			// copy subbands to the featurevector
			feature(row, f++) = *it;
		}
	}

}
