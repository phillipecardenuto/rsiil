/*
 * Copyright (c) 2017, Thibaud Ehret <ehret.thibaud@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef FILTERS_H_INCLUDED
#define FILTERS_H_INCLUDED

#include <vector>
#include <stack>
#include <algorithm>
#include "LibImages.h"
#include "parameters.h"
#include "LibMatrix.h"
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * @brief Filter the image using the median of the element inside the given the radius.
 *
 * @param dispMap: the displacements on which the median filtering is applied 
 * @param imSize: size of the displacement map
 * @param radius: radius of the disk on which the median is computed
 * @param xAxis: specify whether the displacement is along the x axis or the y axis. Necessary to avoid illegal displacement after the median filter
 **/
void medianFilter(std::vector<int>& dispMap, ImageSize imSize, int radius, bool xAxis);

/**
 * @brief Filter the image using the median of the element inside the given the radius.
 *
 * @param detectionMask: detection mask which needs to be filtered
 * @param imSize: size of the detection mask
 * @param minCompSize: minimum accepted size of a connected component in the detection mask
 **/
void sizeFilter(std::vector<bool>& detectionMask, ImageSize imSize, int minCompSize);

/**
 * @brief Compute the error map
 * @param detectionMask: used to save the detection mask
 * @param dispX: displacement map along the x axis representing the PatchMatch matching
 * @param dispY: displacement map along the y axis representing the PatchMatch matching
 * @param imSize: size of the detection mask
 * @param radius: radius of the disk on which the error is computed
 * @param th_e: threshold on the error for the detection of forgeries
 * @param sp: size of the patch on which the features where computed so as to remove the corresponding border
 * @param visual: used to save the errors for visual purpose only
 **/
void errorDetectionFilter(std::vector<bool>& detectionMask, const std::vector<int>& dispX, const std::vector<int>& dispY, ImageSize imSize, int radius, float th_e, int sp, std::vector<float>& visual);

/**
 * @brief Filter the image using dilation process
 *
 * @param detectionMask: detection mask which needs to be filtered
 * @param imSize: size of the detection mask
 * @param radius: radius of the disk on which the dilation is computed
 **/
void dilationFilter(std::vector<bool>& detectionMask, ImageSize imSize, int radius);

/**
 * @brief Filter the image using dilation process
 *
 * @param detectionMask: detection mask which needs to be filtered
 * @param dispX: displacement map along the x axis representing the PatchMatch matching
 * @param dispY: displacement map along the y axis representing the PatchMatch matching
 * @param th_d: minimum accepted distance between two detected matches
 **/
void minDispFilter(std::vector<bool>& detectionMask, const std::vector<int>& dispX, const std::vector<int>& dispY, float th_d);

/**
 * @brief Filter the image using a symmetrization process
 *
 * @param detectionMask: detection mask which needs to be filtered
 * @param imSize: size of the detection mask (and the displacement maps)
 * @param dispX: displacement map along the x axis representing the PatchMatch matching
 * @param dispY: displacement map along the y axis representing the PatchMatch matching
 **/
void symmetrizationFilter(std::vector<bool>& detectionMask, const std::vector<int>& dispX, const std::vector<int>& dispY, ImageSize imSize);


/**
 * @brief Label each detection region inserting the same id to its matched region
 *
 * @param detectionMask: detection mask containing boolean map of the detection region
 * @param matchIds: matchIds contain the Ids of each matched pixel
 **/
void labelDetectionMask(const char* outname, std::vector<float>& detectionMask,std::vector<int>& matchIds, ImageSize imSize, ImageSize visualSize, int radius);

/**
 * @ brief Fill a object region (connected component) with a label
 * 
 * @param detectionMask: detection Mask containing the detected regions
 * @param detectionMapIds: Map in which the label will be propagate to
 * @param label: label of the region
 * @param radius: radius of the considerate neighborhood during the label propagation method
 **/
void fillObjectWithLabel(std::vector<float>& detectionMask, std::vector<int>& detectionMapIds, ImageSize imSize, int position, int label, int radius );

#endif // FILTERS_H_INCLUDED

