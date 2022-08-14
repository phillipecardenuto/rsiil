/*
 * Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
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
 * @file utilities.cpp
 * @brief Utilities functions.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include "Utilities.h"

#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

/**
 * @brief Convenient function to use the sort function provided by the vector library.
 **/
bool comparaisonFirst(
	const pair<float, unsigned> &i_pair1
,	const pair<float, unsigned> &i_pair2
){
	return i_pair1.first < i_pair2.first;
}

/**
 * @brief Convenient function to use the sort function provided by the vector library. Used to ocmpute the inverse order
 **/
bool comparaisonInverseFirst(
	const pair<float, unsigned> &i_pair1
,	const pair<float, unsigned> &i_pair2
){
	return i_pair1.first > i_pair2.first;
}

/**
 * @brief Clip a value between min and max
 *
 * @param i_value: value to clip;
 * @param i_min: minimum value;
 * @param i_max: maximum value.
 *
 * @return value clipped between [min, max].
 **/
float clip(
	const float i_value
,	const float i_min
,	const float i_max
){
	return (i_value < i_min ? i_min : (i_value > i_max ? i_max : i_value));
}

/**
 * @brief Obtain and substract the baricenter of io_group3d.
 *
 * @param io_group3d(p_rows x p_cols) : data to center;
 * @param o_baricenter(p_cols): will contain the baricenter of io_group3d;
 * @param p_rows, p_cols: size of io_group3d.
 *
 * @return none.
 **/
void centerData(
	std::vector<float> &io_group3d
,	std::vector<float> &o_baricenter
,	const unsigned p_rows
,	const unsigned p_cols
){
	const float inv = 1.f / (float) p_rows;
	for (unsigned j = 0; j < p_cols; j++) {
		float sum = 0.f;
		for (unsigned i = 0; i < p_rows; i++) {
			sum += io_group3d[j * p_rows + i];
		}

		o_baricenter[j] = sum * inv;

		for (unsigned i = 0; i < p_rows; i++) {
			io_group3d[j * p_rows + i] -= o_baricenter[j];
		}
	}
}

/**
 * @brief Substract baricenter to io_group3d.
 *
 * @param io_group3d(p_rows x p_cols) : data to center;
 * @param o_baricenter(p_cols): baricenter to remove;
 * @param p_rows, p_cols: size of io_group3d.
 *
 * @return none.
 **/
void centerDataWith(
	std::vector<float> &io_group3d
,	std::vector<float> &o_baricenter
,	const unsigned p_rows
,	const unsigned p_cols
){
	for (unsigned j = 0; j < p_cols; j++) {
		for (unsigned i = 0; i < p_rows; i++) {
			io_group3d[j * p_rows + i] -= o_baricenter[j];
		}
	}
}

/**
 * @brief Obtain and substract the baricenter of io_group3d.
 *
 * @param io_group3d(p_rows x p_cols) : data to center;
 * @param o_baricenter(p_cols): will contain the baricenter of io_group3d;
 * @param p_rows, p_cols: size of io_group3d.
 * @param weights: weights for the baricenter.
 *
 * @return none.
 **/
void centerWeightedData(
	std::vector<float> &io_group3d
,	std::vector<float> &o_baricenter
,	const unsigned p_rows
,	const unsigned p_cols
,	std::vector<double> &weights
){
	double inv = 0.f;
	for(int i = 0; i < weights.size(); ++i)
		inv += weights[i];
	if(inv > 0)
		inv = 1.f / inv;
	else
	{
		for(int i = 0; i < weights.size(); ++i)
			weights[i] = 1;
		inv = 1.f / weights.size();
	}
	for (unsigned j = 0; j < p_cols; j++) {
		double sum = 0.f;
		for (unsigned i = 0; i < p_rows; i++) {
			sum += weights[i] * io_group3d[j * p_rows + i];
		}

		o_baricenter[j] = sum * inv;

		for (unsigned i = 0; i < p_rows; i++) {
			io_group3d[j * p_rows + i] -= o_baricenter[j];
		}
	}
}

/**
 * @brief Compute the average standard deviation of a set of patches.
 *
 * @param i_Set(p_sP, p_nSimP): set of patches;
 * @param p_sP : size of a patch;
 * @param p_nSimP: number of patches in the set;
 * @param p_nChannels: number of channels of the image.
 *
 * @return the average standard deviation of the set
 **/
float computeStdDeviation(
	std::vector<float> const& i_Set
,	const unsigned p_sP
,	const unsigned p_nSimP
,	const unsigned p_nChannels
){
	float sigma = 0.f;

	for (unsigned c = 0; c < p_nChannels; c++) {
		//! Initialization
		float mean = 0.f;
		float std = 0.f;

		//! Compute the sum and the square sum
		for (unsigned n = 0; n < p_nSimP; n++) {
			for (unsigned k = 0; k < p_sP; k++) {
				const float value = i_Set[k + c * p_sP + n * p_sP * p_nChannels];
				mean += value;
				std  += value * value;
			}
		}

		//! Sample standard deviation (Bessel's correction)
		sigma += (std - mean * mean / (float) (p_sP * p_nSimP)) / (float) (p_sP * p_nSimP - 1);
	}

	return sigma / (float) p_nChannels;
}

/**
 * @brief Determine a and b such that : n = a * b, with a and b as greatest as possible
 *
 * @param i_n : number to decompose;
 * @param o_a : will contain a;
 * @param o_b : will contain b.
 *
 * @return none.
 **/
void determineFactor(
    const unsigned i_n
,   unsigned &o_a
,   unsigned &o_b
){
    if (i_n == 1)
	 {
        o_a = 1;
        o_b = 1;
        return;
    }

    o_b = 2;
    while (i_n % o_b > 0) o_b++;
    o_a = i_n / o_b;

    if (o_b > o_a)
	 {
        o_a = o_b;
        o_b = i_n / o_a;
    }
}

/**
 * @brief Write PSNR and RMSE in a .txt for both basic and denoised images.
 *
 * @param p_pathName: path name of the file;
 * @param p_sigma: value of the noise;
 * @param p_psnr: value of the PSNR of the denoised image;
 * @param p_rmse: value of the RMSE of the denoised image;
 * @param p_truncateFile: if true, erase the file when open it. Otherwise
 *        write at the end of the file;
 * @param p_app: in order to specify the image.
 *
 * @return EXIT_FAILURE if the file can't be opened.
 **/
int writingMeasures(
    const char* p_pathName
,   const float p_sigma
,   const float p_psnr
,   const float p_rmse
,   const float p_grps
,   const float p_time
,   const float p_cons
,   const bool  p_truncateFile
,   const char* p_app
){
    //! Open the file
    ofstream file;
    file.open(p_pathName, ios::out | ios::app);

    //! Check if the file is open
    if (!file) {
        return EXIT_FAILURE;
    }

    //! Write measures in the file
    //if (p_truncateFile) {
    //    file << "************" << endl;
    //    file << "-sigma = " << p_sigma << endl;
    //    file << "************" << endl;
    //}
    //file << "-PSNR" << p_app << " = " << p_psnr << endl;
    file << "RMSE = " << p_rmse << endl;
    //file << "-GRPS" << p_app << " = " << p_grps << endl;
    //file << "-TIME" << p_app << " = " << p_time << endl;
    //file << "-CONS" << p_app << " = " << p_cons << endl;
    cout << endl;

    //! Close the file
    file.close();

    return EXIT_SUCCESS;
}
