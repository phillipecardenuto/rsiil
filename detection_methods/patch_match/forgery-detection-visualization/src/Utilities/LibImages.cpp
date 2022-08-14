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
 * @file LibImages.cpp
 * @brief Usefull functions on images
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include "LibImages.h"
#include "Utilities.h"

#include <unistd.h> // getpid
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <math.h>

extern "C" {
#include "iio.h"
}

using namespace std;

/**
 * @brief Load image, check the number of channels.
 *
 * @param p_name : name of the image to read;
 * @param o_im : vector which will contain the image : R, G and B concatenated;
 * @param o_imSize : will contain the size of the image;
 * @param p_verbose : if true, print some informations.
 *
 * @return EXIT_SUCCESS if the image has been loaded, EXIT_FAILURE otherwise.
 **/
int loadImage(
	const char* p_name
,	std::vector<float> &o_im
,	ImageSize &o_imSize
,	const bool p_verbose
){
	//! read input image
	if (p_verbose) cout << endl << "Read input image...";

	float *imTmp = NULL;
	//size_t w, h, c;
	//imTmp = read_png_f32(p_name, &w, &h, &c);
	int w, h, c;
	imTmp =  iio_read_image_float_vec(p_name, &w, &h, &c);


	// Use iio to load either png or tiff
	// Because iio use a different order than the original library used (io_png), the resulting load is shuffled
	// FIXME It is a temporary hack so that nothing break while the transition of library
	
	float* finIm = (float*) malloc(w*c*h*sizeof(float));
	for(int ch = 0, i = 0; ch < c; ++ch)
	for(int y = 0; y < h; ++y)
	for(int x = 0; x < w; ++x, ++i)
		finIm[i] = imTmp[ch + x * c + y * c * w];
	free(imTmp);
	imTmp = finIm;

	if (!imTmp) {
		cout << "error :: " << p_name << " not found or not a correct png image" << endl;
		return EXIT_FAILURE;
	}

	if (p_verbose) cout << "done." << endl;

	//! test if image is really a color image and exclude the alpha channel
	if (c > 2) {
		unsigned k = 0;
		while (k < w * h && imTmp[k] == imTmp[w * h + k] && imTmp[k] == imTmp[2 * w * h + k])
			k++;

		c = (k == w * h ? 1 : 3);
	}

	//! Some image informations
	if (p_verbose) {
		cout << "image size :" << endl;
		cout << " - width          = " << w << endl;
		cout << " - height         = " << h << endl;
		cout << " - nb of channels = " << c << endl;
	}

	//! Initializations
	o_imSize.width      = w;
	o_imSize.height     = h;
	o_imSize.nChannels  = c;
	o_imSize.wh         = w * h;
	o_imSize.whc        = w * h * c;
	o_im.resize(w * h * c);
	for (unsigned k = 0; k < w * h * c; k++)
		o_im[k] = imTmp[k];

	free(imTmp);

	return EXIT_SUCCESS;
}

/**
 * @brief write image.
 *
 * @param p_name : path+name+extension of the image;
 * @param i_im : vector which contains the image;
 * @param p_imSize : size of the image;
 * @param p_min, p_max : range of data (usually [0, 255]).
 *
 * @return EXIT_SUCCESS if the image has been saved, EXIT_FAILURE otherwise
 **/
int saveImage(
    const char* p_name
,   std::vector<float> const& i_im
,   const ImageSize &p_imSize
,   const float p_min
,   const float p_max
){
    //! Allocate Memory
    float* imTmp = new float[p_imSize.whc];

    unsigned c = p_imSize.nChannels;
    unsigned w = p_imSize.width;
    unsigned h = p_imSize.height;

    //! Check for boundary problems
    //for (unsigned k = 0; k < p_imSize.whc; k++) {
    //    imTmp[k] = clip(i_im[k], p_min, p_max);
    //}
    for (unsigned ch = 0, k = 0; ch < c; ch++)
    for (unsigned y = 0; y < h; y++) {
    for (unsigned x = 0; x < w; x++, k++)
        //imTmp[ch + x * c + y * c * w] = clip(i_im[k], p_min, p_max);
        imTmp[ch + x * c + y * c * w] = i_im[k];
    }

    //if (write_png_f32(p_name, imTmp, p_imSize.width, p_imSize.height, p_imSize.nChannels) != 0) {
    //    cout << "... failed to save png image " << p_name << endl;
    //    return EXIT_FAILURE;
    //}

    iio_save_image_float_vec(p_name, imTmp, w, h, c);
    //! Free Memory
    delete[] imTmp;

    return EXIT_SUCCESS;
}

/**
 * @brief Transform the color space of an image, from RGB to YUV, or vice-versa.
 *
 * @param io_im: image on which the transform will be applied;
 * @param p_imSize: size of io_im;
 * @param p_isForward: if true, go from RGB to YUV, otherwise go from YUV to RGB.
 *
 * @return none.
 **/
void transformColorSpace(
	std::vector<float> &io_im
,	const ImageSize p_imSize
,	const bool p_isForward
){
	//! If the image as only one channel, do nothing
	if (p_imSize.nChannels == 1) return;

	//! Initialization
	const unsigned width  = p_imSize.width;
	const unsigned height = p_imSize.height;
	const unsigned chnls  = p_imSize.nChannels;
	const unsigned wh     = width * height;
	vector<float> imTmp(wh * chnls);

	//! RGB to YUV
	if (p_isForward) {
		if (chnls == 3) {
			const unsigned red   = 0;
			const unsigned green = wh;
			const unsigned blue  = wh * 2;
			const float a = 1.f / sqrtf(3.f);
			const float b = 1.f / sqrtf(2.f);
			const float c = 2.f * a * sqrtf(2.f);

			for (unsigned k = 0; k < wh; k++) {
				//! Y channel
				imTmp[k + red  ] = a * (io_im[k + red] + io_im[k + green] + io_im[k + blue]);

				//! U channel
				imTmp[k + green] = b * (io_im[k + red] - io_im[k + blue]);

				//! V channel
				imTmp[k + blue ] = c * (0.25f * io_im[k + red ] - 0.5f * io_im[k + green]
				                      + 0.25f * io_im[k + blue]);
			}
		}
		else { //! chnls == 4
			const unsigned Gr = 0;
			const unsigned R  = wh;
			const unsigned B  = wh * 2;
			const unsigned Gb = wh * 3;
			const float a = 0.5f;
			const float b = 1.f / sqrtf(2.f);

			for (unsigned k = 0; k < wh; k++) {
				imTmp[k + Gr] = a * ( io_im[k + Gr] + io_im[k + R ] +
				                      io_im[k + B ] + io_im[k + Gb]);
				imTmp[k + R ] = b * ( io_im[k + R ] - io_im[k + B ]);
				imTmp[k + B ] = a * (-io_im[k + Gr] + io_im[k + R ] +
				                      io_im[k + B ] - io_im[k + Gb]);
				imTmp[k + Gb] = b * (-io_im[k + Gr] + io_im[k + Gb]);
			}
		}
	}
	//! YUV to RGB
	else {
		if (chnls == 3) {
			const unsigned red   = 0;
			const unsigned green = wh;
			const unsigned blue  = wh * 2;
			const float a = 1.f / sqrtf(3.f);
			const float b = 1.f / sqrtf(2.f);
			const float c = a / b;

			for (unsigned k = 0; k < wh; k++) {
				//! R channel
				imTmp[k + red  ] = a * io_im[k + red] + b * io_im[k + green]
				                               + c * 0.5f * io_im[k + blue];
				//! G channel
				imTmp[k + green] = a * io_im[k + red] - c * io_im[k + blue];

				//! B channel
				imTmp[k + blue ] = a * io_im[k + red] - b * io_im[k + green]
				                               + c * 0.5f * io_im[k + blue];
			}
		}
		else {	//! chnls == 4
			const unsigned Gr = 0;
			const unsigned R  = wh;
			const unsigned B  = wh * 2;
			const unsigned Gb = wh * 3;
			const float a = 0.5f;
			const float b = 1.f / sqrtf(2.f);
			for (unsigned k = 0; k < wh; k++) {
				imTmp[k + Gr] = a * io_im[k + Gr] - a * io_im[k + B] - b * io_im[k + Gb];
				imTmp[k + R ] = a * io_im[k + Gr] + b * io_im[k + R] + a * io_im[k + B];
				imTmp[k + B ] = a * io_im[k + Gr] - b * io_im[k + R] + a * io_im[k + B];
				imTmp[k + Gb] = a * io_im[k + Gr] - a * io_im[k + B] + b * io_im[k + Gb];
			}
		}
	}

	io_im = imTmp;
}
