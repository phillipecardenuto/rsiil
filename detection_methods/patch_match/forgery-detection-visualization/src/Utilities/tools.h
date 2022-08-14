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

#ifndef VISTOOLS_H_INCLUDED
#define VISTOOLS_H_INCLUDED

/**
 * @brief Structure containing size informations of an image.
 *
 * @param width     : width of the image;
 * @param height    : height of the image;
 * @param nChannels : number of channels in the image;
 * @param wh        : equal to width * height. Provided for convenience;
 * @param whc       : equal to width * height * nChannels. Provided for convenience.
 **/
void rescale(std::vector<float>& im, ImageSize sz)
{
	for(unsigned c = 0; c < sz.nChannels; ++c)
	{
		// Compute the min and max
		float mm = im[c*sz.wh];
		float mM = im[c*sz.wh];
		for(unsigned x = 0; x < sz.width; ++x)
			for(unsigned y = 0; y < sz.height; ++y)
			{
				if(im[c*sz.wh + y*sz.width + x] < mm)
					mm = im[c*sz.wh + y*sz.width + x];
				if(im[c*sz.wh + y*sz.width + x] > mM)
					mM = im[c*sz.wh + y*sz.width + x];
			}
		// map the min to 0 and the max to 255.
		float coeff = 255. / (mM - mm);
		for(unsigned x = 0; x < sz.width; ++x)
			for(unsigned y = 0; y < sz.height; ++y)
			{
				im[c*sz.wh + y*sz.width + x] = ((im[c*sz.wh + y*sz.width + x] - mm) * coeff);
			}
	}
}


// This function and the following ones have been inspired by the code from http://www.ipol.im/pub/art/2013/26/
static void hsv_to_rgb_doubles(double *out, double *in)
{
	//assert_hsv(in);
	double r, g, b, h, s, v; r=g=b=h=s=v=0;
	h = in[0]; s = in[1]; v = in[2];
	if (s == 0)
		r = g = b = v;
	else {
		int H = fmod(floor(h/60),6);
		double p, q, t, f = h/60 - H;
		p = v * (1 - s);
		q = v * (1 - f*s);
		t = v * (1 - (1 - f)*s);
		switch (H) {
			case 6:
			case 0: r = v; g = t; b = p; break;
			case 1: r = q; g = v; b = p; break;
			case 2: r = p; g = v; b = t; break;
			case 3: r = p; g = q; b = v; break;
			case 4: r = t; g = p; b = v; break;
			case -1:
			case 5: r = v; g = p; b = q; break;
			default:
				fprintf(stderr, "H=%d\n", H);
		}
	}
	out[0] = r; out[1] = g; out[2] = b;
	//assert_rgb(out);
}

void view_displacement(std::vector<float>& image, ImageSize sz, std::vector<int>& px, std::vector<int>& py)
{
	for(int id = 0; id < sz.wh; ++id)
	{
		double rho = std::sqrt(px[id]*px[id] + py[id]*py[id]);

		if (rho > 1000000) {
			image[id] = 255; 
			image[id + sz.wh] = 255; 
			image[id + 2*sz.wh] = 255; 
			continue;
		}

		rho = std::min(rho, 500.)/500.;
		double theta = atan2(py[id], -px[id]);
		theta = (theta+M_PI)*(180/M_PI);
		theta = fmod(theta, 360);
		double hsv[3], rgb[3];
		hsv[0] = theta;
		hsv[1] = rho;
		hsv[2] = rho;
		hsv_to_rgb_doubles(rgb, hsv);
		image[id] = rgb[0]; 
		image[id + sz.wh] = rgb[1]; 
		image[id + 2*sz.wh] = rgb[2]; 
	}
}

#endif // VISTOOLS_H_INCLUDED
