/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include <iostream>
#include <cmath>
#include <stdexcept>

#include "dct.h"
#include "cv.h"

//#define BLK_SIZE 16
//#define PI 3.14159265

//const double Dct::quantizationMatrix[8][8] = {
//								{ 4, 4, 6, 11, 24, 24, 24, 24 },
//								{ 4, 5, 6, 16, 24, 24, 24, 24 },
//								{ 6, 6, 14, 24, 24, 24, 24, 24 },
//								{ 11, 16, 24, 24, 24, 24, 24, 24 },
//								{ 24, 24, 24, 24, 24, 24, 24, 24 },
//								{ 24, 24, 24, 24, 24, 24, 24, 24 },
//								{ 24, 24, 24, 24, 24, 24, 24, 24 },
//								{ 24, 24, 24, 24, 24, 24, 24, 24 } };

Dct :: Dct(struct ft_cfg _config, const BlockHandling &b, int threadNum)
	: Feature(_config, b, threadNum), PI(acos(-1.0))
{
	quant8 = (cv::Mat_<double>(8,8) <<	4, 4, 6, 11, 24, 24, 24, 24,
										4, 5, 6, 16, 24, 24, 24, 24,
										6, 6, 14, 24, 24, 24, 24, 24,
										11, 16, 24, 24, 24, 24, 24, 24,
										24, 24, 24, 24, 24, 24, 24, 24,
										24, 24, 24, 24, 24, 24, 24, 24,
										24, 24, 24, 24, 24, 24, 24, 24,
										24, 24, 24, 24, 24, 24, 24, 24 );
	q = config.qf;
	if (q == 0.0) {
		throw std::runtime_error("Dct::computeFeatureVec:: q == 0.0 is a wrong qualifier");
	}
	if (block_size >16){
		throw std::runtime_error("DCT::computeFeatureVec:: wrong blocksize for DCT method - has to be 16");
	}

	//quantisierungsmatrix
	buildQuantizationMatrix();

	//Vorberechnen der Cosinusterme,cos_term[i][j] entspricht dabei: cos((pi/(2*16)) * i*(2j +1))
//	PI = acos(-1.0);
	cos_terms = cv::Mat_<double>(block_size, block_size);
	for(int i = 0;i < (int)block_size;i++)
	{
		for(int j = 0; j < (int)block_size;j++)
		{
			cos_terms[i][j]	= cos((PI/(2*block_size)) *(double)i * (2*(double)j+1));
		}
	}

	allocateFeature(block_size*block_size);
	feature.setTo(0.0);
}

void Dct :: computeOne(const cv::Mat & cur, int row)
{	
	double sqrt_116 = 0.25;
	double sqrt_216 = sqrt(1.0/8.0);
	//double pi_sixteen = PI/16;

	// TODO: this could be done easier perhaps
	for(int k = 0; k < (int)block_size; k++)
	{
		for(int l = 0; l < (int)block_size; l++)
		{
			for(int m = 0; m < (int)block_size; m++)
			{
				for(int n = 0; n < (int)block_size; n++)
				{
					feature(row, k*(int)block_size+l) += cur.at<uchar>(n, m)
													* cos_terms[k][m]
													* cos_terms[l][n];
				}
			}
		}
	}

	//quantisieren
	for(int k = 0; k < (int)block_size; k++)
	{
		for(int l = 0; l < (int)block_size; l++)
		{
			feature(row, k*(int)block_size + l) *=  (k==0) ? sqrt_116 : sqrt_216;
			feature(row, k*(int)block_size + l) *=  (l==0) ? sqrt_116 : sqrt_216;
			feature(row, k*(int)block_size + l) /= q;
			feature(row, k*(int)block_size + l) /= quant16[k][l];
			feature(row, k*(int)block_size + l) = round(feature(row, k*block_size + l));
		}
	}
}

void Dct :: buildQuantizationMatrix(void)
{

	quant16 = cv::Mat_<double>((int)block_size,(int)block_size);
	// 8x8 quadrants of the 16x16 matrix
	cv::Mat lo = quant16(cv::Range(0,8), cv::Range(0,8));
	cv::Mat ro = quant16(cv::Range(0,8), cv::Range(8,16));
	cv::Mat lu = quant16(cv::Range(8,16), cv::Range(0,8));
	cv::Mat ru = quant16(cv::Range(8,16), cv::Range(8,16));

	// now set them appropriate:
	cv::Mat_<double> tmp = 2.5*quant8;
	tmp.copyTo(lo);
	lo.at<double>(0.0) = 2.0*quant8(0,0);

	ro.setTo(2.5*quant8(0,7));
	lu.setTo(2.5*quant8(7,0));
	ru.setTo(2.5*quant8(7,7));

/*
	for(int i = 0 ; i<8;i++)
	{
		for(int j = 0; j< 8;j++)
		{
			quant16[i][j] = 2.5*quantizationMatrix[i][j];
		}
	}

	quant16[0][0] = 2.0 * quantizationMatrix[0][0];
		
	for(int i = 8 ; i<16;i++)
	{
		for(int j = 0; j< 8;j++)
		{
			quant16[i][j] = 2.5*quantizationMatrix[0][7];
		}
	}
	for(int i = 8 ; i<16;i++)
	{
		for(int j = 8; j< 16;j++)
		{
			quant16[i][j] = 2.5*quantizationMatrix[7][7];
		}
	}
	for(int i = 0;i<8;i++)
	{
		for(int j = 8;j<16;j++)
		{
			quant16[i][j] = 2.5 * quantizationMatrix[7][0];
		}
	}
*/
}


