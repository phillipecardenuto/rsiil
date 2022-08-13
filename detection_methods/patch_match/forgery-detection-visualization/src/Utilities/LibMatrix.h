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
#ifndef LIB_MATRIX_H_INCLUDED
#define LIB_MATRIX_H_INCLUDED

#include <vector>
#include <string>

/**
 * @brief Invert (in place) a symmetric real matrix, V -> Inv(V).
 *			The input matrix V is symmetric (V[i,j] = V[j,i]).
 *
 * @param io_mat = array containing a symmetric input matrix. This is converted to the inverse
			matrix;
 * @param p_N = dimension of the system (dim(v)=n*n)
 *
 * @return	EXIT_SUCCESS -> normal exit
 *			EXIT_FAILURE -> input matrix not positive definite
 **/
int inverseMatrix(
	std::vector<float> &io_mat
,	const unsigned p_N
);

/**
 * @brief Transpose a real square matrix in place io_mat -> io_mat~.
 *
 * @param io_mat : pointer to an array of p_N by p_N input matrix io_mat. This is overloaded by
			the transpose of io_mat;
 * @param p_N : dimension (dim(io_mat) = p_N * p_N).
 *
 * @return none.
 **/
void transposeMatrix(
	std::vector<float> &io_mat
,	const unsigned p_N
);

/**
 * @brief Multiply two matrix A * B. It uses BLAS SGEMM.
 *
 * @param o_AB = array containing n by l product matrix at exit;
 * @param i_A = input array containing n by m matrix;
 * @param i_B = input array containing m by l matrix;
 * @param p_n, p_m, p_l = dimension parameters of arrays.
 * @param p_transA = true for transposing A.
 * @param p_transA = true for transposing B.
 * @param p_colMajor = true for if matrices should be read by columns.
 *
 * @return  none.
 **/
void productMatrix(
	std::vector<float> &o_AB
,	std::vector<float> const& i_A
,	std::vector<float> const& i_B
,	const unsigned p_n
,	const unsigned p_m
,	const unsigned p_l
,	const bool p_transA
,	const bool p_transB
,	const bool p_colMajor = true
,	unsigned lda = 0
,	unsigned ldb = 0
);


#endif // LIB_MATRIX_H_INCLUDED
