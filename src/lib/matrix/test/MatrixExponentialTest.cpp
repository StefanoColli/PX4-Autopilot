/****************************************************************************
 *
 *   Copyright (C) 2023 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

#include <gtest/gtest.h>
#include <matrix/math.hpp>

using namespace matrix;

TEST(MatrixExponentialTest, Exponential1)
{
	//checking values are computed using matlab expm() function

	// test A
	double data_a[9] = { 0, 1, 3,
			   0, 0, 2,
			   0, 0, 0};
	SquareMatrix<double, 3> A(data_a);

	double data_check_a[9] = { 1, 1, 4,
			  	   0, 1, 2,
			  	   0, 0, 1};
	SquareMatrix<double, 3> A_check(data_check_a);
	EXPECT_TRUE(isEqual(exp(A),A_check, 1e-3d));


	// test B

	double data_b[9] = {0.01, 0.02, 0.03,
			    0.04, 0.05, 0.06,
			    0.07, 0.08, 0.1
			};
	SquareMatrix<double, 3> B(data_b);

	double data_check_b[9] = {
		1.01158503f,  0.02190432f,  0.03238144f,
		0.04349195f,  1.05428524f,  0.06539627f,
		0.07576783f,  0.08708946f,  1.10894048f
	};
	SquareMatrix<double, 3> B_check(data_check_b);

	EXPECT_TRUE(isEqual(exp(B), B_check, 1e-3d));

	// test C
	double data_c[4] = { 1, 2,
		             3, 4};
	SquareMatrix<double, 2> C(data_c);

	double data_check_c[4] = { 51.9689561987050, 74.7365645670033,
				   112.104846850505, 164.073803049210};
	SquareMatrix<double, 2> C_check(data_check_c);

	EXPECT_TRUE(isEqual(exp(C),C_check, 1e-3d));

	// test D

	double data_d[4] = { -1, 0,
		              0, 3};
	SquareMatrix<double, 2> D(data_d);

	double data_check_d[4] = { 0.36787944117144233, 	0,
				   	0, 			20.085536923187664};
	SquareMatrix<double, 2> D_check(data_check_d);

	EXPECT_TRUE(isEqual(exp(D),D_check, 1e-3d));




}

