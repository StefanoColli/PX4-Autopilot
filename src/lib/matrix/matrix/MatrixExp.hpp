/**
 * @file MatrixExp.hpp
 *
 * Implementation of the matrix exponential function, adapted from Eigen implementation (https://eigen.tuxfamily.org).
 * Details of the algorithm can be found in: Nicholas J. Higham, "The scaling and squaring method for the
 * 	matrix exponential revisited," SIAM J. Matrix Anal. Applic., 26:1179–1193, 2005.
 *
 * @author Stefano Colli <stefano2.colli@mail.polimi.it>
 */

#pragma once

#include "math.hpp"

namespace matrix
{

/**
 * @brief Element wise matrix scaling: each element is divided by 2^times
 *
 * @tparam Type
 * @tparam M
 * @param A square matrix to scale
 * @param times scale factor (each element is divided by 2^times)
 */
template<typename Type, size_t M>
void scale(SquareMatrix<Type, M> &A, int times)
{
	for (size_t i = 0; i < M; i++)
	{
		for (size_t j = 0; j < M; j++)
		{
    			A(i,j) = std::ldexp(A(i,j), -times);
		}
	}
}

/**
 * @brief Compute U and V for the (3,3)-Padè approximant to the exponential.
 * After exit, (V+U)(V-U)^{-1} is the Padè approximant of e^A around A = 0
 *
 * @tparam Type
 * @tparam M
 * @param A input matrix
 * @param U output matrix U
 * @param V output matrix V
 */
template<typename Type, size_t M>
void matrix_exp_pade3(const SquareMatrix<Type, M>& A, SquareMatrix<Type, M> &U, SquareMatrix<Type, M> &V)
{
	Type b[] = {120., 60., 12., 1.};
	SquareMatrix<Type, M> I;
	I.setIdentity();
	const SquareMatrix<Type, M> A2 = A * A;
	const SquareMatrix<Type, M> tmp = b[3] * A2 + b[1] * I;
	U = A * tmp;
	V = b[2] * A2 + b[0] * I;
}

/**
 * @brief Compute U and V for the (5,5)-Padè approximant to the exponential.
 * After exit, (V+U)(V-U)^{-1} is the Padè approximant of e^A around A = 0
 *
 * @tparam Type
 * @tparam M
 * @param A input matrix
 * @param U output matrix U
 * @param V output matrix V
 */
template<typename Type, size_t M>
void matrix_exp_pade5(const SquareMatrix<Type, M>& A, SquareMatrix<Type, M> &U,  SquareMatrix<Type, M> &V)
{
  Type b[] = {30240., 15120., 3360., 420., 30., 1.};
  SquareMatrix<Type, M> I;
  I.setIdentity();
  const SquareMatrix<Type, M> A2 = A * A;
  const SquareMatrix<Type, M> A4 = A2 * A2;
  const SquareMatrix<Type, M> tmp = b[5] * A4 + b[3] * A2 + b[1] * I;
  U = A * tmp;
  V = b[4] * A4 + b[2] * A2 + b[0] * I;
}

/**
 * @brief Compute U and V for the (7,7)-Padè approximant to the exponential.
 * After exit, (V+U)(V-U)^{-1} is the Padè approximant of e^A around A = 0
 *
 * @tparam Type
 * @tparam M
 * @param A input matrix
 * @param U output matrix U
 * @param V output matrix V
 */
template<typename Type, size_t M>
void matrix_exp_pade7(const SquareMatrix<Type, M>& A, SquareMatrix<Type, M> &U,  SquareMatrix<Type, M> &V)
{
  Type b[] = {17297280., 8648640., 1995840., 277200., 25200., 1512., 56., 1.};
  SquareMatrix<Type, M> I;
  I.setIdentity();
  const SquareMatrix<Type, M> A2 = A * A;
  const SquareMatrix<Type, M> A4 = A2 * A2;
  const SquareMatrix<Type, M> A6 = A4 * A2;
  const SquareMatrix<Type, M> tmp = b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * I;
  U = A * tmp;
  V = b[6] * A6 + b[4] * A4 + b[2] * A2 + b[0] * I;

}

/**
 * @brief Compute U and V for the (9,9)-Padè approximant to the exponential.
 * After exit, (V+U)(V-U)^{-1} is the Padè approximant of e^A around A = 0
 *
 * @tparam Type
 * @tparam M
 * @param A input matrix
 * @param U output matrix U
 * @param V output matrix V
 */
template<typename Type, size_t M>
void matrix_exp_pade9(const SquareMatrix<Type, M>& A, SquareMatrix<Type, M> &U,  SquareMatrix<Type, M> &V)
{
  Type b[] = {17643225600., 8821612800., 2075673600., 302702400., 30270240.,
                          	2162160., 110880., 3960., 90., 1.};
  SquareMatrix<Type, M> I;
  I.setIdentity();
  const SquareMatrix<Type, M> A2 = A * A;
  const SquareMatrix<Type, M> A4 = A2 * A2;
  const SquareMatrix<Type, M> A6 = A4 * A2;
  const SquareMatrix<Type, M> A8 = A6 * A2;
  const SquareMatrix<Type, M> tmp = b[9] * A8 + b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * I;
  U = A * tmp;
  V = b[8] * A8 + b[6] * A6 + b[4] * A4 + b[2] * A2 + b[0] * I;
}

/**
 * @brief Compute U and V for the (13,13)-Padè approximant to the exponential.
 * After exit, (V+U)(V-U)^{-1} is the Padè approximant of e^A around A = 0
 *
 * @tparam Type
 * @tparam M
 * @param A input matrix
 * @param U output matrix U
 * @param V output matrix V
 */
template<typename Type, size_t M>
void matrix_exp_pade13(const SquareMatrix<Type, M>& A, SquareMatrix<Type, M> &U,  SquareMatrix<Type, M> &V)
{
  Type b[] = {64764752532480000., 32382376266240000., 7771770303897600.,
                          	1187353796428800., 129060195264000., 10559470521600., 670442572800.,
                          	33522128640., 1323241920., 40840800., 960960., 16380., 182., 1.};
  SquareMatrix<Type, M> I;
  I.setIdentity();
  const SquareMatrix<Type, M> A2 = A * A;
  const SquareMatrix<Type, M> A4 = A2 * A2;
  const SquareMatrix<Type, M> A6 = A4 * A2;
  V = b[13] * A6 + b[11] * A4 + b[9] * A2; // used for temporary storage
  SquareMatrix<Type, M> tmp = A6 * V;
  tmp += b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * I;
  U = A * tmp;
  tmp = b[12] * A6 + b[10] * A4 + b[8] * A2;
  V = A6 * tmp;
  V += b[6] * A6 + b[4] * A4 + b[2] * A2 + b[0] * I;
}

/**
 * @brief Compute U and V matrices needed for the exponential computation: (V+U)(V-U)^{-1} is the Padè approximant
 *
 * @tparam Type
 * @tparam M
 * @param A input square matrix
 * @param U output matrix U
 * @param V output matrix V
 * @return int number of squarings done in order to scale input matrix
 */
template<typename Type, size_t M>
int compute_UV(const SquareMatrix<Type, M>& A, SquareMatrix<Type, M>& U, SquareMatrix<Type, M>& V)
{
	int squarings = 0;
	const Type l1norm = A.l1norm();

	SquareMatrix<Type, M> A_local = A; 	//local copy


	if (l1norm < 4.1968497232266989671e-003L)
	{
      		matrix_exp_pade3(A_local, U, V);
    	}
	else if (l1norm < 1.1848116734693823091e-001L)
	{
      		matrix_exp_pade5(A_local, U, V);
    	}
	else if (l1norm < 5.5170388480686700274e-001L)
	{
      		matrix_exp_pade7(A_local, U, V);
    	}
	else if (l1norm < 1.3759868875587845383e+000L)
	{
      		matrix_exp_pade9(A_local, U, V);
    	}
	else
	{
      		const Type maxnorm = 4.0246098906697353063L;
      		std::frexp(l1norm / maxnorm, &squarings);
      		if (squarings < 0) squarings = 0;
		scale(A_local, squarings);
      		matrix_exp_pade13(A_local, U, V);
    	}
	return squarings;
}

/**
 * @brief Compute the exponential e^A of square matrix A
 *
 * @tparam Type
 * @tparam M
 * @param A input matrix
 * @return SquareMatrix<Type, M> matrix exponential e^A
 */
template<typename Type, size_t M>
SquareMatrix<Type, M> exp(const SquareMatrix<Type, M>& A)
{
  	SquareMatrix<Type, M> U, V;
  	int squarings = compute_UV<Type, M>(A, U, V); // (V+U)(V-U)^{-1} is the Padè approximant

  	const SquareMatrix<Type, M> num = V + U;
  	const SquareMatrix<Type, M> denom = V - U;

	SquareMatrix<Type, M> result = inv(denom) * num;
  	for (int i=0; i<squarings; i++)
	{
    		result *= result;   // undo scaling by repeated squaring
	}
	return result;
}

}

