/*!
  \file matrix_diagonal.h
  \brief Jacobi Method to diagonalize real symmetric matrices
 */
#ifndef MATRIX_DIAGONAL_H
#define MATRIX_DIAGONAL_H

#include "alloc.h"
#include "macro.h"

/*!
  \brief Perform diagonalization of real symmetric matrix using
  Jacobi's iterative procedure
  \details \f[
  A \longrightarrow D = V^T A V
  \f]
  where D is diagonal and V is orthogonal. Results are returned in
  decrasing eigenvalue order and eigenvectors are normalized.
  H. Rutihauser, "The Jacobi Method for Real Symmetric
  Matrices", Numerische Mathematic, 9, 1-10, 1966
  \param[in] A_input Symmetric square matrix (only upper triangular
  elements used)
  \param[out] v elements of V matrix ( rows are eigenvectors of V)
  \param[out] d diagonal elements of matrix D (eigenvalues of A)
  \param[out] jrot number of jacobi rotations needed to achieve
  diagonalization 
  \param[in] n dimension of square matrix A
 */
void jacobi(double **A_input, double **V, double *d, int &rot, const int &n);

#endif
