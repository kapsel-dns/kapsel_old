/*!
  \file alloc.h
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief Memory allocation routines (header file)
 */

#include <stdlib.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
  \brief Allocate memory for 1d int array
 */
int *alloc_1d_int(int n1);
/*!
  \brief Free memory for 1d int array
 */
void free_1d_int(int *i);


/*!
  \brief Allocate memory for 1d double array
 */
double *alloc_1d_double(int n1);
/*!
  \brief Free memory for 1d double array
 */
void free_1d_double(double *d);

/*!
  \brief Allocate memory for 2d int array
  \details Allocates linear array of size n1*n2, position of element
  (i,j) is im = i*n2 + j
 */
int **alloc_2d_int(int n1, int n2);
/*!
  \brief Free memory for 2d int array
 */
void free_2d_int(int **ii);

/*!
  \brief Allocate memory for 2d double array
  \details Allocates linear array of size n1*n2, position of element
  (i,j) is im = i*n2 + j
 */
double **alloc_2d_double(int n1, int n2);
/*!
  \brief Free memory for 2d double array
 */
void free_2d_double(double **dd);

/*!
  \brief Allocate memory for 3d int array
  \details Allocates linear array of size n1*n2*n3, position of
  element (i,j,k) is im = i*(n2 * n3) + j*n3 + k
 */
int ***alloc_3d_int(int n1, int n2, int n3);
/*!
  \brief Free memory for 3d int array
 */
void free_3d_int(int ***iii);

/*!
  \brief Allocate memory for 3d double array
  \details Allocates linear array of size n1*n2*n3, position of 
  element (i,j,k) is im = i*(n2*n3) + j*n3 + k
 */
double ***alloc_3d_double(int n1, int n2, int n3);
/*!
  \brief Free memory for 3d double array
 */
void free_3d_double(double ***ddd);

#ifdef __cplusplus
}
#endif
