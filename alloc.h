//
// $Id: alloc.h,v 1.1 2006/06/27 18:41:28 nakayama Exp $
//
/* ---- memory allocation ---- */


#include <stdlib.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

int *alloc_1d_int(int n1);
void free_1d_int(int *i);
double *alloc_1d_double(int n1);
void free_1d_double(double *d);
int **alloc_2d_int(int n1, int n2);
void free_2d_int(int **ii);
double **alloc_2d_double(int n1, int n2);
void free_2d_double(double **dd);
int ***alloc_3d_int(int n1, int n2, int n3);
void free_3d_int(int ***iii);
double ***alloc_3d_double(int n1, int n2, int n3);
void free_3d_double(double ***ddd);

#ifdef __cplusplus
}
#endif
