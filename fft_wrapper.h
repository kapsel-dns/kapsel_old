//
// $Id: fft_wrapper.h,v 1.1 2006/06/27 18:41:28 nakayama Exp $
//
#ifndef FFT_WRAPPER_H
#define FFT_WRAPPER_H

#include <assert.h> 
#include <math.h>

#include "variable.h"
#include "input.h"
#include "alloc.h"
#include "macro.h"

#ifdef _OPENMP
#include <omp.h>
#include <mkl_dfti.h>
#include <complex.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern void rdft3d(int n1, int n2, int n3, int sign, double ***a, double *t, int *ip, double *w);
extern void rdft3dsort(int, int, int, int, double ***);

#ifdef __cplusplus
}
#endif

////////////////////////

extern int *ip;
extern double *w;
extern double *t;

extern double **ucp;
extern double *phi,**up,**u, *rhop;
extern double **work_v3, *work_v1;

extern int *KX_int, *KY_int, *KZ_int;
extern double *K2,*IK2;

extern int (*Calc_KX)( const int &i, const int &j, const int &k);
extern int (*Calc_KY)( const int &i, const int &j, const int &k);
extern int (*Calc_KZ)( const int &i, const int &j, const int &k);
extern void (*Truncate_two_third_rule)(double *a);

//// function prototype
void Init_fft(void);
void A_k2dxa_k(double *a, double *da);
void A_k2dya_k(double *a, double *da);
void A_k2dza_k(double *a, double *da);

void Omega_k2zeta_k(double **omega, double **zetak);
void Omega_k2zeta_k_OBL(double **omega, double **zetak);
void U_k2zeta_k(double **u, double **zeta, double uk_dc[DIM]);
void Zeta_k2u_k(double **zeta, double uk_dc[DIM], double **u);
void Zeta_k2omega_k_OBL(double **zeta, double **omega);
void U_k2zeta_k_OBL(double **u, double **zeta, double uk_dc[DIM]);
void Zeta_k2u_k_OBL(double **zeta, double uk_dc[DIM], double **u);
void Zeta_k2Strain_k(double **zeta, double *strain_k[QDIM]);
void U_k2divergence_k(double **u, double *div);
void U_k2rotation_k(double **u);

inline void U2u_oblique(double **uu) {
    
    int im;
    int im_ob;
    int im_ob_p;

    for (int i = 0; i <NX; i++) {
	for (int j = 0; j <NY; j++) {
	    for (int k = 0; k <NZ; k++) {
		im = (i*NY*NZ_) + (j*NZ_)+k;

		work_v3[0][im] = uu[0][im];
		work_v3[1][im] = uu[1][im];
		work_v3[2][im] = uu[2][im];
	    }
	}
    }
    
    for (int i = 0; i < NX; i++) {
	for (int j = 0; j < NY; j++) {

	    double sign = j - NY/2;
	    //int sign = j - NY/2;
	    if (!(sign == 0)) {
		sign = sign/fabs(sign);
	    }

	    int i_oblique = (int)(sign*degree_oblique*(j - NY/2))*sign;
	    double alpha = (degree_oblique*(j - NY/2) - i_oblique)*sign;
	    double beta  = 1. - alpha;

	    i_oblique      = (int) fmod(i + i_oblique + 4.*NX, NX);
	    int i_oblique_plus = (int) fmod(i_oblique + sign + 4.*NX, NX);

	    
	    for (int k = 0; k < NZ; k++) {
		im      = (i*NY*NZ_) + (j*NZ_) + k;
		im_ob   = (i_oblique*NY*NZ_) + (j*NZ_) + k;
		im_ob_p = (i_oblique_plus*NY*NZ_) + (j*NZ_) + k;

		uu[0][im] =
		    (beta*work_v3[0][im_ob] +
		     alpha*work_v3[0][im_ob_p])
		    - degree_oblique*(beta*work_v3[1][im_ob] +
				      alpha*work_v3[1][im_ob_p]);

		uu[1][im] =
		    beta*work_v3[1][im_ob] +
		    alpha*work_v3[1][im_ob_p];

		uu[2][im] =
		    beta*work_v3[2][im_ob] +
		    alpha*work_v3[2][im_ob_p];
	    }
	}
    }
    
}

inline void U_oblique2u(double **uu) {

    int im;
    int im_ob;
    int im_p;

    for (int i = 0; i <NX; i++) {
	for (int j = 0; j <NY; j++) {
	    for (int k = 0; k <NZ_; k++) {
		im = (i*NY*NZ_) + (j*NZ_)+k;

		work_v3[0][im] = uu[0][im];
		work_v3[1][im] = uu[1][im];
		work_v3[2][im] = uu[2][im];
	    }
	}
    }

    for (int i = 0; i < NX; i++) {
	for (int j = 0; j < NY; j++) {

	    double sign = j - NY/2;
	    //int sign = j - NY/2;
	    if (!(sign == 0)) {
		sign = sign/fabs(sign);
	    }

	    int i_oblique = (int)(sign*degree_oblique*(j - NY/2.))*sign + sign;
	    double alpha = (i_oblique - degree_oblique*(j - NY/2.))*sign;
	    double beta  = 1. - alpha;
	    
	    i_oblique      = (int) fmod(i + i_oblique + 4.*NX, NX);

	    int i_plus = (int) fmod(i + sign + 2*NX, NX);

	    for (int k = 0; k < NZ; k++) {
		im      = (i*NY*NZ_) + (j*NZ_) + k;
		im_ob   = (i_oblique*NY*NZ_) + (j*NZ_) + k;
		im_p    = (i_plus*NY*NZ_) + (j*NZ_) + k;

		uu[0][im_ob] =
		    (beta*work_v3[0][im] +
		     alpha*work_v3[0][im_p])
		    + degree_oblique*(beta*work_v3[1][im] +
				      alpha*work_v3[1][im_p])
		    + Shear_rate_eff*(j - NY/2.);

		uu[1][im_ob] =
		    beta*work_v3[1][im] +
		    alpha*work_v3[1][im_p];

		uu[2][im_ob] =
		    beta*work_v3[2][im] +
		    alpha*work_v3[2][im_p];
	    }
	}
    }
    
}

inline void contra2co(double **contra) {

    int im;

    for (int i = 0; i < NX; i++) {
	for (int j = 0; j < NY; j++) {
	    for (int k = 0; k < NZ_; k++) {
		im      = (i*NY*NZ_) + (j*NZ_) + k;
		
		work_v3[0][im] = contra[0][im];
		work_v3[1][im] = contra[1][im];
		work_v3[2][im] = contra[2][im];
	    }
	}
    }

    for (int i = 0; i < NX; i++) {
	for (int j = 0; j < NY; j++) {
	    for (int k = 0; k < NZ_; k++) {
		im      = (i*NY*NZ_) + (j*NZ_) + k;

		contra[0][im] =
		    work_v3[0][im] +
		    degree_oblique*work_v3[1][im];
		contra[1][im] =
		    degree_oblique*work_v3[0][im] +
		    (1. + degree_oblique*degree_oblique)*work_v3[1][im];
		contra[2][im] =
		    work_v3[2][im];
	    }
	}
    }
}

inline void co2contra(double **contra) {

    int im;

    for (int i = 0; i < NX; i++) {
	for (int j = 0; j < NY; j++) {
	    for (int k = 0; k < NZ_; k++) {
		im      = (i*NY*NZ_) + (j*NZ_) + k;

		work_v3[0][im] = contra[0][im];
		work_v3[1][im] = contra[1][im];
		work_v3[2][im] = contra[2][im];
	    }
	}
    }

    for (int i = 0; i < NX; i++) {
	for (int j = 0; j < NY; j++) {
	    for (int k = 0; k < NZ_; k++) {
		im      = (i*NY*NZ_) + (j*NZ_) + k;

		contra[0][im] =
		    (1. + degree_oblique*degree_oblique)*work_v3[0][im] -
		    degree_oblique*work_v3[1][im];
		contra[1][im] =
		    -degree_oblique*work_v3[0][im] +
		    work_v3[1][im];
		contra[2][im] =
		    work_v3[2][im];
	    }
	}
    }
}

inline void contra2co_single(double contra[]) {
    double dmy[DIM];

    for (int d = 0; d < DIM; d++) {
	dmy[d] = contra[d];
    }

    contra[0] =
	dmy[0] + degree_oblique*dmy[1];
    contra[1] =
	degree_oblique*dmy[0] + (1. + degree_oblique*degree_oblique)*dmy[1];
    contra[2] = dmy[2];
}

inline void co2contra_single(double co[]) {
    double dmy[DIM];

    for (int d = 0; d < DIM; d++) {
	dmy[d] = co[d];
    }

    co[0] =
	(1. + degree_oblique*degree_oblique)*dmy[0] - degree_oblique*dmy[1];
    co[1] =
	-degree_oblique*dmy[0] + dmy[1];
    co[2] = dmy[2];
}

inline void A2a_k(double *a){
  int im;
  
#ifndef _OPENMP
  double ***a_cp;
  a_cp = alloc_3d_double(NX, NY, NZ_);
#pragma omp parallel for schedule(dynamic,1) private(im)
  for (int i = 0; i< NX; i++){
  for (int j = 0; j< NY; j++){
  for (int l = 0; l< NZ; l++){
  im = (i*NY*NZ_)+(j*NZ_)+l; 
  a_cp[i][j][l]=a[im];
  }
  }
  }
  
  rdft3d(NX, NY, NZ, 1, a_cp, t, ip, w);
  rdft3dsort(NX, NY, NZ, 1, a_cp);
  
#pragma omp parallel for schedule(dynamic, 1) private(im) 
  for (int i = 0; i< NX; i++){
  for (int j = 0; j< NY; j++){
  for (int l = 0; l< NZ/2+1; l++){
  im = (i*NY*NZ_)+(j*NZ_)+2*l; 
  a[im]=a_cp[i][j][2*l];
  a[im+1]=a_cp[i][j][2*l+1];
  }
  }
  }
  free_3d_double(a_cp);
#endif

#ifdef _OPENMP
double x_in[NX][NY][NZ];
double _Complex x_out[NX][NY][NZ/2+1];

{
#pragma omp parallel for schedule(dynamic, 1) private(im)
for (int i = 0; i< NX; i++){
for (int j = 0; j< NY; j++){
for (int l = 0; l< NZ; l++){
im=(i*NY*NZ_)+(j*NZ_)+l;
x_in[i][j][l]=a[im];
}
}
}

{
DFTI_DESCRIPTOR_HANDLE Desc_Handle = 0;
long    m;
long    n;
long    k;
long    Status;
double  Scale;
long    lengths[3];
long    strides_in[4]; 
long    strides_out[4]; 

lengths[0] = (NX);
lengths[1] =(NY);
lengths[2] = (NZ);

strides_in[0] = 0;
strides_in[1] = (NZ)*NY;
strides_in[2] = NZ;
strides_in[3] = 1;

strides_out[0] = 0;
strides_out[1] = (NZ/2+1)*NY;
strides_out[2] = NZ/2+1;
strides_out[3] = 1;
Status = DftiCreateDescriptor(&Desc_Handle, DFTI_DOUBLE,
                                    DFTI_REAL, 3, lengths);
Status = DftiSetValue(Desc_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
Status = DftiSetValue(Desc_Handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
Status = DftiSetValue(Desc_Handle, DFTI_INPUT_STRIDES, strides_in);
Status = DftiSetValue(Desc_Handle, DFTI_OUTPUT_STRIDES, strides_out);
Status = DftiCommitDescriptor(Desc_Handle);
Status = DftiComputeForward(Desc_Handle, x_in, x_out);
Status = DftiFreeDescriptor( &Desc_Handle);
}


#pragma omp parallel for schedule(dynamic, 1) private(im)
for (int i = 0; i< NX; i++){
for (int j = 0; j< NY; j++){
for (int l = 0; l< NZ/2+1; l++){
im=(i*NY*NZ_)+(j*NZ_)+2*l;
a[im]=__real__(x_out[i][j][l]);
a[im+1]=-(__imag__(x_out[i][j][l]));
}
}
}
}
#endif 
}

inline void A_k2a(double *a){

int im;

#ifndef _OPENMP
 double ***a_cp;
 a_cp = alloc_3d_double(NX, NY, NZ_);
#pragma omp parallel for schedule(dynamic,1) private(im) 
 for (int i = 0; i< NX; i++){
 for (int j = 0; j< NY; j++){
 for (int l = 0; l< NZ/2+1; l++){
 im = (i*NY*NZ_)+(j*NZ_)+2*l;
 a_cp[i][j][2*l]=a[im];
 a_cp[i][j][2*l+1]=a[im+1];
 }
 }
 }

 static double scale = 2.0/(NX * NY * NZ);
 rdft3dsort(NX, NY, NZ, -1, a_cp);
 rdft3d(NX, NY, NZ, -1, a_cp, t, ip, w);

#pragma omp parallel for schedule(dynamic, 1) private(im) 
   for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	  for(int k=0; k<NZ/2+1; k++){
      im = (i*NY*NZ_)+(j*NZ_)+2*k;
	 a[im] = a_cp[i][j][2*k]*scale;
	 a[im+1] = a_cp[i][j][2*k+1]*scale;
	}
      }
   }
free_3d_double(a_cp);
#endif

#ifdef _OPENMP
static double scale = 1.0/(NX * NY * NZ);
double x_in[NX][NY][NZ];
double _Complex x_out[NX][NY][NZ/2+1];


{
#pragma omp parallel for schedule(dynamic,1) private(im)
for (int i = 0; i< NX; i++){
for (int j = 0; j< NY; j++){
for (int l = 0; l< NZ/2+1; l++){
im = (i*NY*NZ_)+(j*NZ_)+2*l;
__real__(x_out[i][j][l])=a[im];
__imag__(x_out[i][j][l])=-a[im+1];
}
}
}

{
DFTI_DESCRIPTOR_HANDLE Desc_Handle = 0;
long    Status;
long    lengths[3];
long    strides_in[4]; 
long    strides_out[4]; 

lengths[0] = (NX);
lengths[1] =(NY);
lengths[2] = (NZ);

strides_in[0] = 0;
strides_in[1] = NZ*NY;
strides_in[2] = NZ;
strides_in[3] = 1;

strides_out[0] = 0;
strides_out[1] = (NZ/2+1)*NY;
strides_out[2] = NZ/2+1;
strides_out[3] = 1;

Status = DftiCreateDescriptor( &Desc_Handle, DFTI_DOUBLE,
                                    DFTI_REAL, 3, lengths);
Status = DftiSetValue( Desc_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
Status = DftiSetValue(Desc_Handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
Status = DftiSetValue(Desc_Handle, DFTI_OUTPUT_STRIDES, strides_in);
Status = DftiSetValue(Desc_Handle, DFTI_INPUT_STRIDES, strides_out);
Status = DftiCommitDescriptor( Desc_Handle );
Status = DftiComputeBackward( Desc_Handle, x_out, x_in);
Status = DftiFreeDescriptor(&Desc_Handle);
}

#pragma omp parallel for schedule(dynamic, 1) private(im)
for (int i = 0; i< NX; i++){
for (int j = 0; j< NY; j++){
for (int l = 0; l< NZ; l++){
im = (i*NY*NZ_)+(j*NZ_)+l;
a[im]=x_in[i][j][l]*scale;
}
}
}

}
#endif
}

inline void A_k2a_out(double *a_k, double *a_x){ 
  int im;
#pragma omp parallel for schedule(dynamic,1) private(im)
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ_; k++){
		im=(i*NY*NZ_)+(j*NZ_)+k;
	    a_x[im] = a_k[im]; 
      }
    }
  }
  A_k2a(a_x);
}
inline void A2a_k_out(double *a_x, double *a_k){ 
  int im;
#pragma omp parallel for schedule(dynamic, 1) private(im) 
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ; k++){
	    im=(i*NY*NZ_)+(j*NZ_)+k;
	    a_k[im] = a_x[im]; 
      }
    }
  }
  A2a_k(a_k);
}
inline void A_k2da_k(double *a, double **da){
    A_k2dxa_k(a,da[0]);
    A_k2dya_k(a,da[1]);
    A_k2dza_k(a,da[2]);
}

inline int Calc_KX_Ooura(const int &i, const int &j, const int &k){
  return (i>HNX) ? i-NX:i;
}
inline int Calc_KY_Ooura(const int &i, const int &j, const int &k){
  assert(i < NX);
  assert(j < NY);
  return (j>HNY) ? j-NY:j;
}
inline int Calc_KZ_Ooura(const int &i, const int &j, const int &k){
  assert(i < NX);
  assert(j < NY);
  return k/2;
}
inline void Truncate_general(double *a, const Index_range &ijk_range){
  int im;
#pragma omp parallel for schedule(dynamic, 1) private(im) 
  for(int i=ijk_range.istart; i<=ijk_range.iend; i++){
    for(int j=ijk_range.jstart; j<=ijk_range.jend; j++){
      for(int k=ijk_range.kstart; k<=ijk_range.kend; k++){
	assert( (abs(Calc_KY_Ooura(i,j,k))>= TRN_Y || abs(Calc_KX_Ooura(i,j,k))>= TRN_X || Calc_KZ_Ooura(i,j,k)>= TRN_Z)); 
	im=(i*NY*NZ_)+(j*NZ_)+k;
	a[im] = 0.0;
      }
    }
  }
}
inline void Truncate_two_third_rule_ooura(double *a){
  static Index_range dmy_range;
  const int trn_z2=2*TRN_Z;
  int range[2];
  {
    dmy_range.istart=0;
    dmy_range.iend=NX-1;
    dmy_range.jstart=0;
    dmy_range.jend=NY-1;
    dmy_range.kstart=trn_z2;
    dmy_range.kend=NZ_-1;
    Truncate_general(a, dmy_range);
  }

  {
    dmy_range.istart=0;
    dmy_range.iend=NX-1;
    dmy_range.jstart=TRN_Y;
    dmy_range.jend=NY-TRN_Y;
    dmy_range.kstart=0;
    dmy_range.kend=trn_z2-1;
    Truncate_general(a, dmy_range);
  }

  {
    dmy_range.istart=TRN_X;
    dmy_range.iend=NX-TRN_X;
    dmy_range.jstart=0;
    dmy_range.jend=TRN_Y-1;
    dmy_range.kstart=0;
    dmy_range.kend=trn_z2-1;
    Truncate_general(a, dmy_range);
  }

  {
    dmy_range.istart=TRN_X;
    dmy_range.iend=NX-TRN_X;
    dmy_range.jstart=NY-TRN_Y+1;
    dmy_range.jend=NY-1;
    dmy_range.kstart=0;
    dmy_range.kend=trn_z2-1;
    Truncate_general(a, dmy_range);
  }
}

#endif
