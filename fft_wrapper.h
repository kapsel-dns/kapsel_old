//
// $Id: fft_wrapper.h,v 1.12 2005/08/12 15:25:46 nakayama Exp $
//
#ifndef FFT_WRAPPER_H
#define FFT_WRAPPER_H

#include <assert.h> 
#include "variable.h"
#include "input.h"
#include "alloc.h"
#include "macro.h"

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

extern int *ip_extended;
extern double *w_extended;
extern double *t_extended;

extern Value u[DIM], up[DIM], phi;
extern Value_int KX_int, KY_int, KZ_int;
extern Value K2,IK2;

extern int (*Calc_KX)( const int &i, const int &j, const int &k);
extern int (*Calc_KY)( const int &i, const int &j, const int &k);
extern int (*Calc_KZ)( const int &i, const int &j, const int &k);
extern void (*Truncate_two_third_rule)(Value &a);
extern void (*Truncate_four_fifths_rule)(Value &a);

//// function prototype
void Init_fft(void);
void A_k2dxa_k(const Value &a, Value &da);
void A_k2dya_k(const Value &a, Value &da);
void A_k2dza_k(const Value &a, Value &da);
void A_k2dxa_k_extended(const Value &a, Value &da);
void A_k2dya_k_extended(const Value &a, Value &da);
void A_k2dza_k_extended(const Value &a, Value &da);
void Omega_k2zeta_k(const Value omega[DIM], Value zeta[DIM-1]);
void U_k2zeta_k(const Value u[DIM], Value zeta[DIM-1], double uk_dc[DIM]);
void Zeta_k2u_k(const Value zeta[DIM-1], const double uk_dc[DIM], Value u[DIM]);
void Zeta_k2omega_k(const Value zeta[DIM-1], Value omega[DIM]);
void U_k2omega_k_inplace(Value u[DIM]);
void Zeta_k2Strain_k(const Value zeta[DIM-1], Value strain_k[QDIM]);
void U_k2divergence_k(const Value u[DIM], Value &div);


inline void A2a_k(Value &a){
  rdft3d(NX, NY, NZ, 1, a, t, ip, w);
  rdft3dsort(NX, NY, NZ, 1, a);
}
inline void A_k2a(Value &a){
  static double scale = 2.0/(NX * NY * NZ);
  rdft3dsort(NX, NY, NZ, -1, a);
  rdft3d(NX, NY, NZ, -1, a, t, ip, w);
  {
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ; k++){
	  a[i][j][k] *= scale;
	}
      }
    }
  }
}
inline void A_k2a_out(const Value &a_k, Value &a_x){ 
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ_; k++){
	a_x[i][j][k] = a_k[i][j][k]; 
      }
    }
  }
  A_k2a(a_x);
}
inline void A2a_k_out(const Value &a_x, Value &a_k){ 
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ; k++){
	  a_k[i][j][k] = a_x[i][j][k]; 
      }
    }
  }
  A2a_k(a_k);
}
inline void A_k2da_k(const Value &a, Value da[DIM]){
    A_k2dxa_k(a,da[0]);
    A_k2dya_k(a,da[1]);
    A_k2dza_k(a,da[2]);
}
inline void A_k2da_k_extended(const Value &a, Value da[DIM]){
    A_k2dxa_k_extended(a,da[0]);
    A_k2dya_k_extended(a,da[1]);
    A_k2dza_k_extended(a,da[2]);
}

inline void A2a_k_extended(Value &a){
  rdft3d(N2X, N2Y, N2Z, 1, a, t_extended, ip_extended, w_extended);
  rdft3dsort(N2X, N2Y, N2Z, 1, a);
}
inline void A_k2a_extended(Value &a){
  static double scale = 2./(N2X*N2Y*N2Z);
  rdft3dsort(N2X, N2Y, N2Z, -1, a);
  rdft3d(N2X, N2Y, N2Z, -1, a, t_extended, ip_extended, w_extended);
  {
    for(int i=0; i<N2X; i++){
      for(int j=0; j<N2Y; j++){
	for(int k=0; k<N2Z; k++){
	  a[i][j][k] *= scale;
	}
      }
    }
  }
}


//// function definition
inline int Calc_KX_extension_Ooura( const int &i, const int &j, const int &k){
  // i,j,k: 配列の添字, 
  assert(i < N2X);
  assert(j < N2Y);

  return (i>NX) ? i-N2X:i;
}
inline int Calc_KY_extension_Ooura( const int &i, const int &j, const int &k){
  // i,j,k: 配列の添字, 
  assert(i < N2X);
  assert(j < N2Y);
  
  return (j>NY) ? j-N2Y:j;
}
inline int Calc_KZ_extension_Ooura( const int &i, const int &j, const int &k){
  // i,j,k: 配列の添字, 
  assert(i < N2X);
  assert(j < N2Y);
  
  return k/2;
}


inline void K_space_extension(Value &a){
    { /// k < NZ_
	for(int i=HNX+1; i<NX; i++){ // right
	    for(int j=0; j<HNY; j++){ // bottom
		for(int k=0;k<NZ; k++){
		    int i2=NX+i;
		    a[i2][j][k] = a[i][j][k];
		}
	    }
	}
	for(int i=0; i<HNX; i++){ // left
	    for(int j=HNY+1; j<NY; j++){ // top
		for(int k=0;k<NZ; k++){
		    int j2=j + NY;
		    a[i][j2][k] = a[i][j][k];
		}
	    }
	}
	for(int i=HNX+1; i<NX; i++){ // right
	    for(int j=HNY+1; j<NY; j++){ // top
		for(int k=0;k<NZ; k++){
		    int i2=NX+i;
		    int j2=NY+j;
		    a[i2][j2][k] = a[i][j][k];
		}
	    }
	}
    }
    {// zero clear
	// zero clear  for kz >= NZ
	for(int i=0; i<N2X; i++){
	    for(int j=0; j<N2Y; j++){
		for(int k=NZ; k<N2Z_; k++){
		    assert(Calc_KZ_extension_Ooura(i,j,k) >= HNZ);
		    a[i][j][k] = 0.0;
		}
	    }
	}
	{// zero clear  for kz < NZ
	    // zero clear  
	    for(int i=0;i<HNX; i++){
		for(int j=HNY; j<=N2Y-HNY; j++){
		    for(int k=0;k<NZ; k++){
			assert((ABS(Calc_KY_extension_Ooura(i,j,k)) >=HNY));
			a[i][j][k] = 0.0;
		    }
		}
	    }
	    for(int i=HNX; i<=N2X-HNX; i++){
		for(int j=0; j<N2Y; j++){
		    for(int k=0;k<NZ; k++){
			assert((ABS(Calc_KX_extension_Ooura(i,j,k))>= HNX));
			a[i][j][k] = 0.0;
		    }
		}
	    }
	    for(int i=N2X-HNX; i<N2X; i++){
		for(int j=HNY; j<=N2Y-HNY; j++){
		    for(int k=0;k<NZ; k++){
			assert((ABS(Calc_KY_extension_Ooura(i,j,k)) >=HNY));
			a[i][j][k] = 0.0;
		    }
		}
	    }
	}
    }
}
inline void K_space_extension_reverse(Value &a){
    for(int i=HNX+1; i<NX; i++){
	for(int j=0; j<HNY; j++){
	    for(int k=0;k<NZ; k++){
		int kx = Calc_KX(i,j,k);
		int i2=N2X + kx;
		a[i][j][k] = a[i2][j][k];
	    }
	}
    }
    for(int i=0; i<HNX; i++){
	for(int j=HNY+1; j<NY; j++){
	    for(int k=0;k<NZ; k++){
		int ky = Calc_KY(i,j,k);
		int j2=N2Y + ky;
		a[i][j][k] = a[i][j2][k];
	    }
	}
    }
    for(int i=HNX+1; i<NX; i++){
	for(int j=HNY+1; j<NY; j++){
	    for(int k=0;k<NZ; k++){
		int kx = Calc_KX(i,j,k);
		int i2=N2X + kx;
		int ky = Calc_KY(i,j,k);
		int j2=N2Y + ky;
		a[i][j][k] = a[i2][j2][k];
	    }
	}
    }
    int i = HNX;
    for(int j=0; j<NY; j++){
	for(int k=0;k<NZ; k++){
		a[i][j][k] = 0.0;
	}
    }
    for(int i=0; i<NX; i++){
	int j = HNY;
	for(int k=0;k<NZ; k++){
		a[i][j][k] = 0.0;
	}
    }
    for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	    for(int k=NZ;k<NZ_; k++){
		a[i][j][k] = 0.0;
	    }
	}
    }
}
#endif
