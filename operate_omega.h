//
// $Id: operate_omega.h,v 1.1 2006/06/27 18:41:28 nakayama Exp $
//
#ifndef OPERATE_OMEGA_H
#define OPERATE_OMEGA_H

#include "variable.h"
#include "fft_wrapper.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void U2advection_k(double **u, double **advection);
void Zeta_k2advection_k(double **zeta, double uk_dc[DIM], double **advection);
void Zeta_k2advection_k_OBL(double **zeta, double uk_dc[DIM], double **advection);
void Add_zeta_viscous_term(double **zeta, double **f, const Index_range &ijk_range);
void Solenoidal_uk(double **u);
void Solenoidal_uk_OBL(double **u);

inline void U2u_k(double **u){
    A2a_k(u[0]);
    A2a_k(u[1]);
    A2a_k(u[2]);
}

inline void U_k2u(double **u){
    A_k2a(u[0]);
    A_k2a(u[1]);
    A_k2a(u[2]);
}

inline void Solenoidal_u(double **u){
  U2u_k(u);
  for(int d=0; d<DIM; d++){
      Truncate_two_third_rule(u[d]);
  }
  Solenoidal_uk(u);
  U_k2u(u);
}

inline void Truncate_vector_two_third_rule(double **vector,const int &dim){
  for(int d=0; d<dim; d++){
      Truncate_two_third_rule(vector[d]);
  }
}

inline void Zeta_k2u(double **zeta, double uk_dc[DIM], double **u){
  Zeta_k2u_k(zeta, uk_dc, u);
  U_k2u(u);
}

inline void Zeta_k2u_cpuky(double **zeta, double uk_dc[DIM], double **u, double *uk_cp){
  Zeta_k2u_k_OBL(zeta, uk_dc, u);

  int im;
  for(int i = 0; i < NX; i++){
      for(int j = 0; j < NY; j++){
	  for(int k = 0; k < NZ_; k++){
	      im=(i*NY*NZ_)+(j*NZ_)+k;
	      uk_cp[im] = u[1][im];
	  }
      }
  }
  U_k2u(u);
}

inline void Zeta_k2omega_OBL(double **zeta, double **omega){
    Zeta_k2omega_k_OBL(zeta, omega);
    
    for(int d = 0; d < DIM;d ++){
	A_k2a(omega[d]);
    }
}

inline void U2zeta_k(double **zeta, double uk_dc[DIM], double **u){
  U2u_k(u);
  U_k2zeta_k(u, zeta, uk_dc);
}
#endif
