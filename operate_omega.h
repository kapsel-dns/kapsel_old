//
// $Id: operate_omega.h,v 1.11 2005/07/27 09:48:06 nakayama Exp $
//
#ifndef OPERATE_OMEGA_H
#define OPERATE_OMEGA_H

#include "variable.h"
#include "fft_wrapper.h"

void Add_advection_pressurek(Value &pk
			     ,const Value u[DIM]
			     ,Value *uu // working memory [DIM*2]
			     );
void U_k2omega_k(const Value u[DIM], Value omega[DIM], double uk_dc[DIM]);
void Omega_k2u_k(const Value omega[DIM], const double uk_dc[DIM], Value u[DIM]);
void U2zeta_k(Value zeta[DIM-1], double uk_dc[DIM], Value u[DIM]);
void U2advection_k(Value u[DIM], Value advection[DIM-1]);
void Zeta_k2advection_k(Value zeta[DIM-1], const double uk_dc[DIM], Value advection[DIM-1]);
void Add_zeta_viscous_term(const Value zeta[DIM-1], Value f[DIM-1], const Index_range &ijk_range);
void Solenoidal_uk(Value *u);


inline void U2u_k(Value u[DIM]){
  for(int d=0;d<DIM;d++){
    A2a_k(u[d]);
  }
}
inline void U_k2u(Value u[DIM]){
  for(int d=0;d<DIM;d++){
    A_k2a(u[d]);
  }
}
inline void Solenoidal_u(Value *u){
  U2u_k(u);
  Solenoidal_uk(u);
  U_k2u(u);
}


inline void Truncate_vector_two_third_rule(Value *vector,const int &dim){
  for(int d=0; d<dim; d++){
      Truncate_two_third_rule(vector[d]);
  }
}
inline void Zeta_k2omega(const Value zeta[DIM-1], Value omega[DIM]){
  Zeta_k2omega_k(zeta, omega);
  U_k2u(omega);
}
inline void Zeta_k2u(const Value zeta[DIM-1], const double uk_dc[DIM], Value u[DIM]){
  Zeta_k2u_k(zeta, uk_dc, u);
  U_k2u(u);
}
inline void U2zeta_k(Value zeta[DIM-1], double uk_dc[DIM], Value u[DIM]){
  U2u_k(u);
  U_k2zeta_k(u, zeta, uk_dc);
}
inline void Symmetrize_fourier(Value &f){
  for(int i=0;i<NX;i++){
    for(int j=1;j<NY_shear;j++){
      for(int k=0;k<NZ_;k++){
	double dmy = (f[i][j][k]+f[i][NY-j][k])*.5;
	f[i][j][k] = dmy;
	f[i][NY-j][k] = dmy;
      }
    }
  }
}
inline void Antisymmetrize_fourier(Value &f){
  for(int i=0;i<NX;i++){
    for(int k=0;k<NZ_;k++){
	f[i][0][k] = 0.0;
	f[i][NY_shear][k] = 0.0;
    }
    for(int j=1;j<NY_shear;j++){
      for(int k=0;k<NZ_;k++){
	double dmy = (f[i][j][k]-f[i][NY-j][k])*.5;
	f[i][j][k] = dmy;
	f[i][NY-j][k] = -dmy;
      }
    }
  }
}
inline void Symmetrize_uk(Value u[DIM]){
    Symmetrize_fourier(u[0]);
    Antisymmetrize_fourier(u[1]);
    Symmetrize_fourier(u[2]);
}
inline void Symmetrize_zetak(Value zeta[DIM-1]){
    for(int i=0;i<NX;i++){
	for(int k=0;k<NZ_;k++){
	    if(KX_int[i][0][k] != 0){
		//omega[2], sin
		zeta[1][i][0][k] = 0.0;
	    }else if(KY_int[i][0][k] != 0){
		//omega[2], sin
		zeta[0][i][0][k] = 0.0;
		//omega[0], sin
		zeta[1][i][0][k] = 0.0;
	    }else{
		//omega[0], sin
		zeta[0][i][0][k] = 0.0;
	    }

	    if(KX_int[i][NY_shear][k] != 0){
		//omega[2], sin
		zeta[1][i][NY_shear][k] = 0.0;
	    }else if(KY_int[i][NY_shear][k] != 0){
		//omega[2], sin
		zeta[0][i][NY_shear][k] = 0.0;
		//omega[0], sin
		zeta[1][i][NY_shear][k] = 0.0;
	    }else{
		//omega[0], sin
		zeta[0][i][NY_shear][k] = 0.0;
	    }
	}
    }
    for(int i=0; i<NX; i++){
	for(int j=1; j<NY_shear; j++){
	    for(int k=0; k<NZ_; k++){
		double dmy;
		if(KX_int[i][j][k] != 0){
		    //omega[1], cos
		    dmy = (zeta[0][i][j][k]+zeta[0][i][NY-j][k])*.5;
		    zeta[0][i][j][k] = dmy;
		    zeta[0][i][NY-j][k] = dmy;
		    //omega[2], sin
		    dmy = (zeta[1][i][j][k]-zeta[1][i][NY-j][k])*.5;
		    zeta[1][i][j][k] = dmy;
		    zeta[1][i][NY-j][k] = -dmy;
		}else if(KY_int[i][j][k] != 0){
		    //omega[2], sin
		    dmy = (zeta[0][i][j][k]-zeta[0][i][NY-j][k])*.5;
		    zeta[0][i][j][k] = dmy;
		    zeta[0][i][NY-j][k] = -dmy;
		    //omega[0], sin
		    dmy = (zeta[1][i][j][k]-zeta[1][i][NY-j][k])*.5;
		    zeta[1][i][j][k] = dmy;
		    zeta[1][i][NY-j][k] = -dmy;
		}else{
		    //omega[0], sin
		    dmy = (zeta[0][i][j][k]-zeta[0][i][NY-j][k])*.5;
		    zeta[0][i][j][k] = dmy;
		    zeta[0][i][NY-j][k] = -dmy;
		    //omega[1], cos
		    dmy = (zeta[1][i][j][k]+zeta[1][i][NY-j][k])*.5;
		    zeta[1][i][j][k] = dmy;
		    zeta[1][i][NY-j][k] = dmy;
		}
	    }
	}
    }
    assert(zeta[0][0][0][0] == 0.0); 
    assert(zeta[1][0][0][0] == 0.0); 
}
inline void Add_pressurek(Value &pk, const Value &srck){
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ_; k++){
	pk[i][j][k] += (-IK2[i][j][k] * srck[i][j][k]);
      }
    }
  }
}
#endif
