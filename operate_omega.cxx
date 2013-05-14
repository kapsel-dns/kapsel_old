//
// $Id: operate_omega.cxx,v 1.10 2005/07/27 09:48:44 nakayama Exp $
//

#include "operate_omega.h"

const int DIM2 = 2*DIM;
void Add_advection_pressurek(Value &pk
			     ,const Value u[DIM]
			     ,Value *uu // working memory [DIM*2]
			     ){
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ; k++){
	double ux = u[0][i][j][k];
	double uy = u[1][i][j][k];
	double uz = u[2][i][j][k];

	uu[0][i][j][k] = SQ(ux);
	uu[1][i][j][k] = SQ(uy);
	uu[2][i][j][k] = SQ(uz);

	uu[3][i][j][k] = ux*uy;
	uu[4][i][j][k] = uy*uz;
	uu[5][i][j][k] = uz*ux;
      }
    }
  }
  for(int d=0; d<DIM2; d++){
    A2a_k(uu[d]);
  }

  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ_; k++){

	double k1 =KX_int[i][j][k] * WAVE_X;
	double k2 =KY_int[i][j][k] * WAVE_Y;
	double k3 =KZ_int[i][j][k] * WAVE_Z;
	
	double u1u1 = uu[0][i][j][k];
	double u2u2 = uu[1][i][j][k];
	double u3u3 = uu[2][i][j][k];

	double u1u2 = uu[3][i][j][k];
	double u2u3 = uu[4][i][j][k];
	double u3u1 = uu[5][i][j][k];

	uu[0][i][j][k] = RHO *(SQ(k1)*u1u1 
			       + SQ(k2)*u2u2
			       + SQ(k3)*u3u3
			       + 2.* (k1*k2 * u1u2
				      + k2*k3 * u2u3
				      + k3*k1 * u3u1
				      )
			       );
      }
    }
  }
  Add_pressurek(pk, uu[0]);
}

void U_k2omega_k(const Value u[DIM], Value omega[DIM], double uk_dc[DIM]){
  for(int d=0; d<DIM; d++){
    uk_dc[d] = u[d][0][0][0]; 
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ_; k++){
	  omega[d][i][j][k] = u[d][i][j][k]; 
	}
      }
    }
  }
  U_k2omega_k_inplace(omega);
}
void Omega_k2u_k(const Value omega[DIM], const double uk_dc[DIM], Value u[DIM]){
  double admy[DIM];
  U_k2omega_k(omega, u, admy);
  for(int d=0; d<DIM; d++){
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<HNZ_; k++){
	  int k2 = 2 * k;
	  double ik2 = IK2[i][j][k2];
	  u[d][i][j][k2] *= ik2;
	  u[d][i][j][k2+1] *= ik2;
	}
      }
    }
    u[d][0][0][0] = uk_dc[d];  
  }
}

void U2advection_k(Value u[DIM], Value advection[DIM-1]){

  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ; k++){
	double u1=u[0][i][j][k];
	double u2=u[1][i][j][k];
	double u3=u[2][i][j][k];

	u[0][i][j][k] = u2 * u3;
	u[1][i][j][k] = u3 * u1;
	u[2][i][j][k] = u1 * u2;
	advection[0][i][j][k] = (u2-u3) * (u2+u3);
	advection[1][i][j][k] = (u3-u1) * (u3+u1);
      }
    }
  }
  U2u_k(u);
  A2a_k(advection[0]);
  A2a_k(advection[1]);

  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ_; k++){
	double k1 =KX_int[i][j][k] * WAVE_X;
	double k2 =KY_int[i][j][k] * WAVE_Y;
	double k3 =KZ_int[i][j][k] * WAVE_Z;

	double u2u3 = u[0][i][j][k];
	double u3u1 = u[1][i][j][k];
	double u1u2 = u[2][i][j][k];
	double u22_u32 =advection[0][i][j][k];
	double u32_u12 =advection[1][i][j][k];
	double u12_u22 = -u32_u12 - u22_u32;

	double k1k2 = k1*k2;
	double k2k3 = k2*k3;
	double k3k1 = k3*k1;

	u[0][i][j][k] = -(
			  k2k3 * u22_u32
			  +(k3-k2)*(k3+k2) *u2u3
			  +k3k1 * u1u2 - k1k2 * u3u1
			  );
	u[1][i][j][k] = -(
			  k3k1 * u32_u12
			  +(k1-k3)*(k1+k3) *u3u1
			  + k1k2 * u2u3 - k2k3 * u1u2
			  );
	u[2][i][j][k] = -(
			  k1k2 * u12_u22
			  +(k2-k1)*(k2+k1) *u1u2
			  + k2k3 * u3u1 - k3k1 * u2u3
			  );
      }
    }
  }
  Omega_k2zeta_k(u, advection);
}


void Zeta_k2advection_k(Value zeta[DIM-1], const double uk_dc[DIM], Value advection[DIM-1]){
  Truncate_vector_two_third_rule(zeta, DIM-1);
  Zeta_k2u(zeta, uk_dc, u);

  U2advection_k(u, advection);
}

void Add_zeta_viscous_term(const Value zeta[DIM-1], Value f[DIM-1], const Index_range &ijk_range){
  for(int d=0; d<DIM-1; d++){
    for(int i=ijk_range.istart; i<=ijk_range.iend; i++){
      for(int j=ijk_range.jstart; j<=ijk_range.jend; j++){
    	for(int k=ijk_range.kstart; k<=ijk_range.kend; k++){
	  f[d][i][j][k] += -(NU * K2[i][j][k]* zeta[d][i][j][k]);
	}
      }
    }
  }
}
void Solenoidal_uk(Value *u){
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ_; k++){
	double kx =KX_int[i][j][k] * WAVE_X;
	double ky =KY_int[i][j][k] * WAVE_Y;
	double kz =KZ_int[i][j][k] * WAVE_Z;
	
	double dmy =
	  IK2[i][j][k] 
	  * (
	     u[0][i][j][k] * kx
	     +u[1][i][j][k] * ky
	     +u[2][i][j][k] * kz
	     );
	u[0][i][j][k] -= dmy * kx;
	u[1][i][j][k] -= dmy * ky;
	u[2][i][j][k] -= dmy * kz;
	
      }
    }
  }
}
