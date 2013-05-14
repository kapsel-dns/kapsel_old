//
// $Id: operate_omega.cxx,v 1.1 2006/06/27 18:41:28 nakayama Exp $
//

#include "operate_omega.h"

const int DIM2 = 2*DIM;

void U2advection_k(double **u, double **advection){
    double u1;
    double u2;
    double u3;
    double k1 ;
    double k2 ;
    double k3 ;
    double u2u3;
    double u3u1;
    double u1u2;
    double u22_u32;
    double u32_u12;
    double u12_u22;
    double k1k2;
    double k2k3;
    double k3k1;
    int im;
    
    {
#pragma omp parallel for schedule(dynamic, 1) private(u1, u2, u3, im)
	for(int i=0; i<NX; i++){
	    for(int j=0; j<NY; j++){
		for(int k=0; k<NZ; k++){
		    im=(i*NY*NZ_)+(j*NZ_)+k;
		    u1=u[0][im];
		    u2=u[1][im];
		    u3=u[2][im];
		    
		    u[0][im] = u2 * u3;
		    u[1][im] = u3 * u1;
		    u[2][im] = u1 * u2;
		    advection[0][im] = (u2-u3) * (u2+u3);
		    advection[1][im] = (u3-u1) * (u3+u1);
		}
	    }
	}
	
	
	{
	    for(int d=0;d<DIM;d++){
		A2a_k(u[d]);
	    }
	    A2a_k(advection[0]);
	    A2a_k(advection[1]);
	}
	
#pragma omp parallel for schedule(dynamic, 1) private(k1,k2,k3,u2u3,u3u1,u1u2,u22_u32,u32_u12,u12_u22,k1k2,k2k3,k3k1,im)
	for(int i=0; i<NX; i++){
	    for(int j=0; j<NY; j++){
		for(int k=0; k<NZ_; k++){
		    im=(i*NY*NZ_)+(j*NZ_)+k;
		    k1 =KX_int[im] * WAVE_X;
		    k2 =KY_int[im] * WAVE_Y;
		    k3 =KZ_int[im] * WAVE_Z;
		    
		    u2u3 = u[0][im];
		    u3u1 = u[1][im];
		    u1u2 = u[2][im];
		    u22_u32 =advection[0][im];
		    u32_u12 =advection[1][im];
		    u12_u22 = -u32_u12 - u22_u32;
		    
		    k1k2 = k1*k2;
		    k2k3 = k2*k3;
		    k3k1 = k3*k1;
		    
		    u[0][im] = -(
			k2k3 * u22_u32
			+(k3-k2)*(k3+k2) *u2u3
			+k3k1 * u1u2 - k1k2 * u3u1
			);
		    u[1][im] = -(
			k3k1 * u32_u12
			+(k1-k3)*(k1+k3) *u3u1
			+ k1k2 * u2u3 - k2k3 * u1u2
			);
		    u[2][im] = -(
			k1k2 * u12_u22
			+(k2-k1)*(k2+k1) *u1u2
			+ k2k3 * u3u1 - k3k1 * u2u3
			);
		}
	    }
	}
	Omega_k2zeta_k(u, advection);
    }
}

void Zeta_k2advection_k(double **zeta, double uk_dc[DIM], double **advection){
    //Truncate_vector_two_third_rule(zeta, DIM-1);
    //Zeta_k2u(zeta, uk_dc, u);
    for(int d=0; d<(DIM-1); d++){
	Truncate_two_third_rule_ooura(zeta[d]);
    }
    Zeta_k2u_k(zeta, uk_dc, u);
    for(int d=0;d<DIM;d++){
	A_k2a(u[d]);
    }
    U2advection_k(u, advection);
}

void Zeta_k2advection_k_OBL(double **zeta, double uk_dc[DIM], double **advection){
    //Truncate_vector_two_third_rule(zeta, DIM-1);
    //Zeta_k2u(zeta, uk_dc, u);
    for(int d=0; d<(DIM-1); d++){
	Truncate_two_third_rule_ooura(zeta[d]);
    }
    Zeta_k2u_cpuky(zeta, uk_dc, u, work_v1);//contra
    Zeta_k2omega_OBL(zeta, work_v3);//contra

    int im;
    double u1;
    double u2;
    double u3;
    double w1;
    double w2;
    double w3;
#pragma omp parallel for schedule(dynamic, 1) private(u1,u2,u3,w1,w2,w3,im)
    for(int i = 0; i < NX; i++){
	for(int j = 0; j < NY; j++){
	    for(int k = 0; k < NZ; k++){
		im=(i*NY*NZ_)+(j*NZ_)+k;
		
		u1 = u[0][im];
		u2 = u[1][im];
		u3 = u[2][im];
		
		w1 = work_v3[0][im];
		w2 = work_v3[1][im];
		w3 = work_v3[2][im];
		
		u[0][im] = u2*w3 - u3*w2;//co
		u[1][im] = u3*w1 - u1*w3;//co
		u[2][im] = u1*w2 - u2*w1;//co
	    }
	}
    }
    
    for(int d=0;d<DIM;d++){
	A2a_k(u[d]);
    }
    
    for(int i = 0; i < NX; i++){
	for(int j = 0; j < NY; j++){
	    for(int k = 0; k < NZ_; k++) {
		im=(i*NY*NZ_)+(j*NZ_)+k;
		u[0][im] -= 2.*Shear_rate_eff*work_v1[im];//co
		u[1][im] -= 2.*Shear_rate_eff*degree_oblique*work_v1[im];//co
	    }
	}
    }
    
    
    U_k2rotation_k(u);//contra
    Omega_k2zeta_k_OBL(u, advection);//contra
    
}

void Add_zeta_viscous_term(double ** zeta, double **f, const Index_range &ijk_range){
  int im;
#pragma omp parallel for schedule(dynamic, 1) private(im)
    for(int i=ijk_range.istart; i<=ijk_range.iend; i++){
      for(int j=ijk_range.jstart; j<=ijk_range.jend; j++){
    	for(int k=ijk_range.kstart; k<=ijk_range.kend; k++){
	  im=(i*NY*NZ_)+(j*NZ_)+k;
	  f[0][im] += -(NU * K2[im]* zeta[0][im]);
	  f[1][im] += -(NU * K2[im]* zeta[1][im]);
	}
      }
    }
}

void Solenoidal_uk(double **u){
    
    double kx;
    double ky;
    double kz;
    double dmy;
    int im;
    
#pragma omp parallel for schedule(dynamic, 1) private(kx,ky,kz,dmy,im) 
    for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	    for(int k=0; k<NZ_; k++){
		
		im=(i*NY*NZ_)+(j*NZ_)+k;
		kx =KX_int[im] * WAVE_X;
		ky =KY_int[im] * WAVE_Y;
		kz =KZ_int[im] * WAVE_Z;
		
		dmy =
		    IK2[im] 
		    * (
			u[0][im] * kx
			+u[1][im] * ky
			+u[2][im] * kz
			);
		
		u[0][im] -= dmy * kx;
		u[1][im] -= dmy * ky;
		u[2][im] -= dmy * kz;
		
	    }
	}
    }
}


void Solenoidal_uk_OBL(double **u){

double kx;
double ky;
double kz;
double kx_contra;
double ky_contra;
double dmy;
int im;

#pragma omp parallel for schedule(dynamic, 1) private(kx,ky,kz,kx_contra,ky_contra,dmy,im) 
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ_; k++){

	im=(i*NY*NZ_)+(j*NZ_)+k;
	kx =KX_int[im] * WAVE_X;
	ky =KY_int[im] * WAVE_Y;
	kz =KZ_int[im] * WAVE_Z;

	kx_contra =
	    (1. + degree_oblique*degree_oblique)*KX_int[im]*WAVE_X -
	    degree_oblique*KY_int[im]*WAVE_Y;
	ky_contra =
	    -degree_oblique*KX_int[im]*WAVE_X +
	    KY_int[im]*WAVE_Y;
	
	dmy =
	  IK2[im] 
	  * (
	     u[0][im] * kx
	     +u[1][im] * ky
	     +u[2][im] * kz
	     );

	u[0][im] -= dmy * kx_contra;
	u[1][im] -= dmy * ky_contra;
	u[2][im] -= dmy * kz;
	
      }
    }
  }
}

