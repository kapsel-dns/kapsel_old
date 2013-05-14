//
// $Id: fft_wrapper.cxx,v 1.1 2006/06/27 18:41:28 nakayama Exp $
//

#include "fft_wrapper.h"

int *ip;
double *w;
double *t;


double **ucp;
double *phi, **up, **u, *rhop;
double **work_v3, *work_v1;

int *KX_int, *KY_int, *KZ_int;
double *K2, *IK2;

int (*Calc_KX)( const int &i, const int &j, const int &k);
int (*Calc_KY)( const int &i, const int &j, const int &k);
int (*Calc_KZ)( const int &i, const int &j, const int &k);
void (*Truncate_two_third_rule)(double *a);


inline void Init_K(void){
  KX_int = alloc_1d_int(NX*NY*NZ_);
  KY_int = alloc_1d_int(NX*NY*NZ_);
  KZ_int = alloc_1d_int(NX*NY*NZ_);
  K2 = alloc_1d_double(NX*NY*NZ_);
  IK2 = alloc_1d_double(NX*NY*NZ_);

  int kx, ky, kz;
  int im;
#pragma omp parallel for schedule(dynamic, 1) private(kx,ky,kz,im) 
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ_; k++){
	im=(i*NY*NZ_)+(j*NZ_)+k;
	kx = Calc_KX(i,j,k);
	ky = Calc_KY(i,j,k);
	kz = Calc_KZ(i,j,k);
	
	KX_int[im] = kx;
	KY_int[im] = ky;
	KZ_int[im] = kz;
	
	K2[im] = SQ(WAVE_X*(double)kx) 
	  + SQ(WAVE_Y*(double)ky)
	  + SQ(WAVE_Z*(double)kz);
	if(K2[im] > 0.0){
	  IK2[im] = 1./K2[im];
	}else{
	  IK2[im] = 0.0;
	}
      }
    }
  }
}

inline void Init_fft_ooura(void){
#ifndef _OPENMP
  fprintf(stderr,"# Ooura rdft3d is selected.\n");
#endif
#ifdef _OPENMP
  fprintf(stderr,"# Intel Math Kernel Library is selected.\n");
  fprintf(stderr,"# %d Threads\n", omp_get_max_threads());
#endif
  
  {
    int n, nt;
    nt=MAX(NX,NY);
    n=MAX(nt,NZ/2);
    ip = alloc_1d_int(2 + (int) sqrt((double)n + 0.5));
    t = alloc_1d_double( 8 * nt );
    w = alloc_1d_double( n/2 + NZ/4 ) ; 
    ip[0] = 0;
  }
  
  Calc_KX = Calc_KX_Ooura;
  Calc_KY = Calc_KY_Ooura;
  Calc_KZ = Calc_KZ_Ooura;

  Truncate_two_third_rule = Truncate_two_third_rule_ooura;

  Init_K();
}

void Init_fft(void){
  if(SW_FFT==Ooura){ // Ooura fft
    Init_fft_ooura();
  }else if(SW_FFT==IMKL_FFT){
    Init_fft_ooura();
  }else{
    fprintf(stderr,"specify SW_FFT correctly.\n");
    exit_job(EXIT_FAILURE);
  }
}

inline void A_k2dja_k_primitive(double *a
			 ,double *da
			 ,const int &nx
			 ,const int &ny
			 ,const int &hnz_
			 ,const int *k_int
			 ,const double &wave_base
			 ){
  int k2;
  int im;
  double wavenumber;
#pragma omp parallel for schedule(dynamic, 1) private(k2, wavenumber,im)
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      for(int k=0; k<hnz_; k++){
	k2=2*k;
	im=(i*NY*NZ_)+(j*NZ_)+k2; 
	wavenumber = k_int[im] * wave_base;
	da[im] = wavenumber * a[im+1];
	da[im+1] = -wavenumber * a[im];
      }
    }
  }
}
void A_k2dxa_k(double *a, double *da){
  A_k2dja_k_primitive(a,da,NX,NY,HNZ_,KX_int, WAVE_X);
}
void A_k2dya_k(double *a, double *da){
  A_k2dja_k_primitive(a,da,NX,NY,HNZ_,KY_int, WAVE_Y);
}
void A_k2dza_k(double *a, double *da){
  A_k2dja_k_primitive(a,da,NX,NY,HNZ_,KZ_int, WAVE_Z);
}

void Omega_k2zeta_k(double **omega, double **zetak){
  int im;
{
#pragma omp parallel for schedule(dynamic, 1) private(im)
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ_; k++){
	im=(i*NY*NZ_)+(j*NZ_)+k;
	/*
	if(KX_int[im] != 0){
	  zetak[0][im] = omega[1][im];
	  zetak[1][im] = omega[2][im];
	}else if(KZ_int[im] != 0){
	  zetak[0][im] = omega[1][im];
	  zetak[1][im] = omega[0][im];
	}else{
	  zetak[0][im] = omega[0][im];
	  zetak[1][im] = omega[2][im];
	}
	*/
	if(KX_int[im] != 0){
	  zetak[0][im] = omega[1][im];
	  zetak[1][im] = omega[2][im];
	}else if(KY_int[im] != 0){
	  zetak[0][im] = omega[2][im];
	  zetak[1][im] = omega[0][im];
	}else{
	  zetak[0][im] = omega[0][im];
	  zetak[1][im] = omega[1][im];
	}
      }
    }
  }
}
  assert(zetak[0][0] == 0.0); 
  assert(zetak[1][0] == 0.0); 
}

void Omega_k2zeta_k_OBL(double **omega, double **zetak){
  int im;
{
#pragma omp parallel for schedule(dynamic, 1) private(im)
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ_; k++){
	im=(i*NY*NZ_)+(j*NZ_)+k;
	if(KX_int[im] != 0){
	  zetak[0][im] = omega[1][im];
	  zetak[1][im] = omega[2][im];
	}else if(KZ_int[im] != 0){
	  zetak[0][im] = omega[1][im];
	  zetak[1][im] = omega[0][im];
	}else{
	  zetak[0][im] = omega[0][im];
	  zetak[1][im] = omega[2][im];
	}
      }
    }
  }
}
  assert(zetak[0][0] == 0.0); 
  assert(zetak[1][0] == 0.0); 
}

void Zeta_k2Strain_k(double **zeta, double *strain_k[QDIM]){
  // Strain_k は 5成分
  double dmy[DIM]={0.,0.,0.};
  int k2;
  int im0;
  int im1;
  double ks[DIM];
  double u_dmy[DIM][2];

  Zeta_k2u_k(zeta, dmy, strain_k);

#pragma omp parallel for schedule(dynamic, 1) private(im0, im1, k2, ks, u_dmy)
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<HNZ_; k++){
	k2=2*k;
	im0=(i*NY*NZ_)+(j*NZ_)+k2;
	im1=(i*NY*NZ_)+(j*NZ_)+k2+1;
	ks[0] = KX_int[(i*NY*NZ_)+(j*NZ_)+k2] * WAVE_X * .5;
	ks[1] = KY_int[(i*NY*NZ_)+(j*NZ_)+k2] * WAVE_Y * .5;
	ks[2] = KZ_int[(i*NY*NZ_)+(j*NZ_)+k2] * WAVE_Z * .5;
	for(int d=0;d<DIM;d++){
	  u_dmy[d][0] = u[d][im0];
	  u_dmy[d][1] = u[d][im1];
	}

	strain_k[0][im0] = ks[0] * u_dmy[0][1] * 2.;
	strain_k[0][im1] = -ks[0] * u_dmy[0][0] * 2.;

	strain_k[1][im0] = 
	  ks[0] * u_dmy[1][1]
	  +ks[1] * u_dmy[0][1];
	strain_k[1][im1] = 
	  - ks[0] * u_dmy[1][0] 
	  - ks[1] * u_dmy[0][0];

	strain_k[2][im0] = 
	  ks[0] * u_dmy[2][1]
	  +ks[2] * u_dmy[0][1];
	strain_k[2][im1] = 
	  - ks[0] * u_dmy[2][0] 
	  - ks[2] * u_dmy[0][0];

	strain_k[3][im0] = ks[1] * u_dmy[1][1] * 2.;
	strain_k[3][im1] = -ks[1] * u_dmy[1][0] * 2.;

	strain_k[4][im0] = 
	  ks[1] * u_dmy[2][1]
	  +ks[2] * u_dmy[1][1];
	strain_k[4][im1] = 
	  -ks[1] * u_dmy[2][0] 
	  -ks[2] * u_dmy[1][0];
      }
    }
  }
}

void U_k2zeta_k(double **u, double **zeta, double uk_dc[DIM]){
    double ks[DIM];
    double omega_re[DIM];
    double omega_im[DIM];
    int k2;
    int im;
    {
	uk_dc[0] = u[0][0]; 
	uk_dc[1] = u[1][0]; 
	uk_dc[2] = u[2][0]; 
#pragma omp parallel for schedule(dynamic, 1) private(ks, omega_re, omega_im, k2, im)
	for(int i=0; i<NX; i++){
	    for(int j=0; j<NY; j++){
		for(int k=0; k<HNZ_; k++){
		    k2=2*k;
		    im=(i*NY*NZ_)+(j*NZ_)+k2;
		    ks[0] = KX_int[im] * WAVE_X;
		    ks[1] = KY_int[im] * WAVE_Y;
		    ks[2] = KZ_int[im] * WAVE_Z;
		    omega_re[0] = 
			(ks[1] * u[2][im+1] 
			 - ks[2] * u[1][im+1]);
		    omega_im[0] = 
			-(ks[1] * u[2][im] 
			  - ks[2] * u[1][im]
			    );
		    
		    omega_re[1] = 
			(ks[2] * u[0][im+1] 
			 - ks[0] * u[2][im+1]);
		    omega_im[1] = 
			-(ks[2] * u[0][im] 
			  - ks[0] * u[2][im]
			    );
		    
		    omega_re[2] = 
			(ks[0] * u[1][im+1] 
			 - ks[1] * u[0][im+1]);
		    omega_im[2] = 
			-(ks[0] * u[1][im] 
			  - ks[1] * u[0][im]
			    );
		    
		    if(KX_int[im] != 0){
			zeta[0][im] = omega_re[1];
			zeta[0][im+1] = omega_im[1];
			zeta[1][im] = omega_re[2];
			zeta[1][im+1] = omega_im[2];
		    }else if(KY_int[im] != 0){
			zeta[0][im] = omega_re[2];
			zeta[0][im+1] = omega_im[2];
			zeta[1][im] = omega_re[0];
			zeta[1][im+1] = omega_im[0];
		    }else{
			zeta[0][im] = omega_re[0];
			zeta[0][im+1] = omega_im[0];
			zeta[1][im] = omega_re[1];
			zeta[1][im+1] = omega_im[1];
		    }
		    
		    /*
		    if(KX_int[im] != 0){
			zeta[0][im] = omega_re[1];
			zeta[0][im + 1] = omega_im[1];
			zeta[1][im] = omega_re[2];
			zeta[1][im + 1] = omega_im[2];
		    }else if(KZ_int[im] != 0){
			zeta[0][im] = omega_re[1];
			zeta[0][im + 1] = omega_im[1];
			zeta[1][im] = omega_re[0];
			zeta[1][im + 1] = omega_im[0];
		    }else{
			zeta[0][im] = omega_re[0];
			zeta[0][im + 1] = omega_im[0];
			zeta[1][im] = omega_re[2];
			zeta[1][im + 1] = omega_im[2];
		    }
		    */
		}
	    }
	}
    }
    assert(zeta[0][0] == 0.0); 
    assert(zeta[1][0] == 0.0); 
}

void U_k2zeta_k_OBL(double **u, double **zeta, double uk_dc[DIM]){
    double ks[DIM];
    double omega_re[DIM];
    double omega_im[DIM];
    int k2;
    int im;
    {
	uk_dc[0] = u[0][0]; 
	uk_dc[1] = u[1][0]; 
	uk_dc[2] = u[2][0]; 
	co2contra_single(uk_dc);
#pragma omp parallel for schedule(dynamic, 1) private(ks, omega_re, omega_im, k2, im)
	for(int i=0; i<NX; i++){
	    for(int j=0; j<NY; j++){
		for(int k=0; k<HNZ_; k++){
		    k2=2*k;
		    im=(i*NY*NZ_)+(j*NZ_)+k2;
		    ks[0] = KX_int[im] * WAVE_X;
		    ks[1] = KY_int[im] * WAVE_Y;
		    ks[2] = KZ_int[im] * WAVE_Z;
		    omega_re[0] = 
			(ks[1] * u[2][im+1] 
			 - ks[2] * u[1][im+1]);
		    omega_im[0] = 
			-(ks[1] * u[2][im] 
			  - ks[2] * u[1][im]
			    );
		    
		    omega_re[1] = 
			(ks[2] * u[0][im+1] 
			 - ks[0] * u[2][im+1]);
		    omega_im[1] = 
			-(ks[2] * u[0][im] 
			  - ks[0] * u[2][im]
			    );
		    
		    omega_re[2] = 
			(ks[0] * u[1][im+1] 
			 - ks[1] * u[0][im+1]);
		    omega_im[2] = 
			-(ks[0] * u[1][im] 
			  - ks[1] * u[0][im]
			    );
		    if(KX_int[im] != 0){
			zeta[0][im] = omega_re[1];
			zeta[0][im + 1] = omega_im[1];
			zeta[1][im] = omega_re[2];
			zeta[1][im + 1] = omega_im[2];
		    }else if(KZ_int[im] != 0){
			zeta[0][im] = omega_re[1];
			zeta[0][im + 1] = omega_im[1];
			zeta[1][im] = omega_re[0];
			zeta[1][im + 1] = omega_im[0];
		    }else{
			zeta[0][im] = omega_re[0];
			zeta[0][im + 1] = omega_im[0];
			zeta[1][im] = omega_re[2];
			zeta[1][im + 1] = omega_im[2];
		    }
		}
	    }
	}
    }
    assert(zeta[0][0] == 0.0); 
    assert(zeta[1][0] == 0.0); 
}

void Zeta_k2u_k(double **zeta, double uk_dc[DIM], double **u){
  double omega_re[DIM];
  double omega_im[DIM];
  double dmy1_re, dmy2_re;
  double dmy1_im, dmy2_im;
  double kx;
  double ky;
  double kz;
  double ik2;
  int im;
  int k2;
  double dmy;

  {
#pragma omp parallel for schedule(dynamic, 1) private(omega_re,omega_im, dmy1_re,dmy2_re,dmy1_im,dmy2_im, kx, ky, kz,ik2,im,k2,dmy)
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<HNZ_; k++){
	k2=2*k;
	im=(i*NY*NZ_)+(j*NZ_)+k2;
	kx = WAVE_X * KX_int[im];
	ky = WAVE_Y * KY_int[im];
	kz = WAVE_Z * KZ_int[im];
	ik2 = IK2[im];

	dmy1_re =zeta[0][im];
	dmy1_im =zeta[0][im+1];
	dmy2_re =zeta[1][im];
	dmy2_im =zeta[1][im+1];
	
	if(KX_int[im] != 0){
	  omega_re[0] = -(1./kx)*(ky * dmy1_re + kz * dmy2_re);
	  omega_im[0] = -(1./kx)*(ky * dmy1_im + kz * dmy2_im);
	  omega_re[1] = dmy1_re;
	  omega_im[1] = dmy1_im;
	  omega_re[2] = dmy2_re;
	  omega_im[2] = dmy2_im;
	}else if(KY_int[im] != 0){
	  dmy = -(kz/ky);
	  omega_re[0] = dmy2_re;
	  omega_im[0] = dmy2_im;
	  omega_re[1] = dmy * dmy1_re;
	  omega_im[1] = dmy * dmy1_im;
	  omega_re[2] = dmy1_re;
	  omega_im[2] = dmy1_im;
	}else{
	  omega_re[0] = dmy1_re;
	  omega_im[0] = dmy1_im;
	  omega_re[1] = dmy2_re;
	  omega_im[1] = dmy2_im;
	  omega_re[2] = 0.0;
	  omega_im[2] = 0.0;
	}
	
	/*
	if(KX_int[im] != 0){
	    omega_re[0] = -(1./kx)*(ky * dmy1_re + kz * dmy2_re);
	    omega_im[0] = -(1./kx)*(ky * dmy1_im + kz * dmy2_im);
	    omega_re[1] = dmy1_re;
	    omega_im[1] = dmy1_im;
	    omega_re[2] = dmy2_re;
	    omega_im[2] = dmy2_im;
	    
	    
	}else if(KZ_int[im] != 0){
	    dmy = -(ky/kz);
	    omega_re[0] = dmy2_re;
	    omega_im[0] = dmy2_im;
	    omega_re[1] = dmy1_re;
	    omega_im[1] = dmy1_im;
	    omega_re[2] = dmy*dmy1_re;
	    omega_im[2] = dmy*dmy1_im;
	}else{
	    omega_re[0] = dmy1_re;
	    omega_im[0] = dmy1_im;
	    omega_re[1] = 0.0;
	    omega_im[1] = 0.0;
	    omega_re[2] = dmy2_re;
	    omega_im[2] = dmy2_im;
	}
	*/
	kx *= ik2;
	ky *= ik2;
	kz *= ik2;
	u[0][im] = 
	  (ky * omega_im[2] - kz * omega_im[1]); 
	u[0][im+1] = 
	  -(ky * omega_re[2] - kz * omega_re[1]); 
	u[1][im] = 
	  (kz * omega_im[0] - kx * omega_im[2]); 
	u[1][im+1] = 
	  -(kz * omega_re[0] - kx * omega_re[2]); 
	u[2][im] = 
	  (kx * omega_im[1] - ky * omega_im[0]); 
	u[2][im+1] = 
	  -(kx * omega_re[1] - ky * omega_re[0]); 
      }
    }
  }
  u[0][0] = uk_dc[0];
  u[1][0] = uk_dc[1];
  u[2][0] = uk_dc[2];
}
}

void Zeta_k2u_k_OBL(double **zeta, double uk_dc[DIM], double **u){
  double omega_re[DIM];
  double omega_im[DIM];
  double dmy1_re, dmy2_re;
  double dmy1_im, dmy2_im;
  double kx;
  double ky;
  double kz;
  double ik2;
  int im;
  int k2;
  double dmy;

  {
#pragma omp parallel for schedule(dynamic, 1) private(omega_re,omega_im, dmy1_re,dmy2_re,dmy1_im,dmy2_im, kx, ky, kz,ik2,im,k2,dmy)
      for(int i=0; i<NX; i++){
	  for(int j=0; j<NY; j++){
	      for(int k=0; k<HNZ_; k++){
		  k2=2*k;
		  im=(i*NY*NZ_)+(j*NZ_)+k2;
		  kx = WAVE_X * KX_int[im];
		  ky = WAVE_Y * KY_int[im];
		  kz = WAVE_Z * KZ_int[im];
		  ik2 = IK2[im];
		  
		  dmy1_re =zeta[0][im];
		  dmy1_im =zeta[0][im+1];
		  dmy2_re =zeta[1][im];
		  dmy2_im =zeta[1][im+1];
		  
		  if(KX_int[im] != 0){
		      omega_re[0] = -(1./kx)*(ky * dmy1_re + kz * dmy2_re);
		      omega_im[0] = -(1./kx)*(ky * dmy1_im + kz * dmy2_im);
		      omega_re[1] = dmy1_re;
		      omega_im[1] = dmy1_im;
		      omega_re[2] = dmy2_re;
		      omega_im[2] = dmy2_im;
		      
		      
		  }else if(KZ_int[im] != 0){
		      dmy = -(ky/kz);
		      omega_re[0] = dmy2_re;
		      omega_im[0] = dmy2_im;
		      omega_re[1] = dmy1_re;
		      omega_im[1] = dmy1_im;
		      omega_re[2] = dmy*dmy1_re;
		      omega_im[2] = dmy*dmy1_im;
		  }else{
		      omega_re[0] = dmy1_re;
		      omega_im[0] = dmy1_im;
		      omega_re[1] = 0.0;
		      omega_im[1] = 0.0;
		      omega_re[2] = dmy2_re;
		      omega_im[2] = dmy2_im;
		  }
		  
		  //From contravariant to covariant
		  contra2co_single(omega_re);
		  contra2co_single(omega_im);
		  
		  kx *= ik2;
		  ky *= ik2;
		  kz *= ik2;
		  u[0][im] = 
		      (ky * omega_im[2] - kz * omega_im[1]);//contra
		  u[0][im+1] = 
		      -(ky * omega_re[2] - kz * omega_re[1]);//contra
		  u[1][im] = 
		      (kz * omega_im[0] - kx * omega_im[2]);//contra 
		  u[1][im+1] = 
		      -(kz * omega_re[0] - kx * omega_re[2]);//contra
		  u[2][im] = 
		      (kx * omega_im[1] - ky * omega_im[0]);//contra
		  u[2][im+1] = 
		      -(kx * omega_re[1] - ky * omega_re[0]);//contra
	      }
	  }
      }
      u[0][0] = uk_dc[0];
      u[1][0] = uk_dc[1];
      u[2][0] = uk_dc[2];
  }
}

void Zeta_k2omega_k_OBL(double **zeta, double **omega){
    double dmy1, dmy2;
    double kx;
    double ky;
    double kz;
    int im;
    int k2;
    double dmy;
    
    {
#pragma omp parallel for schedule(dynamic, 1) private(dmy1, dmy2, kx, ky, kz, im, k2, dmy)
	for(int i=0; i<NX; i++){
	    for(int j=0; j<NY; j++){
		for(int k=0; k<NZ_; k++){
		    im=(i*NY*NZ_)+(j*NZ_)+k;
		    kx = WAVE_X * KX_int[im];
		    ky = WAVE_Y * KY_int[im];
		    kz = WAVE_Z * KZ_int[im];
		    
		    dmy1 = zeta[0][im];
		    dmy2 = zeta[1][im];
		    
		    if(KX_int[im] != 0){

			omega[0][im]   = -(1./kx)*(ky * dmy1 + kz * dmy2);
			omega[1][im]   = dmy1;
			omega[2][im]   = dmy2;

		    }else if(KZ_int[im] != 0){

			dmy = -(ky/kz);
			omega[0][im]   = dmy2;
			omega[1][im]   = dmy1;
			omega[2][im]   = dmy*dmy1;

		    }else{

			omega[0][im]   = dmy1;
			omega[1][im]   = 0.0;
			omega[2][im]   = dmy2;

		    }
		}
	    }
	}
	
	for(int d=0; d<DIM; d++){
	    omega[d][0] = 0.0; 
	}
    }
}

void U_k2divergence_k(double **u, double *div){
    double ks[DIM];
    int k2;
    int im0, im1;
  {
#pragma omp parallel for schedule(dynamic, 1) private(ks, k2, im0, im1)
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<HNZ_; k++){
	k2=2*k;
	ks[DIM];
	im0=(i*NY*NZ_)+(j*NZ_)+k2;
	im1=im0+1;

	ks[0] = KX_int[im0] * WAVE_X;
	ks[1] = KY_int[im0] * WAVE_Y;
	ks[2] = KZ_int[im0] * WAVE_Z;
	//div[im0] = 0.0;
	//div[im1] = 0.0;
	  div[im0] = ks[0] * u[0][im1]
	            +ks[1] * u[1][im1]
	            +ks[2] * u[2][im1];
	  div[im1] = -ks[0] * u[0][im0]
	             -ks[1] * u[1][im0]
	             -ks[2] * u[2][im0];
      }
    }
  }
  }
}

void U_k2rotation_k(double **u){
    double ks[DIM];
    double dmy_u_re[DIM];
    double dmy_u_im[DIM];
    double dmy_rot_re[DIM];
    double dmy_rot_im[DIM];
    
    int k2;
    int im;
    {
#pragma omp parallel for schedule(dynamic, 1) private(ks, dmy_rot_re, dmy_rot_im, k2, im)
	for(int i = 0; i < NX; i++){
	    for(int j = 0; j < NY; j++){
		for(int k = 0; k < HNZ_; k++){
		    k2 = 2*k;
		    im = (i*NY*NZ_) + (j*NZ_) + k2;
		    
		    ks[0] = KX_int[im] * WAVE_X;
		    ks[1] = KY_int[im] * WAVE_Y;
		    ks[2] = KZ_int[im] * WAVE_Z;
		    
		    for(int d=0; d<DIM; d++){
			dmy_u_re[d] =u[d][im];
			dmy_u_im[d] =u[d][im + 1];
		    }
		    
		    dmy_rot_re[0] = (ks[1] * dmy_u_im[2] - ks[2] * dmy_u_im[1]);
		    dmy_rot_im[0] = -(ks[1] * dmy_u_re[2] - ks[2] * dmy_u_re[1]);
		    
		    dmy_rot_re[1] = (ks[2] * dmy_u_im[0] - ks[0] * dmy_u_im[2]);
		    dmy_rot_im[1] = -(ks[2] * dmy_u_re[0] - ks[0] * dmy_u_re[2]);
		    
		    dmy_rot_re[2] = (ks[0] * dmy_u_im[1] - ks[1] * dmy_u_im[0]);
		    dmy_rot_im[2] = -(ks[0] * dmy_u_re[1] - ks[1] * dmy_u_re[0]);
		    
		    for(int d=0; d<DIM; d++){
			u[d][im] = dmy_rot_re[d];
			u[d][im + 1] = dmy_rot_im[d];
		    }
		}
	    }
	}
    }
}
