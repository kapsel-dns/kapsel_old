//
// $Id: f_particle.h,v 1.2 2006/08/15 15:01:35 nakayama Exp $
//
#ifndef F_PARTICLE_H
#define F_PARTICLE_H

#include "variable.h"
#include "operate_omega.h"
#include "input.h"
#include "make_phi.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern double **f_particle;

void Mem_alloc_f_particle(void);

inline void Make_f_particle_dt_nonsole(double **f
				       ,double **u 
				       ,double **up
				       ,double *phi
				       ){
  // !! これを呼ぶ前に、 up, phi を計算しておくこと。
  //
  // f = up - phi *u
  //
  int im;
  {
#pragma omp parallel for schedule(dynamic, 1) private(im)
    for(int i=0; i<NX; i++){
     for(int j=0; j<NY; j++){
	  for(int k=0; k<NZ; k++){
	  im=(i*NY*NZ_)+(j*NZ_)+k;
	  f[0][im] = up[0][im] - phi[im] * u[0][im];
	  f[1][im] = up[1][im] - phi[im] * u[1][im];
	  f[2][im] = up[2][im] - phi[im] * u[2][im];
	  }
     }
    }	
  }
}
inline void Make_f_particle_dt_sole(double **f
				    ,double **u
				    ,double **up
				    ,double *phi
				    ){
    // !! これを呼ぶ前に、 up, phi を計算しておくこと。
    Make_f_particle_dt_nonsole(f, u, up, phi);
  {
    // f = up - phi *u
     Solenoidal_u(f); // div f = 0
  }
}
inline void Add_f_particle(double **u, double **f, const int dim=DIM){
    // !! これを呼ぶ前に、 f を計算しておくこと。
int im;
   // for(int d=0; d<dim; d++){
#pragma omp parallel for schedule(dynamic, 1) private(im) 
    	for(int i=0; i<NX; i++){
	    for(int j=0; j<NY; j++){
		for(int k=0; k<NZ; k++){
		 im=(i*NY*NZ_)+(j*NZ_)+k;
		 u[0][im] += f[0][im];
		 u[1][im] += f[1][im];
		 u[2][im] += f[2][im];
		}
	    }
	    }	
  //  }
}

#endif
