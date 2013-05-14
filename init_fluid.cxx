/*!
  \file init_fluid.cxx
  \brief Initialize fluid velocity fields
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
 */
#include "init_fluid.h"

void Init_zeta_k(double **zeta, double *uk_dc){
    
#pragma omp parallel for schedule(dynamic, 1)
	  for(int i=0; i<NX; i++){
		for(int j=0; j<NY; j++){
		    for(int k=0; k<NZ; k++){
			 int im=(i*NY*NZ_)+(j*NZ_)+k;
			 u[0][im] = 0.e0;
			 u[1][im] = 0.e0;
			 u[2][im] = 0.e0;
		    }
		}
	    }
	U2u_k(u);
	Solenoidal_uk(u);
	U_k2zeta_k(u, zeta, uk_dc);
}


