//
// $Id: init_fluid.cxx,v 1.5 2005/04/21 05:26:10 nakayama Exp $
//

#include "init_fluid.h"

void Init_zeta_k(Value *zeta, double *uk_dc){

  SRA(GIVEN_SEED,10);
  if(SW_EQ == Shear_Navier_Stokes){
    for(int d=0;d<DIM;d++){
      for(int i=0; i<NX; i++){
	for(int j=1; j<NY_shear; j++){
	  for(int k=0; k<NZ; k++){
	    if(d==1){ // sin
	      double dmy = 0.e0;
	      //double dmy = RA() * 1.e-1;
	      u[d][i][j][k] = dmy;
	      u[d][i][NY-j][k] = -dmy;
	    }else { // cos
	      double dmy = 0.e0;
	      //double dmy = RA() * 1.e-1;
	      u[d][i][j][k] = dmy;
	      u[d][i][NY-j][k] = dmy;
	    }
	  }
	}
      }
    }
    for(int d=0;d<DIM;d++){
      for(int i=0; i<NX; i++){
	for(int k=0; k<NZ; k++){
	  if(d==1){ // sin
	    u[d][i][NY_shear][k]
	    = u[d][i][0][k] = 0.e0;
	  }else { // cos
	    double dmy = 0.e0;
	    //double dmy = RA() * 1.e-1;
	    u[d][i][0][k] = dmy;
	    u[d][i][NY_shear][k] = dmy;
	  }
	}
      }
    }
    U2u_k(u);
    {
      Symmetrize_fourier(u[0]);
      Antisymmetrize_fourier(u[1]);
      Symmetrize_fourier(u[2]);
    }
    Solenoidal_uk(u);
    U_k2zeta_k(u, zeta, uk_dc);
  }else {
    for(int d=0;d<DIM;d++){
      for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	  for(int k=0; k<NZ; k++){
	    u[d][i][j][k] = 0.e0;
	    //u[d][i][j][k] = 1.e1 * RA();
	  }
	}
      }
    }
    U2u_k(u);
    Solenoidal_uk(u);
    U_k2zeta_k(u, zeta, uk_dc);
  }
}

