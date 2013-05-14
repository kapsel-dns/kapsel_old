//
// $Id: f_particle.cxx,v 1.1 2006/06/27 18:41:28 nakayama Exp $
//
#include "f_particle.h"

double **f_particle;

void Mem_alloc_f_particle(void){
  f_particle = (double **) malloc(sizeof(double *) * DIM);
  for(int d=0;d<DIM;d++){
    f_particle[d] = alloc_1d_double(NX*NY*NZ_);
  }
}

 
