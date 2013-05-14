/*!
  \file f_particle.cxx
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief Compute Hydrodynamic force on particle
 */
#include "f_particle.h"

double **f_particle;

void Mem_alloc_f_particle(void){
  f_particle = (double **) malloc(sizeof(double *) * DIM);
  for(int d=0;d<DIM;d++){
    f_particle[d] = alloc_1d_double(NX*NY*NZ_);
  }
}

 
