//
// $Id: f_particle.cxx,v 1.2 2004/11/24 05:02:38 nakayama Exp $
//
#include "f_particle.h"

Value *f_particle;

void Mem_alloc_f_particle(void){
  f_particle = (Value *) malloc(sizeof(Value *) * DIM);
  for(int d=0;d<DIM;d++){
    f_particle[d] = alloc_3d_double(NX, NY, NZ_);
  }

}

 
