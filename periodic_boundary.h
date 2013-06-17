#ifndef PERIODIC_BOUNDARY_H
#define PERIODIC_BOUNDARY_H

#include "input.h"
/*!
  \brief Enforce pbc on position 
  \param[in,out] x position
 */
inline void PBC(double *x){
  for(int d = 0; d < DIM; d++){
    x[d] = fmod(x[d] + L_particle[d], L_particle[d]);
    assert(x[d] >= 0);
    assert(x[d] < L[d]);
  }
}

/*!
  \brief Enforce Lees-Edwards pbc on position
  \param[in,out] x position
  \param[out] delta_vx change in x-velocity due to boundary crossing
 */
inline int PBC_OBL(double *x, double &delta_vx){
  double signY = x[1];
  x[1] = fmod(x[1] + L_particle[1], L_particle[1]);
  signY -= x[1];
  int sign = (int) signY;
  if(!(sign == 0)){
    sign = sign / abs(sign);
  }
  
  
  x[0] -= (double)sign * degree_oblique * L_particle[1];
  x[0] = fmod(x[0] + L_particle[0], L_particle[0]);
  x[2] = fmod(x[2] + L_particle[2], L_particle[2]);
  for(int d = 0; d < DIM; d++){
    assert(x[d] >= 0);
    assert(x[d] < L[d]);
  }
  
  delta_vx = -(double)sign * Shear_rate_eff * L_particle[1];
  return sign;
}
#endif
