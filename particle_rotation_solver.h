/*!
  \file particle_rotation_solver.h
  \brief Solver routines for particle orientations
  \author J. Molina
  \date 2013/05/15
  \version 1.0
 */
#ifndef PARTICLE_ROTATION_SOLVER_H
#define PARTICLE_ROTATION_SOLVER_H

#include "variable.h"
#include "rigid_body.h"
//Orientation solvers
//First-order Euler
inline void MD_solver_orientation_Euler(Particle &p, const double &dt){
  quaternion dqdt;
  qtn_init(p.q_old, p.q);
  qdot(dqdt, p.q, p.omega, SPACE_FRAME);
  qtn_add(p.q, dqdt, dt);
  qtn_normalize(p.q);
}

//Simo & Wong second-order scheme
inline void MD_solver_orientation_AB2(Particle &p, const double &dt){
  double wb[DIM];
  double wb_old[DIM];
  //only add angular velocity vectors in body coordinates !
  rigid_body_rotation(wb, p.omega, p.q, SPACE2BODY);
  rigid_body_rotation(wb_old, p.omega_old, p.q_old, SPACE2BODY);
  for(int d = 0; d < DIM; d++){
    wb[d] = 3.0*wb[d] - wb_old[d];
  }
  
  quaternion dqdt;
  qtn_init(p.q_old, p.q);
  qdot(dqdt, p.q, wb, BODY_FRAME);
  qtn_add(p.q, dqdt, dt);
  qtn_normalize(p.q);
}

// Samuel Buss' second-order scheme
// untested
inline void MD_solver_orientation_SB2(Particle &p, const double &dt){
  double wb[DIM];
  for(int d = 0; d < DIM; d++){
    wb[d] = p.omega[d] + dt/2.0 * IMOI[p.spec]*p.torque_hydro[d];
  }
  
  quaternion dqdt;
  qtn_init(p.q_old, p.q);
  qdot(dqdt, p.q, wb, SPACE_FRAME);
  qtn_add(p.q, dqdt, dt);
  qtn_normalize(p.q);
}

#endif
