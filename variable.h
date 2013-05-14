/*!
  \file variable.h
  \brief Defines the global structs (CTime, Particle, Index_range)
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
 */
#ifndef VARIABLE_H
#define VARIABLE_H

#include "parameter_define.h"
#include "quaternion.h"

// omega,ux,uy,phi, phi_up など場の変数を格納する
typedef double ** *Value;
typedef int ** *Value_int;

typedef struct CTime {
  int ts; // time step
  double time; // time
  double dt_fluid;  // time increment
  double hdt_fluid; // 1/2 * dt
  double dt_md;  // time increment
  double hdt_md; // 1/2 * dt
} CTime;

typedef struct Particle {
  int spec;
  double x[DIM];
  double x_previous[DIM];
  double x_nopbc[DIM];

  double v[DIM];
  double v_old[DIM];
  double v_slip[DIM];

  double f_hydro[DIM];
  double f_hydro_previous[DIM];
  double f_hydro1[DIM];
  double f_slip[DIM];
  double f_slip_previous[DIM];

  double fr[DIM];
  double fr_previous[DIM];

  double omega[DIM];
  double omega_old[DIM];
  double omega_slip[DIM];

  double torque_hydro[DIM];
  double torque_hydro_previous[DIM];
  double torque_hydro1[DIM];
  double torque_slip[DIM];
  double torque_slip_previous[DIM];

  double momentum_depend_fr[DIM];

  double QR[DIM][DIM];
  double QR_old[DIM][DIM];

  quaternion q;
  quaternion q_old;

  //working memory for iterative slip implementation
  double mass;                       //fluid particle mass
  double mass_center[DIM];           //center of mass of fluid particle
  double inertia[DIM][DIM];          //moment of inertia of fluid particle

  double surface_mass;               //surface fluid mass
  double surface_mass_center[DIM];   //surface fluid center of mass
  double surface_inertia[DIM][DIM];  //surface fluid moment of inertia

  double surface_dv[DIM];            //momentum change due to slip
  double surface_dw[DIM];        //ang. momentum change due to slip
} Particle;

typedef struct Index_range {
  int istart;
  int iend;
  int jstart;
  int jend;
  int kstart;
  int kend;
} Index_range;


#endif
