//
// $Id: variable.h,v 1.7 2005/12/01 09:31:43 nakayama Exp $
//
#ifndef VARIABLE_H
#define VARIABLE_H

#include "parameter_define.h"

// omega,ux,uy,phi, phi_up など場の変数を格納する
typedef double *** Value;
typedef int *** Value_int;

typedef struct CTime{
    int ts; // time step
    double time; // time
    double dt_fluid;  // time increment
    double hdt_fluid; // 1/2 * dt 
    double dt_md;  // time increment
    double hdt_md; // 1/2 * dt 
} CTime;

typedef struct Particle{
    double eff_mass_ratio;
    int spec;
    double x[DIM];

    double v[DIM];
    double v_old[DIM];

    double f_hydro[DIM];
    double f_hydro1[DIM];
  //double f_hydro_previous[DIM];
    double fr[DIM];
    double fr_previous[DIM];
    double fv[DIM];
    double fv_previous[DIM];

    double f_collison[DIM];
    double f_collison_previous[DIM];

    double omega[DIM];
    double omega_old[DIM];

    double torque_hydro[DIM];
    double torque_hydro1[DIM];
  //double torque_hydro_previous[DIM];
    double torquer[DIM];
    double torquer_previous[DIM];
    double torquev[DIM];
    double torquev_previous[DIM];
} Particle;

typedef struct Index_range{
    int istart;
    int iend;
    int jstart;
    int jend;
    int kstart;
    int kend;
} Index_range;
#endif
