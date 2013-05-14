//
// $Id: fluct.cxx,v 1.1 2006/06/27 18:41:28 nakayama Exp $
//
#include "fluct.h"

void Add_random_force_thermostat(Particle *p, const CTime &jikan){
  static const double Zeta_drag = 6.* M_PI * ETA * RADIUS;
  static const double Zeta_drag_rot = 8.* M_PI * ETA * POW3(RADIUS);
  
  const double sdv_v = sqrt(Zeta_drag * jikan.dt_md * kBT * alpha_v);
  const double sdv_omega = sqrt(Zeta_drag_rot * jikan.dt_md * kBT * alpha_o);
   
  static double noise_intensity_v = kT_snap_v * sdv_v;
  static double noise_intensity_o = kT_snap_o * sdv_omega;

  for(int n=0; n<Particle_Number; n++){
    double dmy[6];
    Gauss2(dmy);
    Gauss2(dmy+2);
    Gauss2(dmy+4);

    double imass = IMASS[p[n].spec];
    double imoi = IMOI[p[n].spec];
    
    for(int d=0; d<DIM; d++){
      p[n].v[d] +=    dmy[d]*noise_intensity_v*imass;
      if(ROTATION){
      p[n].omega[d] += dmy[d+3]*noise_intensity_o*imoi;
      }
    }
  }
}
