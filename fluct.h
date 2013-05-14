//
// $Id: fluct.h,v 1.1 2006/06/27 18:41:28 nakayama Exp $
//
#ifndef FLUCT_H
#define FLUCT_H

#include <stdio.h>
#include <math.h>
#include "input.h"
#include "variable.h"
#include "operate_omega.h"
#include "make_phi.h"
#include "particle_solver.h"
//#include "SFMT.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void init_genrand(unsigned long s);
extern double genrand_real3(void);

#ifdef __cplusplus
}
#endif

void Add_random_force_thermostat(Particle *p, const CTime &jikan);

inline void MT_seed(const int &SW_seed, const unsigned long &seed){
    if(SW_seed == RANDOM_SEED){
      init_genrand(time(NULL));
      fprintf(stderr, "# MT_seed= time(NULL)\n");
    }else if(SW_seed == GIVEN_SEED){
      init_genrand(seed);
      fprintf(stderr, "# MT_seed= %lu\n",seed);
    }else {
      fprintf(stderr, "MT_seed(): invalid SW_seed.\n");
      exit_job(EXIT_FAILURE);
    }

}

inline void Gauss2(double random[2]){
  static double x1,x2;
  static double rsq, factor;
  do{
    x1 = 2.*genrand_real3()-1.;
    x2 = 2.*genrand_real3()-1.;

    rsq = x1 * x1 + x2 * x2;
  }while( rsq >= 1.0 || rsq == 0.0);
  factor = sqrt(-2.*log(rsq)/rsq);
  random[0] = factor * x1;
  random[1] = factor * x2;
}

inline void Force_random_walk(Particle *p){
  Calc_f_Lennard_Jones(p);
}

inline void Random_Walk(Particle *p
			,const double dr_factor = 0.5e-1*.5
			){
  const double dr = RADIUS * dr_factor;

  for(int n=0; n<Particle_Number; n++){
    double dmy[4];
    Gauss2(dmy);
    Gauss2(dmy+2);

    double h = 0.;
    for(int d=0; d<DIM; d++){  
      h += SQ(dmy[d]);
    }
    h = sqrt(h);
    if(h> 0.){    
      h = dr/h;
    }
    for(int d=0; d<DIM; d++){      
	  p[n].x[d] += h * dmy[d];
      p[n].x[d] -= floor(p[n].x[d]*iL_particle[d])*L_particle[d];
      
	  assert(p[n].x[d] >= 0.);
      assert(p[n].x[d] < L_particle[d]);
	}
  }
}
inline void Steepest_descent(Particle *p
			,const double dr_factor = 0.5e-1
			){
  const double dr = RADIUS * dr_factor;

  Force_random_walk(p);
  for(int n=0; n<Particle_Number; n++){
    double h = 0.;
    for(int d=0; d<DIM; d++){  
      h += SQ(p[n].fr[d]);
    }
    h = sqrt(h);
    if(h> 0.){    
      h = dr/h;
    }
    for(int d=0; d<DIM; d++){      
	  p[n].x[d] += h * p[n].fr[d];
      p[n].x[d] -= floor(p[n].x[d]*iL_particle[d])*L_particle[d];
   
      assert(p[n].x[d] >= 0.);
      assert(p[n].x[d] < L_particle[d]);
	}
  }
}
#endif
