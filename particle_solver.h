//
// $Id: particle_solver.h,v 1.22 2006/04/10 11:32:21 kin Exp $
//
#ifndef PARTICLE_SOLVER_H
#define PARTICLE_SOLVER_H

#include <math.h>
#include "variable.h"
#include "input.h"
#include "md_force.h"

void MD_solver_position_Euler(Particle *p, const CTime &jikan);
void MD_solver_position_CN(Particle *p, const CTime &jikan);
void MD_solver_position_CN_shear(Particle *p, const CTime &jikan);
void MD_solver_position_AB2(Particle *p, const CTime &jikan);
void MD_solver_velocity_Euler(Particle *p, const CTime &jikan);
void MD_solver_velocity_AB2(Particle *p, const CTime &jikan);
void MD_solver_velocity_Euler_hydro(Particle *p, const CTime &jikan);
void MD_solver_velocity_AB2_hydro(Particle *p, const CTime &jikan);
void MD_solver_velocity_sAB2_nohydro(Particle *p, const CTime &jikan);

inline void Reset_particle_velocity(Particle *p){
  for(int n=0; n< Particle_Number; n++){
    for(int d=0; d<DIM; d++){
      p[n].v[d] = 0.;
      p[n].v_old[d] = 0.;

      p[n].omega[d] = 0.;
      p[n].omega_old[d] = 0.;
    }
  }
}
inline void Reset_particle_force_current(Particle *p){
    for(int n=0; n< Particle_Number; n++){
	for(int d=0; d<DIM; d++){  
	  p[n].f_hydro[d] = 0.;
	  p[n].f_hydro1[d] = 0.;
	  p[n].fr[d] = 0.;
	  p[n].fv[d] = 0.;
	  p[n].f_collison[d] = 0.0;
	  
	  p[n].torque_hydro[d] = 0.;
	  p[n].torque_hydro1[d] = 0.;
	  p[n].torquer[d] = 0.;
	  p[n].torquev[d] = 0.;
	}
    }
}
inline void Reset_particle_force_previous(Particle *p){
    for(int n=0; n< Particle_Number; n++){
	for(int d=0; d<DIM; d++){
	  p[n].fr_previous[d] = 0.;
	  p[n].fv_previous[d] = 0.;
	  p[n].f_collison_previous[d] = 0.0;

	  p[n].torquer_previous[d] = 0.;
	  p[n].torquev_previous[d] = 0.;
	}
    }
}
inline void Reset_particle_force(Particle *p){
  Reset_particle_force_previous(p);
  Reset_particle_force_current(p);
}

inline void Force_collison(Particle *p){
    {
	if(SW_EQ == Shear_Navier_Stokes){
	    Calc_f_Lennard_Jones_shear_hydro_cap(p, Srate_depend_LJ_cap);
	    Calc_f_Lennard_Jones_shear_wall(p);
	    //Calc_f_Lennard_Jones_shear_wall_cap(p, Srate_depend_LJ_cap);
	    {
		char line[1<<10];
		{
		    int d=1;
		    sprintf(line,"%g %g %g %g %g %g#LJ"
			    ,SIGMA
			    ,Min_rij
			    ,Min_rij_wall
			    ,Max_force
			    ,Max_force_wall
			    ,Collison_time_shear_hydro(p)
			    );
		    sprintf(line,"%s %d:dia",line, d++);
		    sprintf(line,"%s %d:rij",line, d++);
		    sprintf(line,"%s %d:rij_wall",line, d++);
		    sprintf(line,"%s %d:fmax",line, d++);
		    sprintf(line,"%s %d:fmax_wall",line, d++);
		    sprintf(line,"%s %d:t_collision",line, d++);
		}
		fprintf(stdout, "%s\n", line);
		//fflush(stdout);
	    }
	}else {
	    //Calc_f_Lennard_Jones_shear(p);
	    Calc_f_Lennard_Jones_shear_cap(p);
	}
    }
    {
	Calc_f_Lennard_Jones(p);
    }
}
inline void Force(Particle *p){
  dev_shear_stress_total = 0.0;
  dev_shear_stress_lub = 0.0;
  dev_shear_stress_lj = 0.0;

  if(LJ_truncate >= 0){
      if(SW_EQ == Shear_Navier_Stokes){
	  if(HYDRO_int > 0){
	    dev_shear_stress_lj += Calc_f_Lennard_Jones_shear_hydro_cap(p, Srate_depend_LJ_cap);
	    
	    dev_shear_stress_lub += Calc_f_hydro_lubrication_shear_hydro(p);
	    Calc_f_Lennard_Jones_shear_wall(p);
	    //Calc_f_Lennard_Jones_shear_wall_cap(p, Srate_depend_LJ_cap);
	    {
	      char line[1<<10];
	      {
		int d=1;
		sprintf(line,"%g %g %g %g %g %g#LJ"
			,SIGMA
			,Min_rij
			,Min_rij_wall
			,Max_force
			,Max_force_wall
			,Collison_time_shear_hydro(p)
			);
		sprintf(line,"%s %d:dia",line, d++);
		sprintf(line,"%s %d:rij",line, d++);
		sprintf(line,"%s %d:rij_wall",line, d++);
		sprintf(line,"%s %d:fmax",line, d++);
		sprintf(line,"%s %d:fmax_wall",line, d++);
		sprintf(line,"%s %d:t_collision",line, d++);
	      }
	      fprintf(stdout, "%s\n", line);
	      //fflush(stdout);
	    }

	  }else {
	    //dev_shear_stress_lj += Calc_f_Lennard_Jones_shear(p);
	    dev_shear_stress_lj += Calc_f_Lennard_Jones_shear_cap(p);
	    if(HYDRO_int < 0){
	      dev_shear_stress_lub += Calc_f_hydro_lubrication_shear(p);
	    }
	  }
      }else{
	  Calc_f_Lennard_Jones(p);
      }
  }
  if(G != 0.0){
      Add_f_gravity(p);
  }
  if(0){
    Calc_f_contanct_nonslip(p);
  }
  if(0){
    Calc_f_Coulomb_friction(p);
  }
  if(kBT > 0.0 && Particle_Number == 1 
     && !(SW_EQ == Electrolyte || SW_EQ == Two_fluid)
     ){
      Calc_harmonic_force(p);
  }
  if(SW_EQ == Shear_Navier_Stokes){
      dev_shear_stress_lub *= Ivolume;
      dev_shear_stress_lj *= Ivolume;
      dev_shear_stress_total = dev_shear_stress_lub + dev_shear_stress_lj;
  }
}

#endif
