//
// $Id: particle_solver.h,v 1.2 2006/10/16 19:05:02 nakayama Exp $
//
#ifndef PARTICLE_SOLVER_H
#define PARTICLE_SOLVER_H

#include <math.h>
#include "variable.h"
#include "input.h"
#include "md_force.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void MD_solver_position_Euler(Particle *p, const CTime &jikan);
void MD_solver_position_AB2(Particle *p, const CTime &jikan);
void MD_solver_velocity_Euler(Particle *p, const CTime &jikan);
void MD_solver_velocity_Euler_hydro(Particle *p, const CTime &jikan);
void MD_solver_velocity_AB2_hydro(Particle *p, const CTime &jikan);
void MD_solver_position_Euler_OBL(Particle *p, const CTime &jikan);
void MD_solver_position_AB2_OBL(Particle *p, const CTime &jikan);
void MD_solver_velocity_Euler_OBL(Particle *p, const CTime &jikan);
void MD_solver_velocity_AB2_hydro_OBL(Particle *p, const CTime &jikan);

inline void Force(Particle *p){
    
    if(LJ_truncate >= 0){
	Calc_f_Lennard_Jones(p);
    }
    
    if(G != 0.0){
	Add_f_gravity(p);
    }
    if(SW_PT == chain){
	Calc_anharmonic_force_chain(p);
    }
}

inline void Force_OBL(Particle *p){
    
    dev_shear_stress_lj = 0.0;

    if(LJ_truncate >= 0){
	double dummy_lj = Calc_f_Lennard_Jones_OBL(p);
	dev_shear_stress_lj += dummy_lj;
    }
    
    if(G != 0.0){
	Add_f_gravity(p);
    }
    if(SW_PT == chain){
	dev_shear_stress_lj +=
	    Calc_anharmonic_force_chain_OBL(p);
    }
    dev_shear_stress_lj *= Ivolume;
}
#endif
