/*!
  \file particle_solver.h
  \brief Solver routines for particle position and velocity (header file)
  \author Y. Nakayama
  \date 2006/10/16
  \version 1.2
  \todo documentation
 */
#ifndef PARTICLE_SOLVER_H
#define PARTICLE_SOLVER_H

#include <math.h>
#include "variable.h"
#include "input.h"
#include "md_force.h"
#include "rigid_body.h"
#include "rigid.h"

#ifdef _OPENMP
#include <omp.h>
#endif

enum ITER {start_iter, new_iter, reset_iter, end_iter};

/*!
  \brief Update particle positions and orientations using the Euler
  forward method
  \details \f{align*}{
  \vec{R}_i^{n+1} &= \vec{R}_i^{n} + h \vec{V}_i^{n} \\
  \qtn{q}_i^{n+1} &= \qtn{q}_i^{n} + \frac{h}{2}(0,\vec{\omega}_i^{n})\circ\qtn{q}_i^{n}
  \f}
  \param[in,out] p particle data
  \param[in] jikan time data
 */
void MD_solver_position_Euler(Particle *p, const CTime &jikan);

/*!
  \brief Update particle positions and orientations using a
  second-order Adams-Bashforth scheme
  \details \f{align*}{
  \vec{R}_i^{n+1} &= \vec{R}_i^{n} + \frac{h}{2}\left(3\vec{V}_i^{n} -
  \vec{V}_i^{n-1}\right) \\
  \qtn{q}_i^{n+1} &= \qtn{q}_i^n +
  \frac{h}{4}\qtn{q}_i^{n}\circ(0,3\vec{\omega}_i^{\prime
  n}-\vec{\omega}_i^{\prime n-1})
  \f}
  Note that addition of angular velocity vectors only makes sense in
  body (primed) coordinates.
  \param[in,out] p particle data
  \param[in] jikan time data
 */
void MD_solver_position_AB2(Particle *p, const CTime &jikan);

/*!
  \brief  Update particle velocities using Euler forward method
  \details Saves data for old velocities and forces
  \param[in,out] p particle data
  \param[in] jikan time data
 */
void MD_solver_velocity_Euler(Particle *p, const CTime &jikan);
/*!
  \brief Update particle velocities using a second-order
  Adams-Bashforth scheme
  \details Saves data for old velocities and forces
  \param[in,out] p particle data
  \param[in] jikan time data
 */
void MD_solver_velocity_AB2_hydro(Particle *p, const CTime &jikan);

/*!
  \brief Update particle velocities for swimming particles 
  \details Automatically chooses appropriate Euler/Adams-Bashforth
  scheme depending on the jikan step. Part of iterative solution for
  particle velocities. Only slip force changes, other quantities
  should only be computed at the first iteration.
 */
void MD_solver_velocity_slip_iter(Particle *p, const CTime &jikan, const ITER &iter_flag);


// Oblique coordinates
/*!
  \brief Update particle positions and orientations using the Euler
  forward method for sheared systems with Lees-Edwards PBC
 */
void MD_solver_position_Euler_OBL(Particle *p, const CTime &jikan);
/*!
  \brief Update particle positions and orientations using a
  second-order Adams-Bashforth scheme for sheared systems with Lees-Edwards PBC
 */
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
