//
// $Id: md_force.h,v 1.16 2005/09/13 06:38:58 nakayama Exp $
//
#ifndef MD_FORCE_H
#define MD_FORCE_H

#include <assert.h> 
#include "variable.h"
#include "input.h"
#include "interaction.h"
#include "make_phi.h"

extern double Min_rij;
extern double Max_force;
extern double Min_rij_wall;
extern double Max_force_wall;

enum Particle_BC {
  PBC_particle
  ,Lees_Edwards
  ,Shear_hydro
};
double Calc_f_Lennard_Jones_shear_cap_primitive(Particle *p
				       ,void (*distance0_func)(const double *x1,const double *x2,double &r12,double *x12)
				      ,const double cap
				       );
double Calc_f_hydro_lubrication_shear_primitive(Particle *p
						,void (*distance0_func)(const double *x1,const double *x2,double &r12,double *x12)
						,const Particle_BC SW_BC);
double Collison_time(Particle *p
		     ,void (*distance0_func)(const double *x1,const double *x2,double &r12,double *x12)
		     ,const Particle_BC sw_bc
		     );

void Add_f_gravity(Particle *p);
void Calc_f_hydro_draining(Particle *p, const Value *fp, const CTime &jikan);
void Calc_f_hydro_draining_shear(Particle *p, const Value *fp, const CTime &jikan);
void Calc_f_hydro_correct(Particle *p, const Value *fp, const CTime &jikan);
void Calc_particle_velocity(const Value &phi, const Value *up, Particle *p);

void Calc_f_contanct_nonslip(Particle *p);
void Calc_f_Coulomb_friction(Particle *p);
void Calc_f_Lennard_Jones_wall_cap(Particle *p
				   ,const int &dwall
				   ,const double &a_radius_wall
				   ,const double cap = 5.e2
				   );

/////
inline double Collison_time_shear_hydro(Particle *p){
  return
    Collison_time(p,Distance0_shear_hydro,Shear_hydro);
}
inline void Calc_f_hydro_lubrication(Particle *p){
  Calc_f_hydro_lubrication_shear_primitive(p,Distance0,PBC_particle);
}
inline double Calc_f_hydro_lubrication_shear(Particle *p){
  return
    Calc_f_hydro_lubrication_shear_primitive(p,Distance0_shear,Lees_Edwards);
}
inline double Calc_f_hydro_lubrication_shear_hydro(Particle *p){
  return
    Calc_f_hydro_lubrication_shear_primitive(p,Distance0_shear_hydro,Shear_hydro);
}
/////
inline void Calc_f_Lennard_Jones_cap(Particle *p
				     ,const double cap = 5.e2){
  Calc_f_Lennard_Jones_shear_cap_primitive(p,Distance0,cap);
}
inline double Calc_f_Lennard_Jones_shear_cap(Particle *p
					     ,const double cap = 5.e2){
  return 
    Calc_f_Lennard_Jones_shear_cap_primitive(p,Distance0_shear,cap);
}

inline double Calc_f_Lennard_Jones_shear_hydro_cap(Particle *p
						   ,const double cap = 5.e2){
  return 
    Calc_f_Lennard_Jones_shear_cap_primitive(p, Distance0_shear_hydro, cap);
}
inline void Calc_f_Lennard_Jones(Particle *p){
  Calc_f_Lennard_Jones_shear_cap_primitive(p,Distance0,DBL_MAX);
}
inline double Calc_f_Lennard_Jones_shear(Particle *p){
  return 
    Calc_f_Lennard_Jones_shear_cap_primitive(p,Distance0_shear,DBL_MAX);
}
inline double Calc_f_Lennard_Jones_shear_hydro(Particle *p){
  return 
    Calc_f_Lennard_Jones_shear_cap_primitive(p, Distance0_shear_hydro, DBL_MAX);
}
/////
inline void Calc_f_Lennard_Jones_wall(Particle *p
				      ,const int &dwall
				      ,const double &a_radius_wall
				      ){
  Calc_f_Lennard_Jones_wall_cap(p,dwall,a_radius_wall,DBL_MAX);
}
inline void Calc_f_Lennard_Jones_zwall(Particle *p){
  Calc_f_Lennard_Jones_wall_cap(p,2,Radius_wall,DBL_MAX);
}
inline void Calc_f_Lennard_Jones_shear_wall(Particle *p){
  Calc_f_Lennard_Jones_wall_cap(p,1,0.,DBL_MAX);
}
inline void Calc_f_Lennard_Jones_shear_wall_cap(Particle *p
						,const double cap = 5.e2
						){
  Calc_f_Lennard_Jones_wall_cap(p,1,0.,cap);
}
/////

const double harmonic_spring_cst=1.0; 
inline void Calc_harmonic_force(Particle *p){
    const double harmonic_eq[DIM]={HL_particle[0], HL_particle[1], HL_particle[2]};
    int n=0;
    for(int d=0; d < DIM; d++ ){ 
	p[n].fr[d] += -harmonic_spring_cst * (p[n].x[d]- harmonic_eq[d]);
    }
} 
inline void Calc_harmonic_energy(const Particle *p, double energy[DIM]){
    const double harmonic_eq[DIM]={HL_particle[0], HL_particle[1], HL_particle[2]};

    const double dmy = harmonic_spring_cst * 0.5;
    int n=0;
    for(int d=0; d < DIM; d++ ){ 
	energy[d] = dmy * SQ(p[n].x[d]- harmonic_eq[d]);
	//fprintf(stderr, "%g %g\n", p[n].x[d], harmonic_eq[d]);
    }
} 
#endif
