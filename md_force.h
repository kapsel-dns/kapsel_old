//
// $Id: md_force.h,v 1.1 2006/06/27 18:41:28 nakayama Exp $
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
extern double *Hydro_force;

enum Particle_BC {
  PBC_particle
  ,Lees_Edwards
  ,Shear_hydro
};

void Add_f_gravity(Particle *p);
void Calc_f_hydro_correct_precision(Particle *p, double **fp, const CTime &jikan);
void Calc_f_hydro_correct_precision_OBL(Particle *p, double **fp, const CTime &jikan);
double Calc_f_Lennard_Jones_shear_cap_primitive_lnk(Particle *p
				       ,void (*distance0_func)(const double *x1,const double *x2,double &r12,double *x12)
				      ,const double cap
				       );
double Calc_f_Lennard_Jones_shear_cap_primitive(Particle *p
				       ,void (*distance0_func)(const double *x1,const double *x2,double &r12,double *x12)
				      ,const double cap
				       );
inline void Calc_f_Lennard_Jones(Particle *p){
//  Calc_f_Lennard_Jones_shear_cap_primitive(p,Distance0,DBL_MAX);
    Calc_f_Lennard_Jones_shear_cap_primitive_lnk(p, Distance0, DBL_MAX);
}
inline double Calc_f_Lennard_Jones_OBL(Particle *p){
    return Calc_f_Lennard_Jones_shear_cap_primitive(p, Distance0_OBL, DBL_MAX);
}
inline void Calc_anharmonic_force_chain(Particle *p){
    double anharmonic_spring_cst=30.*EPSILON/SQ(SIGMA);
    const double R0=1.5*SIGMA;
    const double iR0=1./R0;
    const double iR02=SQ(iR0);
    
    int rest_Chain_Number[Component_Number];
    for (int i = 0; i < Component_Number; i++) {
	rest_Chain_Number[i] = Chain_Numbers[i];
//	fprintf(stderr,"### rest_p_num %d %d\n", i, rest_Chain_Number[i]);
    }
    
    int n = 0;

    while (n < Particle_Number - 1) {
	for (int i = 0; i < Component_Number; i++) {
	    if (rest_Chain_Number[i] > 0) {
		for (int j = 0; j < Beads_Numbers[i]; j++) {
		    //fprintf(stderr,"### anharmonic %d\n", j);
		    if (j < Beads_Numbers[i] - 1) {
			int m = n + 1;
			double dmy_r1[3];
			for(int d=0; d<DIM ;d++){
			    dmy_r1[d]=p[n].x[d]-p[m].x[d];
			    dmy_r1[d] -= (double)Nint(dmy_r1[d]/L_particle[d])*L_particle[d];
			}
			double dm_r1=SQ(dmy_r1[0])+SQ(dmy_r1[1])+SQ(dmy_r1[2]);
			double dm1=1./(1. - dm_r1*iR02);
			
			if(dm1 < 0.0){
			    fprintf(stderr,"### anharmonic error\n");
			}
			
			for(int d=0; d<DIM; d++){
			    double dmy=dm1*dmy_r1[d];
			    p[n].fr[d] += -anharmonic_spring_cst*dmy;
			    p[m].fr[d] += anharmonic_spring_cst*dmy;
			}
		    }
		    n++;
		}
		rest_Chain_Number[i]--;
	    }
	}
    }
}

inline double Calc_anharmonic_force_chain_OBL(Particle *p){
    double anharmonic_spring_cst=30.*EPSILON/SQ(SIGMA);
    const double R0=1.5*SIGMA;
    const double iR0=1./R0;
    const double iR02=SQ(iR0);
    
    int rest_Chain_Number[Component_Number];
    for (int i = 0; i < Component_Number; i++) {
	rest_Chain_Number[i] = Chain_Numbers[i];
//	fprintf(stderr,"### rest_p_num %d %d\n", i, rest_Chain_Number[i]);
    }
    
    int n = 0;
    double shear_stress = 0.0;
    while (n < Particle_Number - 1) {
	for (int i = 0; i < Component_Number; i++) {
	    if (rest_Chain_Number[i] > 0) {
		for (int j = 0; j < Beads_Numbers[i]; j++) {
		    //fprintf(stderr,"### anharmonic %d\n", j);
		    if (j < Beads_Numbers[i] - 1) {
			int m = n + 1;
			double dmy_r1[DIM];
			double dm_r1 = 0.0;
			Distance0_OBL(p[m].x, p[n].x, dm_r1, dmy_r1);
			/*
			for(int d=0; d<DIM ;d++){
			    dmy_r1[d]=p[n].x[d]-p[m].x[d];
			    dmy_r1[d] -= (double)Nint(dmy_r1[d]/L_particle[d])*L_particle[d];
			}
			double dm_r1=SQ(dmy_r1[0])+SQ(dmy_r1[1])+SQ(dmy_r1[2]);
			*/
			double dm1=1./(1. - dm_r1*iR02);
			
			if(dm1 < 0.0){
			    fprintf(stderr,"### anharmonic error\n");
			}
			
			for(int d=0; d<DIM; d++){
			    double dmy=dm1*dmy_r1[d];
			    p[n].fr[d] += -anharmonic_spring_cst*dmy;
			    p[m].fr[d] += anharmonic_spring_cst*dmy;
			}
			shear_stress +=
			    ((-anharmonic_spring_cst * dm1 * dmy_r1[0])* (dmy_r1[1]));
		    }
		    n++;
		}
		rest_Chain_Number[i]--;
	    }
	}
    }
    return shear_stress;
}
#endif
