//
// $Id: operate_electrolyte.h,v 1.1 2006/06/27 18:41:28 nakayama Exp $
//
#ifndef OPERATE_ELECTROLYTE_H
#define OPERATE_ELECTROLYTE_H

#include "fluid_solver.h"
#include "solute_rhs.h"

///////////////////////////////////////
extern double Bjerrum_length;
extern double Surface_ion_number;
extern double Counterion_number;
extern double *Valency;
extern double *Onsager_coeff;

const double TOL = 1.e-4;
///////////////////////////////////////


extern double *Total_solute;

void Init_rho_ion(double **Concentration, Particle *p, CTime &jikan);

void Mem_alloc_charge(void);
void Calc_free_energy_PB(double **conc_k
			  ,Particle *p
			  ,double *free_energy
			  ,double *phi // working memory
			  ,double *charge_density // working memory
			  ,double *dmy_value // working memory
			  ,const CTime &jikan
			  );
void Make_phi_qq_particle(double *phi
			   ,double *surface
			   ,Particle *p);
void Make_phi_qq_fixed_particle(double *phi
			   ,double *surface
			   ,Particle *p);
void Make_Coulomb_force_x_on_fluid(double **force
				    ,Particle *p
				    ,double **conc_k
				    ,double *charge_density // working memory
				    ,double *potential // working memory
				    ,const CTime &jikan
				    );

void Conc_k2charge_field(Particle *p
			 ,double **conc_k
			 ,double *charge_density
			 ,double *phi // working memory
			 ,double *dmy_value // working memory
			 );

void Charge_field_k2Coulomb_potential_k_PBC(double *potential);

#endif

