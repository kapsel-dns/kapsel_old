//
// $Id: operate_electrolyte.h,v 1.10 2006/05/15 09:43:39 kin Exp $
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

void Init_rho_ion(Value *Concentration, Particle *p, CTime &jikan);

void Mem_alloc_charge(void);
void Calc_free_energy_PB(const Value *conc_k
			  ,Particle *p
			  ,double *free_energy
			  ,Value &phi // working memory
			  ,Value &charge_density // working memory
			  ,Value &dmy_value // working memory
			  ,const CTime &jikan
			  );
void Make_phi_qq_particle(Value &phi
			   ,Value &surface
			   ,const Particle *p);
void Make_phi_qq_fixed_particle(Value &phi
			   ,Value &surface
			   ,const Particle *p);
void Calc_Coulomb_potential_k_PBC(const Particle *p
				   ,const Value *conc_k
				   ,Value &potential
				   ,Value &charge_density
				   ,Value &phi // working memory
				   );
double Total_ion_charge(const Value *conc_k
			 ,Particle *p
			 ,Value &phi // working memory
			 ,Value &dmy_value // working memory
			 );
void Make_Coulomb_force_x_on_fluid(Value force[DIM]
				    ,const Particle *p
				    ,const Value *conc_k
				    ,Value &charge_density // working memory
				    ,Value &potential // working memory
				    ,const CTime &jikan
				    );

void Conc_k2charge_field(const Particle *p
			 ,const Value *conc_k
			 ,Value &charge_density
			 ,Value &phi // working memory
			 ,Value &dmy_value // working memory
			 );
void Charge_field_k2Coulomb_potential_k_PBC(Value &potential);

#endif

