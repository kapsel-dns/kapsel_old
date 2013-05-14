//
// $Id: solute_rhs.h,v 1.7 2006/02/01 05:11:48 kin Exp $
//
#ifndef SOLUTE_RHS_H
#define SOLUTE_RHS_H

#include "operate_omega.h"
#include "operate_electrolyte.h"
#include "f_particle.h"

extern int NS_source;

extern double *Valency_e;

extern int NP_domain_exponential;
extern int **Sekibun_cell_exponential;

extern double *Total_solute;
extern Value *Concentration;
extern Value *Concentration_rhs0;
extern Value *Concentration_rhs1;

extern Value Surface_normal[DIM];

void Mem_alloc_solute(void);
void Make_phi_exclude(Value &phi
		      ,Particle *p
		      ,const double &dx
		      ,const int &np_domain
		      ,int **sekibun_cell
		      ,const int Nlattice[DIM]
		      );
void Init_solute(Value &Concentration
		 ,Particle *p
		 );
void Solute_solver_rhs_nonlinear_x_single(const Value grad_potential[DIM]
					  ,const Value &concentration_x
					  ,Value solute_flux[DIM]
					  ,const double &valency_e
					  ,const double &Onsager_coeff
					  ,Value &omega_rhs
					  );
void Add_advection_flux(Value solute_flux[DIM]
			,const Value u_solvent[DIM]
			,const Value &concentration_x
			);
void Solute_impermeability(const Particle *p
			   ,Value solute_flux_x[DIM]
			   ,const Value surface_normal[DIM]
			   );
double Count_single_solute(const Value &conc_x
			  ,const Value &phi_p 
			  );
void Rescale_solute(double *rescale_factor
		    ,const double *total_solute
		    ,Value *conc_k
		    ,Particle *p
		    ,Value &phi_p // working memory
		    ,Value &conc_x // working memory
		    );
void Diffusion_flux_single(Value diff_flux_x[DIM]
			    ,const Value &conc_k
			    ,const double &onsager_coeff
			    ,Value &dmy_value // working memory
			   );
inline void Count_solute_each(double *n_solute
			      ,const Value *conc_k
			       ,Particle *p
			       ,Value &phi // working memory
			       ,Value &dmy_value // working memory
			       ){
    Reset_phi(phi);
    Make_phi_particle(phi, p);
    
    for(int n=0;n < N_spec;n++){
	A_k2a_out(conc_k[n], dmy_value);
	n_solute[n] = Count_single_solute(dmy_value, phi);
    }
}
#endif
