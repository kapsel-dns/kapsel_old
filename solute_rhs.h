//
// $Id: solute_rhs.h,v 1.1 2006/06/27 18:41:29 nakayama Exp $
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

extern double **Surface_normal;
extern double **Concentration;
extern double **Concentration_rhs0;
extern double **Concentration_rhs1;

void Mem_alloc_solute(void);
void Solute_solver_rhs_nonlinear_x_single(double **grad_potential
					  ,double *concentration_x
					  ,double **solute_flux
					  ,double &valency_e
					  ,double &Onsager_coeff
					  ,double *omega_rhs
					  );
void Add_advection_flux(double **solute_flux
			,double **u_solvent
			,double *concentration_x
			);
void Solute_impermeability(Particle *p
			   ,double **solute_flux_x
			   ,double **surface_normal
			   );
double Count_single_solute(double *conc_x
			  ,double *phi_p 
			  );
void Rescale_solute(double *rescale_factor
		    ,double *total_solute
		    ,double **conc_k
		    ,Particle *p
		    ,double *phi // working memory
		    ,double *conc_x // working memory
		    );
void Diffusion_flux_single(double **diff_flux_x
			    ,double *conc_k
			    ,double &onsager_coeff
			    ,double *dmy_value // working memory
			   );
inline void Count_solute_each(double *n_solute
			      ,double **conc_k
			       ,Particle *p
			       ,double *phi // working memory
			       ,double *dmy_value // working memory
			       ){
    Reset_phi(phi);
    Make_phi_particle(phi, p);
    
    for(int n=0;n < N_spec;n++){
	A_k2a_out(conc_k[n], dmy_value);
	n_solute[n] = Count_single_solute(dmy_value, phi);
    }
}
#endif
