/*!
  \file solute_rhs.h
  \brief Routines to compute the terms appearing in the right hand side of the solute advection diffusion equation (header file)
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
 */
#ifndef SOLUTE_RHS_H
#define SOLUTE_RHS_H

#include "operate_omega.h"
#include "operate_electrolyte.h"
#include "f_particle.h"

extern int NS_source;//!< (Unused) flag to enable solute source calculations

extern double *Valency_e;//!< Charge valency for each species


extern double *Total_solute;//!< Total amount of solute for each species

extern double **Surface_normal;//!< Surface normal field at particle interface
extern double **Concentration;//!< Concentration field for each species
extern double **Concentration_rhs0;//!< Auxiliary concentration field to solve solute advection-diffusion equation
extern double **Concentration_rhs1;//!< Auxiliary concentration field to sovle solute-advection-diffusion equation

/*!
  \brief Allocate solute concentration memory
 */
void Mem_alloc_solute(void);

/*!
  \brief Add the non-linear diffusive flux (due to the total charge distribution) to the rhs of the solute advection-diffusion equation for a given species
  \details
  \f[
  -\vec{j}_\alpha(\vec{r}) \longrightarrow
  -\left(\vec{j}_\alpha(\vec{r}) + 
  \Gamma_\alpha Z_\alpha e C_\alpha(\vec{r}) \nabla \Phi^{\text{tot}}\right)
  \f]
  \param[in] grad_potential gradient of the total electrostatic potential (solute + colloid charges) in r-space
  \param[in] concentration_x solute concentration field (r-space)
  \param[in,out] solute_flux updated solute flux, including the non-linear diffusive current due to electrostatic interactions
  \param[in] valency_e valence charges of the given solute species
  \param[in] Onsager_coeff Onsager coefficient of the given solute species
 */
void Solute_solver_rhs_nonlinear_x_single(double **grad_potential
					  ,double *concentration_x
					  ,double **solute_flux
					  ,double &valency_e
					  ,double &Onsager_coeff
					  );
/*!
  \brief Adds the advection flux contribution (due to solvent motion) to the rhs of the advection-diffusion equation for a given solute species
  \details 
  \f[
  -\vec{j}_\alpha(\vec{r}) \longrightarrow
  -\left(\vec{j}_\alpha(\vec{r}) + C_\alpha(\vec{r})\vec{u}(\vec{r})\right)
  \f]
  \param[in,out] solute_flux updated (negative) solute flux
  \param[in] u_solvent solvent velocity field
  \param[in] concentration_x solute concentration field (r-space)
 */
void Add_advection_flux(double **solute_flux
			,double **u_solvent
			,double *concentration_x
			);
/*!
  \brief Enforce the no-penetration condition (with the particle domain) on the (total) solute diffusive flux
  \details
  \f[
  \vec{j}_\alpha^{d} \longrightarrow
  (\tensor{I} - \vec{n}\vec{n}) \vec{j}_{\alpha}^{d}
  \f]
  where \f$\vec{n}(\vec{r})\f$ is the unit surface-normal vector field defined on the particle interface domain.
  \note The advection contribution to the solute flux should \b NOT be included
 */
void Solute_impermeability(Particle *p
			   ,double **solute_flux_x
			   ,double **surface_normal
			   );
/*!
  \brief Returns the total amount of solute (outside the particle domain) of a given species from the concentration field
  \details
  \f[
  N_\alpha^{\text{total}} = \int\vdf{r} (1 - \phi(r)) C_{\alpha}(\vec{r})
  \f]
  \param[in] conc_x concentration field of a single species (r-space)
  \param[in] phi_p smooth particle field
 */
double Count_single_solute(double *conc_x
			  ,double *phi_p 
			  );
/*!
  \brief Rescales the solute concentration of all species to ensure total amount is conserved
  \details
  \f[
  C_\alpha(\vec{r}) \longrightarrow \frac{N_\alpha^{\text{total}}}
  {\int\vdf{r} C_\alpha(\vec{r})} C_\alpha(\vec{r})
  \f]
  \param[out] rescale_factor global scaling factors used to normalize each concentration field
  \param[in] total_solute total (target) amount of solute for each species
  \param[in] conc_k concentration field for each species (k-space)
  \param[in] p particle data
  \param[in,out] phi working memory to compute smooth particle field
  \param[in,out] conc_x working memory to compute real-space concentration field for a single species
 */
void Rescale_solute(double *rescale_factor
		    ,double *total_solute
		    ,double **conc_k
		    ,Particle *p
		    ,double *phi // working memory
		    ,double *conc_x // working memory
		    );

/*!
  \brief Compute the (negative) single ion diffusion flux of a given species
  due to its own concentration gradient
  \details
  \f[
  -\vec{j}_\alpha^{d}(\vec{r}) = 
  +(\kbt\Gamma_\alpha) \nabla C_\alpha(\vec{r})
  \f]
  \param[in] diff_flux_x (negative) single ion diffusion flux (r-space)
  \param[in] conc_k concentration field (k-space)
  \param[in] onsager_coeff Onsager coefficient of species
  \param[in,out] dmy_value working memory to compute gradient of the concentration field in reciprocal space
 */
void Diffusion_flux_single(double **diff_flux_x
			    ,double *conc_k
			    ,double &onsager_coeff
			    ,double *dmy_value // working memory
			   );

/*!
  \brief Compute the total amount of solute (outside the particle domain) for  all the species from the concentration fields
  \param[out] n_solute total solute amount for each species
  \param[in] conc_k concentration field for each species (k-space)
  \param[in] p particle data
  \param[in,out] phi working memory to compute smooth particle field
  \param[in,out] dmy_value working memory to compute real-space concentration field for a single species
  \see Count_solute_each
 */
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
