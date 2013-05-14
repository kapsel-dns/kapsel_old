/*!
  \file operate_electrolyte.h
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief Routines to compute the charge distributions and forces (header file)
 */
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

/*!
  \brief Determine how to initialize the solute density
  \param[out] Concentration inital concentration fields specified on input for all solute species
  \param[in] p particle data
  \param[in] jikan time data
 */
void Init_rho_ion(double **Concentration, Particle *p, CTime &jikan);

/*!
  \brief Allocate the requires memory for charged species data
 */
void Mem_alloc_charge(void);

void Calc_free_energy_PB(double **conc_k
			  ,Particle *p
			  ,double *free_energy
			  ,double *phi // working memory
			  ,double *charge_density // working memory
			  ,double *dmy_value // working memory
			  ,const CTime &jikan
			  );

/*!
  \brief Compute current smooth particle and charge field profiles (real space)
  \details Given the current particle data (positions), computes the smooth particle profile \f$\phi\f$ and the colloid surface charge density \f$\sigma_e\f$
  \param[out] phi smooth particle profile field
  \param[out] surface colloid surface charge density
  \param[in] p Current particle data
 */
void Make_phi_qq_particle(double *phi
			   ,double *surface
			   ,Particle *p);

void Make_phi_qq_fixed_particle(double *phi
			   ,double *surface
			   ,Particle *p);


/*!
  Compute the solute source term appearing in the rhs of the Navier-Stokes equation (including possible contributions from the external field)
  \details
  \f[
  -\sum_\alpha Z_\alpha e C_\alpha(\vec{r}) \nabla\Phi^{\text{tot}}(\vec{r})
  \f]
  \param[in] force Coulomb force exerterd on the fluid due to the presence of the charged solute
  \param[in] p particle data
  \param[in] conc_k solute concentration field (reciprocal space)
  \param[out] charge_density auxiliary field to compute total charge density
  \param[out] potential auxiiary field to compute total electrostatic potential
  \param[out] jikan time data
 */
void Make_Coulomb_force_x_on_fluid(double **force
				    ,Particle *p
				    ,double **conc_k
				    ,double *charge_density // working memory
				    ,double *potential // working memory
				    ,const CTime &jikan
				    );

/*!
  \brief Computes the total charge density field (real space)
  \details Computes the total charge density field  by adding the colloid surface charge density (for the current particle positions) to the solute charge density, where the total charge density is given by
  \f[
  \rho_e = (1-\phi)\sum_\alpha Z_\alpha e C_\alpha + \norm{\nabla \phi} e \sigma_e
  \f]
  with \f$\sigma_e\f$ the colloid surface charge density (defined over the particle/fluid interface).
  \param[in] p particle data
  \param[in] conc_k solute concentration field (reciprocal space) 
  \param[out] charge_density total charge density (real space)
  \param[in,out] phi auxiliary field to compute the smooth profile
  \param[in,out] dmy_value auxiliary field to compute the concentration of a single solute species
 */
void Conc_k2charge_field(Particle *p
			 ,double **conc_k
			 ,double *charge_density
			 ,double *phi // working memory
			 ,double *dmy_value // working memory
			 );


/*!
  \brief Compute the electrostatic potential due to all charges in the system (reciprocal space)
  \details Given the Fourier Transform of the total charge density \f$\ft{\rho}_e\f$, returns the Fourier Transform of the electrostatic potential by solving Poisson's Equation
  \f[
  \ft{\Phi}(\vec{k}) = \frac{\ft{\rho}_e}{(2\pi k)^2 \epsilon}
  \f]
  \param[in,out] potential charge density on input, electrostatic potential on output
 */
void Charge_field_k2Coulomb_potential_k_PBC(double *potential);

#endif

