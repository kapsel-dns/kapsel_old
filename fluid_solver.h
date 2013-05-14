/*!
  \file fluid_solver.h
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief Solver for Navier-Stokes equations (header file)
 */
#ifndef FLUID_SOLVER_H
#define FLUID_SOLVER_H

#include <math.h>
#include "variable.h"
#include "input.h"
#include "f_particle.h"
#include "operate_omega.h"
#include "solute_rhs.h"
#include "fluct.h"

extern double *Pressure;
extern double **Shear_force;
extern double **Shear_force_k;
extern double **f_ns0;
extern double **f_ns1;


//Documentation for inline functions defined in fluid_solver.cxx

/*!
  \fn void Field_solver_Euler(const int &dim, double **zeta_k, const CTime &jikan, double **rhs, const Index_range &ijk_range)
  \brief solver euler
 */

/*!
  \fn void Rhs_NS(double **zeta, double uk_dc[DIM], double **rhs, const Index_range *ijk_range, const int &n_ijk_range, Particle *p)
  \brief rhs of ns
 */

/*!
  \fn void Rhs_solvent(double **u_solvent, double **rhs_solvent, double rhs_uk_dc[DIM])
  \brief rhs
*/

/*!
  \fn void Rhs_NS_solute(Particle *p, double **zeta_k, double uk_dc[DIM], double **u, double **concentration_k, double **rhs_ns, double **rhs_solute, const Index_range *ijk_range, const int &n_ijk_range, double **solute_flux, double **grad_potential, double **surface_normal, double rhs_uk_dc[DIM])
  \brief rhs + ns
 */

/*!
  \fn Add_constant_field_k(double **grad_potential_k, double e_ext[DIM], const CTime &jikan)
  \brief field
 */

/*!
  \fn Rhs_NS_Nernst_Planck(Particle *p, double **zeta, double uk_dc[DIM], double **u, double **concentration_k, double **rhs_ns, double **rhs_solute, const Index_range *ijk_range, const int &n_ijk_range, double **solute_flux, double **grad_potential, double **surface_normal, double rhs_uk_dc[DIM], const CTime &jikan)
  \brief ns + np
 */

// End of inline documentation

/*!
  \brief Allocate workspace variables 
 */
void Mem_alloc_NS_solver(void);

/*!
  \brief Solve Navier-Stokes equation to update reduced vorticity field
  \details \f[
  \ft{\vec{\zeta}} \longrightarrow \ft{\vec{\zeta}} + \left(e^{-\nu (2\pi k)^2 h} - 1\right)\left[\ft{\vec{\zeta}} + \frac{\ft{\vec{\Omega}}}{\nu(2\pi k)^2}\right]
  \f]
  Analytic solution is valid in the absence of solute terms, with no shear.
  \param[in,out] zeta reduced vorticity field
  \param[in] jikan time data
  \param[in] uk_dc zero-wavenumber Fourier transform of the velocity field
  \param[in] ijk_range field iterator parameters for update
  \param[in] n_ijk_range field iterator parameters for update
  \param[in] p particle data (unused)
  \see \ref page_design_fsolver section of manual for further details.
  \todo Specify the difference between ijk_range and n_ijk_range
 */
void NS_solver_slavedEuler(double **zeta, const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p);

// Shear_Navier_Stokes
void NS_solver_slavedEuler_Shear_PBC(double **zeta, const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p, double **force);

void NS_solver_slavedEuler_Shear_OBL(double **zeta, const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p, double **force);

//Type of Rheology 
void Mean_shear_sustaining_yforce_PBC(double **u, double **force,const CTime &jikan);

void Mean_shear_sustaining_kforce_PBC(double **u, double **force,const CTime &jikan);

//void Mean_shear_sustaining_kforce_PBC(Value u[DIM], Value force[DIM],const CTime &jikan);
void Mean_shear_sustaining_force_PBC_OBL(double **u);

void Ion_diffusion_solver_Euler(double **zeta
				,const CTime &jikan
				,double uk_dc[DIM]
				,double **concentration_k
				,Particle *p
				,const Index_range *ijk_range
				,const int &n_ijk_range
				);
/*!
  \brief Solves Navier-Stokes and Poisson-Nernst Plank equation to update the fluid velocity and the electrolyte concentration
  \details Calls Rhs_NS_Nernst_Plank to compute the various terms appearing in the right hand side of the relevant equations, then proceeds to call Field_solver_Euler to update the reduced vorticity field \f$\vec{\zeta}\f$ and the solute concentration \f$C_{\alpha}\f$.
  \param[in,out] zeta reduced vorticity field (reciprocal space)
  \param[in] jikan time data
  \param[in,out] uk_dc zero-wavenumber Fourier transform of velocity field
  \param[in,out] concentration_k solute concentration field (reciprocal space)
  \param[in] p particle data
  \param[in] ijk_range field iterator data
  \param[in] n_ijk_range field iterator data
  \note The solute source term appearing in the rhs of the Navier-Stokes equation has not been included. This section of the code should probably be removed for clarity...
  \see \ref page_design_ssolver section of manual for furter details.
 */
void NSsolute_solver_Euler(double **zeta
			   ,const CTime &jikan
			   ,double uk_dc[DIM]
			   ,double **concentration_k
			   ,Particle *p
			   ,const Index_range *ijk_range
			   ,const int &n_ijk_range
			   );

inline void Calc_shear_stress(const CTime &jikan
			      ,Particle *p
			      ,double *phi
			      //,double **force
			      ,double **Shear_force
			      ,double stress[DIM][DIM]
			      ){

    
    {
	Reset_phi(phi);
	Make_rho_field(phi, p);
	//Make_phi(phi, p);
    }

   int im;
   if(Shear_AC){
   static double factor=Ivolume/jikan.dt_fluid;
   int jm;
   Inertia_stress=0.;
#pragma omp parallel for schedule(dynamic, 1) reduction(+:Inertia_stress) private(im,jm)
   for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
       for(int k=0; k<NZ; k++){
        jm=j;   
		im=(i*NY*NZ_)+(j*NZ_)+k;
        Inertia_stress +=  (double)jm*DX*(phi[im]*ucp[0][im]-rhop[im]); 
        rhop[im] = phi[im]*ucp[0][im]; 
       }
     }
   }
   Inertia_stress *= factor;
  }
    
    double mean_force_x = 0.;
#pragma omp parallel for schedule(dynamic, 1) reduction(+:mean_force_x) private(im)
	  for(int i=0; i<NX; i++){
		for(int j=0; j<NY; j++){
	      for(int k=0; k<NZ; k++){
			im=(i*NY*NZ_)+(j*NZ_)+k;
			//force[0][im] *= phi[im];
			Shear_force[0][im] *= phi[im];
			//mean_force_x += force[0][im];
			mean_force_x += Shear_force[0][im];
		   }
	     }
	   }
	static const double ivolume = Ivolume * POW3(DX);
	mean_force_x *= ivolume;

	 double stress_yx = 0.0;

	double y;
	double fx;
#pragma omp parallel for schedule(dynamic,1) reduction(+:stress_yx) private(y,fx,im)
	    for(int i=0; i<NX; i++){
		for(int j=0; j<NY; j++){
		for(int k=0; k<NZ; k++){
			im=(i*NY*NZ_)+(j*NZ_)+k;
			y = (double)j*DX;
			//fx = force[0][im] - mean_force_x;
			fx = Shear_force[0][im] - mean_force_x;
		    stress_yx += y * fx;
		}
		}
	    }
	 const double dmy = ivolume/jikan.dt_fluid; //Select this 
	 stress[1][0]=stress_yx*dmy;
}

inline void Calc_hydro_stress(const CTime &jikan
			      ,const Particle *p
			      ,double *phi
			      ,double *force
			      ,double stress[DIM][DIM]
    ){
    
    double stress_yx = 0.0;
    
    static const double ivolume = Ivolume * POW3(DX);
#pragma omp parallel for reduction(+:stress_yx) 
    for(int i = 0; i < NX; i++){
	for(int j = 0; j < NY; j++){
	    for(int k = 0; k < NZ; k++){
		/*
		int jj = (int)fmod((double)targetJ + j + NY, NY); 

		double y = (double)j*DX;
		double fx = force[i][jj][k] - mean_force_x;
		stress_yx += y * fx;
		*/
		stress_yx += force[i*NY*NZ_ + j*NZ_ + k];
	    }
	}	
    }
    
    static const double dmy = (double)ivolume/jikan.dt_fluid; //Select this 
    stress[1][0]=stress_yx*dmy;
    
}
#endif


