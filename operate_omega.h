/*!
  \file operate_omega.h
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief High-level routines to handle basic operations on velocity field (header file)
  \see \ref page_design_fsolver section of manual for further details
 */
#ifndef OPERATE_OMEGA_H
#define OPERATE_OMEGA_H

#include "variable.h"
#include "fft_wrapper.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/*!
  \brief Compute the reduced advection term appearing on the rhs of the NS equation (reciprocal space) from the velocity field (real space)
  \details \f[
  \vec{u}(\bm{r})\longrightarrow -\ft{\vec{\Omega}}^*(\vec{k})
  \f]
  The \f$\vec{u}\vec{u}\f$ term is computed in real space, the result is Fourier transformed, and the derivatives are computed by multiplying by the appropriate wavevector. Only the two linearly independent components of the final result are kept.
  \param[in,out] u velocity field (input: real space), off-diagonal components of dyadic product (output: reciprocal space) on output
  \param[out] advection (negative) reduced advection term (reciprocal space)
 */
void U2advection_k(double **u, double **advection);

/*!
  \brief Compute the reduced advection term appearing on the rhs of the NS equation from the reduced vorticity field (reciprocal space)
  \details \f[
  \ft{\vec{\zeta}} = \ft{\vec{\omega}}^* \longrightarrow -\ft{\vec{\Omega}}^*
  \f]
  The vorticity field is converted back to real space and U2advection_k is called
  \param[in] zeta reduced vorticity field (reciprocal space)
  \param[in] uk_dc zero-wavenumber Fourier transform of the velocity field
  \param[out] advection (negative) reduced advection term (reciprocal space)
 */
void Zeta_k2advection_k(double **zeta, double uk_dc[DIM], double **advection);

void Zeta_k2advection_k_OBL(double **zeta, double uk_dc[DIM], double **advection);

/*!
  \brief Compute viscous term appearing in the rhs of the Navier-Stokes equation for the reduced vorticity
  \details \f[
  \ft{\vec{f}} \longrightarrow \ft{\vec{f}} - \nu (2\pi k)^2 \ft{\vec{\zeta}}
  \f]
  \param[in,out] f reduced field variable to update with viscous term
  \param[in] zeta current reduced vorticity field
  \param[in] ijk_range range indices for field update
 */
void Add_zeta_viscous_term(double **zeta, double **f, const Index_range &ijk_range);

/*!
  \brief Enforce zero divergence of field u (in reciprocal space)
  \details \f[
  \ft{\vec{u}} \longrightarrow \ft{\vec{u}} - 
  \frac{(\ft{\vec{u}}\cdot{\vec{k}})}{\norm{k}}\vec{k}
  \f]
  \param[in,out] u Fourier transform of vector field
 */
void Solenoidal_uk(double **u);

void Solenoidal_uk_OBL(double **u);

/*!
  \brief Compute Fourier transform of vector field u
  \details \f[\vec{u}(\vec{r}) \longrightarrow \ft{\vec{u}}(\vec{k})\f]
  \param[in,out] u vector field to transform
 */
inline void U2u_k(double **u){
    A2a_k(u[0]);
    A2a_k(u[1]);
    A2a_k(u[2]);
}

/*!
  \brief Compute inverse Fourier transform of vector field u
  \details \f[\ft{\vec{u}}(\vec{k}) \longrightarrow \vec{u}(\vec{r})\f]
  \param[in,out] u Fourier transform of vector field to inverse transform
 */
inline void U_k2u(double **u){
    A_k2a(u[0]);
    A_k2a(u[1]);
    A_k2a(u[2]);
}

/*!
  \brief Enforce zero divergence of field u (in real space)
  \details \f[
  \vec{u} \longrightarrow \vec{u}', \qquad \nabla\cdot\vec{u}' = 0
  \f]
  \param[in,out] u vector field (real space) to transform
 */
inline void Solenoidal_u(double **u){
  U2u_k(u);
  for(int d=0; d<DIM; d++){
      Truncate_two_third_rule(u[d]);
  }
  Solenoidal_uk(u);
  U_k2u(u);
}

/*!
  \brief Filter high frequency components of vector field according to Orzag's 2/3 rule
  \param[in,out] vector vector field to de-alias
  \param[in] dim dimension of vectors
 */
inline void Truncate_vector_two_third_rule(double **vector,const int &dim){
  for(int d=0; d<dim; d++){
      Truncate_two_third_rule(vector[d]);
  }
}

/*!
  \brief Compute velocity field (real-space) from reduced vorticity field (reciprocal space)
  \details \f[
  \ft{\vec{\zeta}}(\bm{k}) \longrightarrow \vec{u}(\vec{r})
  \f]
  \param[in] zeta reduced vorticity field (reciprocal space)
  \param[in] uk_dc zero wavenumber Fourier transform of velocity field
  \param[out] u velocity field (real space) corresponding to zeta
 */
inline void Zeta_k2u(double **zeta, double uk_dc[DIM], double **u){
  Zeta_k2u_k(zeta, uk_dc, u);
  U_k2u(u);
}

/*!
  \brief Compute contravariant velocity field (real-space) from
  reduced vorticity field (reciprocal space). Also returns a copy of
  the y component of the Fourier transform of the velocity field
  \details \f[
  \ft{\zeta}^\alpha(\vec{k})\longrightarrow u^\alpha(\vec{r}), \ft{u}^y(\vec{k})
  \f]
  \param[in] zeta reduced contravariant vorticity field (reciprocal
  space)
  \param[in] uk_dc zero wavenumber Fourier transform of the
  contravariant velocity field
  \param[out] u contravariant velocity field (real space)
  \param[out] uk_cp contravariant y-velocity field (reciprocal space)
 */
inline void Zeta_k2u_cpuky(double **zeta, double uk_dc[DIM], double **u, double *uk_cp){
  Zeta_k2u_k_OBL(zeta, uk_dc, u);//contra

  int im;
  for(int i = 0; i < NX; i++){
      for(int j = 0; j < NY; j++){
	  for(int k = 0; k < NZ_; k++){
	      im=(i*NY*NZ_)+(j*NZ_)+k;
	      uk_cp[im] = u[1][im];
	  }
      }
  }
  U_k2u(u);
}

/*!
  \brief Compute contravariant vorticity field (real space) from
  contravariant reduced vorticity field (reciprocal space)
  \details \f[
  \ft{\zeta}^\alpha(\vec{k}) \longrightarrow \ft{\omega}^\alpha(\vec{r})
  \f]
  \param[in] zeta contravariant reduced vorticity field (reciprocal
  space)
  \param[out] omega contravariant vorticity field (real space)
 */
inline void Zeta_k2omega_OBL(double **zeta, double **omega){
    Zeta_k2omega_k_OBL(zeta, omega);
    
    for(int d = 0; d < DIM;d ++){
	A_k2a(omega[d]);
    }
}

/*!
  \brief Compute reduced vorticity field (reciprocal space) from velocity field (real space)
  \details \f[
  \vec{u}(\vec{r})\longrightarrow \ft{\vec{\zeta}}(\vec{k})
  \f]
  \param[out] zeta reduced vorticity field corresponding to u (reciprocal space)
  \param[in,out] u velocity field (input:real space, output: reciprocal space)
  \param[in] uk_dc zero wavenumber Fourier transform of velocity field u
 */
inline void U2zeta_k(double **zeta, double uk_dc[DIM], double **u){
  U2u_k(u);
  U_k2zeta_k(u, zeta, uk_dc);
}
#endif
