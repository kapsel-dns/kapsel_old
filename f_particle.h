/*!
  \file f_particle.h
  \author Y. Nakayama
  \date 2006/08/15
  \version 1.2
  \brief Compute Hydrodynamic force on particle (header file)
 */
#ifndef F_PARTICLE_H
#define F_PARTICLE_H

#include "variable.h"
#include "operate_omega.h"
#include "input.h"
#include "make_phi.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern double **f_particle;

/*!
  \brief Allocate f_particle working array
 */
void Mem_alloc_f_particle(void);

/*!
  \brief Compute the force density field which imposes the rigidity constraint on the total velocity field
  \details The force density field \f$\phi\vec{f}_p\f$ is derived assuming momentum conservation between colloids and fluid at each time step
  \f[
  \int_{t_n}^{t_n+h} \df{s}\, \phi\vec{f}_p = 
  \phi^{n+1}\left(\vec{v}_p^{n+1} - \vec{u}\right) - \frac{h}{\rho}\nabla p
  \f]
  where \f$\phi^{n+1}\vec{v}_p^{n+1}\f$ represents the particle velocity field and \f$\vec{u}\f$ the intermediate fluid velocity field (obtained from solving the NS equations for the total fluid). 
  \note It is not necessary to include the pressure term, one just needs to enforce the solenoidal condition on the final velocity field
  \param[out] f force density field needed to enforce rigidity constraint on the final field
  \param[in] u total fluid velocity field
  \param[in] up particle velocity field
  \param[in] phi particle concentration field
 */
inline void Make_f_particle_dt_nonsole(double **f
				       ,double **u 
				       ,double **up
				       ,double *phi
				       ){
  // !! これを呼ぶ前に、 up, phi を計算しておくこと。
  //
  // f = up - phi *u
  //
  int im;
  {
#pragma omp parallel for schedule(dynamic, 1) private(im)
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ; k++){
	  im=(i*NY*NZ_)+(j*NZ_)+k;
	  f[0][im] = up[0][im] - phi[im] * u[0][im];
	  f[1][im] = up[1][im] - phi[im] * u[1][im];
	  f[2][im] = up[2][im] - phi[im] * u[2][im];
	}
      }
    }	
  }
}


/*!
  \brief Compute the (solenoidal) force density field which imposes the rigidity constraint on the total velocity field
  \see Make_f_particle_dt_nonsole
*/				
inline void Make_f_particle_dt_sole(double **f
				    ,double **u
				    ,double **up
				    ,double *phi
				    ){
    // !! これを呼ぶ前に、 up, phi を計算しておくこと。
    Make_f_particle_dt_nonsole(f, u, up, phi);
  {
    // f = up - phi *u
     Solenoidal_u(f); // div f = 0
  }
}
inline void Add_f_particle(double **u, double **f, const int dim=DIM){
    // !! これを呼ぶ前に、 f を計算しておくこと。
int im;
   // for(int d=0; d<dim; d++){
#pragma omp parallel for schedule(dynamic, 1) private(im) 
    	for(int i=0; i<NX; i++){
	    for(int j=0; j<NY; j++){
		for(int k=0; k<NZ; k++){
		 im=(i*NY*NZ_)+(j*NZ_)+k;
		 u[0][im] += f[0][im];
		 u[1][im] += f[1][im];
		 u[2][im] += f[2][im];
		}
	    }
	    }	
  //  }
}

#endif
