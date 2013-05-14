/*!
  \file operate_surface.h
  \author J. Molina
  \date 2012/06/22
  \version 1.0
  \brief Routines to control the slip velocity at particle fluid
  boundaries (header file)
  \details \see \ref page_design_swimmer for a detailed description
 */
#ifndef OPERATE_SURFACE_H
#define OPERATE_SURFACE_H

#include "variable.h"
#include "make_phi.h"
#include "particle_solver.h"

/*!                                             
  \fn void slip_droplet(double *vp, double *wp, double
  *delta_v, double *delta_w, Particle p)
  \brief Computes velocity of fluid droplet needed to ensure momentum
  conservation if surface slip profile is enforced
  \details \f[
  \{\vec{V}_i^\prime, \vec{\Omega}_i^\prime\} \longrightarrow
  \{\delta\vec{V}_i, \delta\vec{\Omega}_i\}
  \f]
  \param[out] vp particle velocity used to enforce surface slip
  \param[out] wp particle angular velocity used to enforce surface
  slip
  \param[out] delta_v velocity of coutner-flow needed to maintain
  momentum conservation
  \param[out] delta_w angular velocity of counter-flow needed to
  momentum conservation
  \param[in] p particle data
 */

/*!
  \brief Computes geometric factors related to particle and surface
  fluid mass distributions needed for surface slip implementation
  \param[in] u current velocity field 
  \param[in] p particle data
 */
void Make_particle_momentum_factor(double const* const* u, Particle *p);

/*!
  \brief Computes the change in the velocity field needed to enforce
  the slip profile at the particle/fluid boundaries
  \details \f[
  u^{**} - u^{*} = \sum_{i=1}^N\varphi_i\left(\vec{V}^\prime_i +
  \vec{\Omega}^\prime_i\times\vec{r}_i + \vec{u}^s_i -
  \vec{u}^*\right) + \sum_{i=1}^N\phi_i\left(\delta\vec{V}_i + \delta\vec{\Omega}_i\times\vec{r}_i\right)
  \f]
  \param[out] up force density field needed to impose surface slip
  profiles
  \param[in] u current velocity field
  \param[in] p particle data
  \param[in] jikan time data
 */
void Make_force_u_slip_particle(double **up, double const* const* u,
				Particle *p, const CTime &jikan);

/*!
  \brief Computes convergence factor for iterative slip velocity as
  maximum of squared relative difference between succesive velocities
  \details \f[
  \epsilon =
  \text{Max}_i\left(
  \frac{\norm{\vec{V}_i-\vec{V}_i^\prime}}{\norm{\vec{V}_i}}
  + \frac{\norm{\vec{\Omega}_i - \vec{\Omega}_i^\prime}}{\norm{\vec{\Omega}_i}}
  \right)
  \f]
  \param[in] p particle data
 */
inline double Slip_particle_convergence(Particle *p){
  double eps,dmy_v, dmy_w, nv, nw;
  eps = 0.0;

  for(int n = 0; n < Particle_Number; n++){
    dmy_v = dmy_w = nv = nw = 0.;
    if(janus_propulsion[p[n].spec] == slip){
      for(int d = 0; d < DIM; d++){
	nv += p[n].v[d] * p[n].v[d];
	nw += p[n].omega[d] * p[n].omega[d];
	
	dmy_v += (p[n].v_slip[d] - p[n].v[d]) * (p[n].v_slip[d] - p[n].v[d]);
	dmy_w += (p[n].omega_slip[d] - p[n].omega[d]) * (p[n].omega_slip[d] - p[n].omega[d]);
      }
      if(positive_mp(nv)){
	dmy_v /= nv;
      }
      if(positive_mp(nw)){
	dmy_w /= nw;
      }
    }
    eps = MAX(eps, dmy_v + dmy_w);
  }
  return sqrt(eps);
}

/*!
  \brief Set the velocity to be used for the surface slip at the next iteration
 */
inline void Update_slip_particle_velocity(Particle *p, const int &iter){
  if(iter == 0){
#pragma omp parallel for schedule(dynamic, 1)
    for(int n = 0; n < Particle_Number; n++){
      for(int d = 0; d < DIM; d++){
	p[n].v_slip[d] =  p[n].v[d];
	p[n].omega_slip[d] = p[n].omega[d];
      }
    }
  }else{
    double dmy_new = 0.7;
    double dmy_old = 1.0 - dmy_new;
#pragma omp parallel for schedule(dynamic, 1)
    for(int n = 0; n < Particle_Number; n++){
      for(int d = 0; d < DIM; d++){
	p[n].v_slip[d] = dmy_old*p[n].v_slip[d] + dmy_new*p[n].v[d];
	p[n].omega_slip[d] = dmy_old*p[n].omega_slip[d] + dmy_new*p[n].omega[d];
      }
    }
  }
}
inline void update_angular(double* l, const double *x, const double *v, const double alpha=1.0){
  l[0] += alpha * (x[1] * v[2] + x[2] * v[1]);
  l[1] += alpha * (x[2] * v[0] + x[0] * v[2]);
  l[2] += alpha * (x[0] * v[1] + x[1] * v[0]);
}
inline void momentum_check(double const* const* up, Particle *p, const CTime &jikan){
  ////////////////////////
  const double dx = DX;
  const int np_domain = NP_domain;
  int const* const* sekibun_cell = Sekibun_cell;
  const int* Nlattice = Ns;
  const double radius = RADIUS;
  static const double dmy0 = DX3 * RHO;
  /////////////////////////
  double xp[DIM], vp[DIM], omega_p[DIM], v_rot[DIM], r[DIM], x[DIM], fv[DIM],residue[DIM], fu[DIM];
  int x_int[DIM], r_mesh[DIM];
  int sw_in_cell, pspec;
  double dmy_r, dmy_phi, dmy_phic, dmy_mass, dmy_xi;
  /////////////////////////
  double tot_p[DIM], part_p[DIM], fluid_p[DIM], md_p[DIM], full_p[DIM];
  double tot_int_p[DIM], part_int_p[DIM], fluid_int_p[DIM], md_int_p[DIM];

  double tot_w[DIM], part_w[DIM], fluid_w[DIM], md_w[DIM], full_w[DIM];
  double tot_int_w[DIM], part_int_w[DIM], fluid_int_w[DIM], md_int_w[DIM];
  for(int d = 0; d < DIM; d++){
    tot_p[d] = part_p[d] = fluid_p[d] = md_p[d] = 0.0;
    tot_int_p[d] = part_int_p[d] = fluid_int_p[d] = md_int_p[d] = 0.0;

    tot_w[d] = part_w[d] = fluid_w[d] = md_w[d] = 0.0;
    tot_int_w[d] = part_int_w[d] = fluid_int_w[d] = md_int_w[d] = 0.0;
  }
  dmy_mass = 0.0;
  for(int n = 0; n < Particle_Number; n++){
    pspec = p[n].spec;
    for(int d = 0; d < DIM; d++){
      xp[d] = p[n].x[d];
      vp[d] = p[n].v[d];
      omega_p[d] = p[n].omega[d];
    }

    sw_in_cell = Particle_cell(xp, dx, x_int, residue);
    sw_in_cell = 1;
    for(int mesh = 0; mesh < np_domain; mesh++){
      Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
      Angular2v(omega_p, r, v_rot);
      int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
      for(int d = 0; d < DIM; d++){
	x[d] = r_mesh[d] * dx;
	fu[d] = up[d][im];
      }
      dmy_r = Distance(x, xp);
      dmy_phi = Phi(dmy_r, radius);
      dmy_phic = 1.0 - dmy_phi;
      dmy_mass += (1.0 - dmy_phi);
      dmy_xi = ABS(dmy_r - radius);
      for(int d = 0; d < DIM; d++){
	fv[d] = (vp[d] + v_rot[d]);
      }

      for(int d = 0; d < DIM; d++){
	md_p[d] += dmy_phi * fv[d];
	part_p[d] += dmy_phi * fu[d];
	fluid_p[d] += dmy_phic * fu[d];
      }
      {
	update_angular(md_w, x, fv, dmy_phi);
	update_angular(part_w, x, fu, dmy_phi);
	update_angular(fluid_w, x, fu, dmy_phic);
      }

      if(dmy_r <= radius + dx){
	for(int d = 0; d < DIM; d++){
	  tot_p[d] += fu[d];
	}
	{
	  update_angular(tot_w, x, fu);
	}
	
	if(dmy_xi <= HXI){
	  for(int d = 0; d < DIM; d++){
	    tot_int_p[d] += fu[d];
	    md_int_p[d] += dmy_phi * fv[d];
	    part_int_p[d] += dmy_phi * fu[d];
	    fluid_int_p[d] += dmy_phic * fu[d];
	  }
	  {
	    update_angular(tot_int_w, x, fu);
	    update_angular(md_int_w, x, fv, dmy_phi);
	    update_angular(part_int_w, x, fu, dmy_phi);
	    update_angular(fluid_int_w, x, fu, dmy_phic);
	  }
	}
      }
    }
  }

  for(int d = 0; d < DIM; d++){
    full_p[d] = fluid_p[d] = 0.0;
    full_w[d] = fluid_w[d] = 0.0;
  }
  int outside = 1;
  for(int i = 0; i < NX; i++){
    x[0] = (double)i * dx;
    for(int j = 0; j < NY; j++){
      x[1] = (double)j * dx;
      for(int k = 0; k < NZ; k++){
	x[2] = (double)k * dx;

	int im = (i * NY * NZ_) + (j * NZ_) + k;
	for(int d = 0; d < DIM; d++){
	  fu[d] = up[d][im];
	  full_p[d] += fu[d];
	}
	update_angular(full_w, x, fu);

	for(int n = 0; n < Particle_Number; n++){
	  pspec = p[n].spec;
	  for(int d = 0; d < DIM; d++){
	    xp[d] = p[n].x[d];
	    vp[d] = p[n].v[d];
	    omega_p[d] = p[n].omega[d];
	  }
	  dmy_r = Distance(x, xp);
	  dmy_phi = Phi(dmy_r, radius);
	  dmy_phic = 1.0 - dmy_phi;
	  if(dmy_phi > 0.){
	    outside = 0.;
	    break;
	  }
	}

	if(outside){
	  for(int d = 0; d < DIM; d++){
	    fluid_p[d] += up[d][im];
	  }
	  update_angular(fluid_w, x, fu);
	}
      }
    }
  }
  fprintf(stderr, "%10.8E %10.8E %10.8E %10.8E %10.8E %10.8E %10.8E %10.8E %10.8E %10.8E %10.8E %10.8E\n",
	  tot_p[0], tot_p[1], tot_p[2], tot_w[0], tot_w[1], tot_w[2],
	  full_p[0], full_p[1], full_p[2], full_w[0], full_w[1], full_w[2]);
}

#endif
