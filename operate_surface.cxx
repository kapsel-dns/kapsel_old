#include "operate_surface.h"
/*!
  \file operate_surface.cxx
  \author J. Molina
  \date 2012/06/22
  \version 1.0
  \brief Routines to control the slip velocity at particle fluid
  boundaries
  \details \see \ref page_design_swimmer for a detailed description
 */

void Make_particle_momentum_factor(double const* const* u, Particle *p){
  //////////////////////////////
  const double dx = DX;
  const double dx3 = DX3;
  const int np_domain = NP_domain;
  int const* const* sekibun_cell = Sekibun_cell;
  const int* nlattice = Ns;
  const double radius = RADIUS;
  //////////////////////////////
  int sw_in_cell, pspec;
  int x_int[DIM], r_mesh[DIM];
  double dmy_r, dmy_sr, dmy_phi, dmy_phi_s;
  double xp[DIM], r[DIM], x[DIM], residue[DIM], u_fluid[DIM];

  double dmy_xi, dmy_theta, dmy_tau;
  double n_r[DIM], n_theta[DIM], n_tau[DIM];
  double slip_mode, slip_vel, slip_magnitude;
  double r_x_us[DIM], us[DIM];

  double M0, SM0; // mass
  double M1[DIM], SM1[DIM]; // center of mass
  double M2[DIM][DIM], SM2[DIM][DIM]; // moment of inertia
  double dv_s[DIM], dw_s[DIM]; //momentum change due to slip at surface

#pragma omp parallel for schedule(dynamic, 1) \
  private(sw_in_cell, pspec, x_int, r_mesh, dmy_r, dmy_sr, dmy_phi, dmy_phi_s, xp, \
	  r, x, residue, u_fluid, dmy_xi, dmy_theta, dmy_tau, \
	  n_r, n_theta, n_tau, slip_mode, slip_vel, slip_magnitude, r_x_us, us, \
          M0, SM0, M1, SM1, M2, SM2, dv_s, dw_s)
  for(int n = 0; n < Particle_Number; n++){
    pspec = p[n].spec;

    for(int d = 0; d < DIM; d++){
      xp[d] = p[n].x[d];
     
      M1[d] = SM1[d] = dv_s[d] = dw_s[d] = 0.0;
      for(int l = 0; l < DIM; l++){
	M2[d][l] = SM2[d][l] = 0.0;
      }
    }
    M0 = SM0 = 0.0;
    
    sw_in_cell = Particle_cell(xp, dx, x_int, residue);
    sw_in_cell = 1;
    for(int mesh = 0; mesh < np_domain; mesh++){
      Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, nlattice, dx, r_mesh, r);
      int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
      for(int d = 0; d < DIM; d++){
	x[d] = r_mesh[d] * dx;
	u_fluid[d] = u[d][im];
      }
      dmy_r = Distance(x, xp);
      dmy_phi = Phi(dmy_r, radius);
      dmy_xi = ABS(dmy_r - radius);
      
      //particle properties
      M0 += dmy_phi;
      for(int d = 0; d < DIM; d++){
	M1[d] += dmy_phi * r[d];
	M2[d][d] += dmy_phi * dmy_r * dmy_r;
	for(int l = 0; l < DIM; l++){
	  M2[d][l] -= dmy_phi * r[d] * r[l];
	}
      }

      //fluid surface properties
      if(janus_propulsion[pspec] == slip &&  less_than_mp(dmy_xi, HXI)){

	slip_vel = janus_slip_vel[pspec];  //B1
	slip_mode = janus_slip_mode[pspec];//alpha/2
        dmy_phi_s = (1.0 - dmy_phi) * DPhi_compact_sin_norm(dmy_r, radius);
	
	SM0 += dmy_phi_s;
	Squirmer_coord(r, n_r, n_theta, n_tau, dmy_sr, dmy_theta, dmy_tau, p[n]);
	slip_magnitude = slip_vel * (sin(dmy_theta) + slip_mode * sin(2.0*dmy_theta));
        for(int d = 0; d < DIM; d++){
          us[d] = n_theta[d] * slip_magnitude - u_fluid[d];
        }

	v_cross(r_x_us, r, us);
	for(int d = 0; d < DIM; d++){
	  dv_s[d] += dmy_phi_s * us[d];
	  dw_s[d] += dmy_phi_s * r_x_us[d];
	  
	  SM1[d] += dmy_phi_s * r[d];
	  SM2[d][d] += dmy_phi_s * dmy_r * dmy_r;
	  for(int l = 0; l < DIM; l++){
	    SM2[d][l] -= dmy_phi_s * r[d] * r[l];
	  }
	}
      }//slip surface

    }//mesh

    p[n].mass = M0;
    p[n].surface_mass = SM0;
    for(int d = 0; d < DIM; d++){
      
      p[n].mass_center[d] = M1[d];
      p[n].surface_mass_center[d] = SM1[d];
      p[n].surface_dv[d] = dv_s[d];
      p[n].surface_dw[d] = dw_s[d];
      
      for(int l = 0; l < DIM; l++){
	p[n].inertia[d][l] = M2[d][l];
	p[n].surface_inertia[d][l] = SM2[d][l];
      }
    }

  }//Particle_Number
}

inline void slip_droplet(double *vp, double *wp, double *delta_v, double *delta_w, Particle p){
  double LL[DIM][DIM];
  double dP[DIM], dL[DIM], L0[DIM], mc[DIM];
  double imass, mc2;

  for(int d = 0; d < DIM; d++){
    vp[d] = p.v_slip[d];        //velocity used to enforce surface slip
    wp[d] = p.omega_slip[d];    
    mc[d] = p.mass_center[d];   
    LL[d][0] = LL[d][1] = LL[d][2] = 0.0;
  }
  imass = 1.0 / p.mass;

  //momentum change due to slip at boundary
  double w_x_r[DIM];
  double r_x_v[DIM];
  v_cross(w_x_r, wp, p.surface_mass_center);
  v_cross(r_x_v, p.surface_mass_center, vp);
  for(int d = 0; d < DIM; d++){
    dP[d] = -(p.surface_mass * vp[d] + w_x_r[d] + p.surface_dv[d]);
    dL[d] = -(r_x_v[d] + 
              p.surface_inertia[d][0] * wp[0] + p.surface_inertia[d][1] * wp[1] + p.surface_inertia[d][2] * wp[2] +
              p.surface_dw[d]);
  }

  // Compute angular velocity of droplet
  v_cross(L0, mc, dP, imass);
  mc2 = -(mc[0] * mc[0] + mc[1] * mc[1] + mc[2] * mc[2]) * imass;
  for(int d = 0; d < DIM; d++){
    L0[d] = dL[d] - L0[d];
    LL[d][d] = mc2;

    for(int l = 0; l < DIM; l++){
      LL[d][l] += (p.inertia[d][l] + imass * mc[d] * mc[l]);
    }
  }
  M_inv(LL);
  M_v_prod(delta_w, LL, L0);

  // Compute linear velocity of droplet
  v_cross(L0, delta_w, mc);
  for(int d = 0; d < DIM; d++){
    delta_v[d] = imass * (dP[d] - L0[d]);
  }
}

void Make_force_u_slip_particle(double **up, double const* const* u, Particle *p, const CTime &jikan){
  //////////////////////// System parameters
  const double dx = DX;
  const double dx3 = DX3;
  const int np_domain = NP_domain;
  int const* const* sekibun_cell = Sekibun_cell;
  const int* nlattice = Ns;
  const double radius = RADIUS;

  ////////////////////////  Function variables
  int sw_in_cell, pspec;
  int x_int[DIM], r_mesh[DIM];
  double dmy_r, dmy_sr, dmy_phi, dmy_phi_s;
  double xp[DIM], vp[DIM], omega_p[DIM], v_rot[DIM];
  double delta_v[DIM], delta_w[DIM], delta_v_rot[DIM];
  double r[DIM], x[DIM], residue[DIM], u_fluid[DIM];
  
  double dmy_xi, dmy_theta, dmy_tau;
  double dmy_vslip, slip_mode, slip_vel, dmy_us;
  double n_r[DIM], n_theta[DIM], n_tau[DIM];
  double dmy_fv[DIM], force_s[DIM], torque_s[DIM], force_p[DIM], torque_p[DIM];

#pragma omp parallel for schedule(dynamic, 1) \
  private(sw_in_cell, pspec, x_int, r_mesh, dmy_r, dmy_sr, dmy_phi, dmy_phi_s, xp, vp, omega_p, v_rot, \
	  delta_v, delta_w, delta_v_rot, r, x, residue, u_fluid, \
	  dmy_xi, dmy_theta, dmy_tau, dmy_vslip, slip_mode, slip_vel, dmy_us,  \
	  n_r, n_theta, n_tau, dmy_fv, force_s, torque_s, force_p, torque_p)
  for(int n = 0; n < Particle_Number; n++){
    pspec = p[n].spec;

    for(int d = 0; d < DIM; d++){ //reset slip force/torque
      p[n].f_slip[d] = 0.0;
      p[n].torque_slip[d] = 0.0;
    }
    
    if(janus_propulsion[pspec] == slip){
      
      slip_vel = janus_slip_vel[pspec];
      slip_mode = janus_slip_mode[pspec];
      slip_droplet(vp, omega_p, delta_v, delta_w, p[n]);
      for(int d = 0; d < DIM; d++){
	xp[d] = p[n].x[d];
	force_s[d] = torque_s[d] = 0.0;
	force_p[d] = torque_p[d] = 0.0;
      }
      sw_in_cell = Particle_cell(xp, dx, x_int, residue);
      sw_in_cell = 1;

      for(int mesh = 0; mesh < np_domain; mesh++){
	Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, nlattice, dx, r_mesh, r);
	int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
	for(int d = 0; d < DIM; d++){
	  u_fluid[d] = u[d][im];
	  x[d] = r_mesh[d] * dx;
	}
	dmy_r = Distance(x, xp);
	dmy_phi = Phi(dmy_r, radius);	 
	dmy_xi = ABS(dmy_r - radius);

	{//fluid particle domain (droplet with counter flow) :
         //delta_v, delta_w
          Angular2v(delta_w, r, delta_v_rot);
          for(int d = 0; d < DIM; d++){
            dmy_fv[d] = dmy_phi * (delta_v[d] + delta_v_rot[d]);
            force_p[d] += dmy_fv[d];
            
#pragma omp atomic
            up[d][im] += dmy_fv[d];
          }
          {
            torque_p[0] += (r[1] * dmy_fv[2] - r[2] * dmy_fv[1]);
            torque_p[1] += (r[2] * dmy_fv[0] - r[0] * dmy_fv[2]);
            torque_p[2] += (r[0] * dmy_fv[1] - r[1] * dmy_fv[0]);
          }
        }
	
        dmy_phi_s = (1.0 - dmy_phi) * DPhi_compact_sin_norm(dmy_r, radius);
	if(less_than_mp(dmy_xi, HXI)){//interface domain
          //slip enforced wrt vp, omega_p
	  Angular2v(omega_p, r, v_rot);
	  Squirmer_coord(r, n_r, n_theta, n_tau, dmy_sr, dmy_theta, dmy_tau, p[n]);
	  dmy_vslip = slip_vel * (sin(dmy_theta) + slip_mode * sin(2.0*dmy_theta));

          for(int d = 0; d < DIM; d++){
            dmy_fv[d] = dmy_phi_s * (vp[d] + v_rot[d] + n_theta[d]*dmy_vslip - u_fluid[d]);
          }

	  for(int d = 0; d < DIM; d++){
	    force_s[d] += dmy_fv[d];
#pragma omp atomic
	    up[d][im] += dmy_fv[d];
	  }
	  {
	    torque_s[0] += (r[1] * dmy_fv[2] - r[2] * dmy_fv[1]);
	    torque_s[1] += (r[2] * dmy_fv[0] - r[0] * dmy_fv[2]);
	    torque_s[2] += (r[0] * dmy_fv[1] - r[1] * dmy_fv[0]);
	  }
	}//interface_domain
      }//mesh


      for(int d = 0; d < DIM; d++){
	force_p[d] = -force_p[d];
	torque_p[d] = -torque_p[d];
      }
      if((v_rms(force_p, force_s) > LARGE_TOL_MP || 
          v_rms(torque_p, torque_s) > LARGE_TOL_MP)){
        fprintf(stderr, "###############################");
	fprintf(stderr, "# Momentum Conservation Warning : %10.8E %10.8E\n", v_rms(force_p, force_s), v_rms(torque_p, torque_s));
        fprintf(stderr, "# Force  : %10.8E %10.8E %10.8E %10.8E %10.8E %10.8E\n", force_p[0], force_p[1], force_p[2], force_s[0], force_s[1], force_s[2]);
        fprintf(stderr, "# Torque : %10.8E %10.8E %10.8E %10.8E %10.8E %10.8E\n", torque_p[0], torque_p[1], torque_p[2], torque_s[0], torque_s[1], torque_s[2]);
        fprintf(stderr, "###############################");
	
      }
    }// slip_particle ?
  }// Particle_Number
}

