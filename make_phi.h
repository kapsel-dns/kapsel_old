/*!
  \file make_phi.h
  \brief Routines to compute smooth profile and grid particle properties (header file)
  \author Y. Nakayama
  \date 2006/07/28
  \version 1.3
  \see The details on the \ref sec_design_grid_sekibun in the manual
 */
#ifndef MAKE_PHI_H
#define MAKE_PHI_H

#include <assert.h> 
#include "input.h"
#include "variable.h"
#include "profile.h"
#include "rigid_body.h"
#include "operate_omega.h"

extern void (*Angular2v)(const double *omega, const double *r, double *v);
extern int NP_domain;
extern int **Sekibun_cell;

/*!
  \brief Compute smooth particle position and advection fields
  \details Advection field maps only the linear velocity of the particles, it ignores the angular velocity (in case \c ROTATION is turned on)
  \param[out] phi smooth particle field
  \param[out] up smooth particle advection (linear velocity) field
  \param[in] p particle data
  \see Make_phi_u_particle
 */
void Make_phi_u_advection(double *phi, double **up, Particle *p);

void Make_phi_janus_particle(double *phi,
			     double *id_phi,
			     Particle *p);
void Make_phi_janus_particle_OBL(double *phi,
				 double *id_phi,
				 Particle *p);

/*!
  \brief Compute smooth particle position field
  \details 
  \f[
  \phi_p(\vec{r}) = \sum_i^{N_p} \phi_i(\vec{r})
  \f]
  where \f$\phi_i(\vec{r})\f$ is the profile field of the \f$i\f$-th particle
  \param[out] phi smooth particle position field
  \param[in] p particel data
  \param[in] radius particle radius (domain over which profile field is non-zero)
 */
void Make_phi_particle(double *phi
		       ,Particle *p
		       ,const double radius = RADIUS
		       );

/*!
  \brief Compute smooth particle position and velocity fields
  \details
  \f{eqnarray*}{
  \phi_p(\vec{r}) &=& \sum_i^{N_p} \phi_i(\vec{r}) \\
  \phi\vec{v}_p(\vec{r}) &=& \sum_i^{N_p} \left[
  \vec{v}_i + \vec{\omega}_i\times\left(\vec{r} - \vec{r}_i\right)
  \right]\phi_i(\vec{r})
  \f}
  with \f$\vec{r}_i, \vec{v}_i, \vec{\omega}_i\f$ the position, velocity, and angular velocity of the \f$i\f$-th particle.
  \param[out] phi smooth particle position field
  \param[out] up smooth particle velocity field
  \param[in] p particle data
 */
void Make_phi_u_particle(double *phi, double **up, Particle *p);
void Make_phi_particle_OBL(double *phi
			   ,Particle *p
			   ,const double radius = RADIUS
    );
void Make_phi_u_particle_OBL(double *phi, double **up, Particle *p);

/*!
  \brief Compute surface normal field defined on the interface domain of the particles
  \param[out] surface_normal surface normal field
  \param[in] p particle data
  \warning In case of overlapping interfaces, the surface normal vectors
  are not actually unitary. 
 */
void Make_surface_normal(double **surface_normal
			   ,const Particle *p);


void Make_rho_field(double *phi,Particle *p);

/*!
  \brief Compute the position of the particle on the discrete grid
  \details Returns the coordinates of grid particle, defined as the coordinates of the back-lower-left edge of the bin which contains the particle center
  \f[
  \vec{x}_i = (\text{int}) \vec{x}_p / \df{x}
  \f]
  This is the position of the Sekibun_cell used to compute integrals over the particle domain
  \param[in] xp real space position of the particle
  \param[in] dx grid spacing
  \param[out] x_int position of the particle on the discrete grid
  \param[out] residue displacement vector from the real position to the grid position
 */
inline int Particle_cell(const double *xp
			 ,const double &dx
			 ,int *x_int
			 ,double *residue){
  const double idx = 1./dx;
  {
    assert(xp[0] >= 0);
    assert(xp[0] < L_particle[0]);

    assert(xp[1] >= 0);
    assert(xp[1] < L_particle[1]);

    assert(xp[2] >= 0);
    assert(xp[2] < L_particle[2]);
  }
  
  // r_i の左下のメッシュ座標を計算
  int sw_in_cell=0;
  for(int d=0; d<DIM; d++){
    double dmy = (xp[d] *idx);
    x_int[d] = (int)dmy;

    residue[d] = (dmy - (double)x_int[d]) * dx;
    if( dmy - (double)x_int[d] > 0.0 ){
      sw_in_cell = 1;
    }
  }
  return sw_in_cell;
}

/*!
  \brief Returns global (world) coordinates of a given point of the particle Sekibun integration cell
  \param[in] cell local coordinates of a point in the Sekibun grid
  \param[in] x_int center of the Sekibun grid (origin of local frame)
  \param[in] residue offset between particle position and position of Sekibun cell
  \param[in] sw_in_cell unused
  \param[in] Nlattice (global) number of lattice points in each direction
  \param[in] dx grid spacing
  \param[out] r_mesh global coordinates of \c cell point
  \param[out] r displacement vector from particle center to \c cell point (i.e. radial vector from the particle to the given mesh point)
 */
inline void Relative_coord(const int *cell
			   ,const int *x_int
			   ,const double *residue
			   ,const int &sw_in_cell
			   ,const int Nlattice[DIM]
			   ,const double dx
			   ,int *r_mesh
			   ,double *r){

    for(int d=0;d<DIM;d++){
	r_mesh[d] = (x_int[d] + cell[d] + Nlattice[d] ) % Nlattice[d]; 
	{
	    assert( r_mesh[d] < Nlattice[d] );
	    assert( r_mesh[d] >= 0 );
	}
	r[d] = (double)cell[d] * dx - residue[d]; 
    }
}
inline void Relative_coord_OBL(const int *cell
			       ,const int *x_int
			       ,const double *residue
			       ,const int &sw_in_cell
			       ,const int Nlattice[DIM]
			       ,const double dx
			       ,int *r_mesh
			       ,double *r){
 
    double signY = x_int[1] + cell[1];
    r_mesh[1] = (x_int[1] + cell[1] + Nlattice[1]) % Nlattice[1];
    signY -= r_mesh[1];
    int sign = (int)signY;
    if (!(sign == 0)) {
	sign  = sign/abs(sign);
    }

    r_mesh[0] = (int)fmod(x_int[0] 
			  - sign*degree_oblique*Nlattice[1] + residue[0]/dx +
			  cell[0] + 4*Nlattice[0], Nlattice[0]);
    double x_oblique_residue =
	fmod(x_int[0] - 
	     sign*degree_oblique*Nlattice[1] + residue[0]/dx +
	     cell[0] + 4*Nlattice[0], Nlattice[0]) - r_mesh[0];
    
    r_mesh[2] = (x_int[2] + cell[2] + Nlattice[2]) % Nlattice[2];   

    for(int d=0;d<DIM;d++){
	{
	    assert( r_mesh[d] < Nlattice[d] );
	    assert( r_mesh[d] >= 0 );
	}
    }
    r[0] = (double)cell[0] * dx - x_oblique_residue*dx;
    r[1] = (double)cell[1] * dx - residue[1];
    r[2] = (double)cell[2] * dx - residue[2];   
}
inline int Relative_coord_check_stepover_Y(const int *cell
					   ,const int *x_int
					   ,const double *residue
					   ,const int &sw_in_cell
					   ,const int Nlattice[DIM]
					       ,const double dx
					   ,int *r_mesh
					   ,double *r){
    
    double signY = x_int[1] + cell[1];
    r_mesh[1] = (x_int[1] + cell[1] + Nlattice[1]) % Nlattice[1];
    signY -= r_mesh[1];
    int sign = (int)signY;
    if (!(sign == 0)) {
	sign  = sign/abs(sign);
    }

    r_mesh[0] = (int)fmod(x_int[0] 
			  - sign*degree_oblique*Nlattice[1] + residue[0]/dx +
			  cell[0] + 4*Nlattice[0], Nlattice[0]);
    double x_oblique_residue =
	fmod(x_int[0] - 
	     sign*degree_oblique*Nlattice[1] + residue[0]/dx +
	     cell[0] + 4*Nlattice[0], Nlattice[0]) - r_mesh[0];
    
    r_mesh[2] = (x_int[2] + cell[2] + Nlattice[2]) % Nlattice[2];   

    for(int d=0;d<DIM;d++){
	{
	    assert( r_mesh[d] < Nlattice[d] );
	    assert( r_mesh[d] >= 0 );
	}
    }
    r[0] = (double)cell[0] * dx - x_oblique_residue*dx;
    r[1] = (double)cell[1] * dx - residue[1];
    r[2] = (double)cell[2] * dx - residue[2];

    return sign;
}

/*!
  \brief Returns the orthonormal spherical coordinate basis vectors (in space
  coordinates) for a given point
  \param[in] x distance vector of the target point, measured in space
  frame, with particle center as origin
  \param[out] r unit radial basis vector
  \param[out] theta unit tangential (polar) basis vector
  \param[out] phi unit tangential (azimuthal) basis vector
  \param[out] norm_r distance from center to the target point
  \param[out] theta_angle polar angle [0,+pi], measured with respect to the
  swimming axis of the particle
  \param[out] phi_angle azimuthal angle [-pi,+pi], measured with respect to the
  x (e1) particle axis
  \param[in] p reference particle data
 */
inline void Squirmer_coord(const double *x, double *r, double *theta, double *phi, 
		      double &norm_r, double &theta_angle, double &phi_angle, const Particle &p){
  double ex[DIM] = {1.0, 0.0, 0.0};
  double ey[DIM] = {0.0, 1.0, 0.0};
  double ez[DIM] = {0.0, 0.0, 1.0};
  double *e1, *e2, *e3;

  // determine janus (e3) axis
  int dmy_axis = janus_axis[p.spec];
  if(dmy_axis == z_axis){
    e1 = ex;
    e2 = ey;
    e3 = ez;
  }else if(dmy_axis == y_axis){
    e1 = ez;
    e2 = ex;
    e3 = ey;
  }else if(dmy_axis == x_axis){
    e1 = ey;
    e2 = ez;
    e3 = ex;
  }else{
    fprintf(stderr, "Error: Invalid Janus axis for Spherical_coord\n");
    exit_job(EXIT_FAILURE);
  }

  //r normal vector
  rigid_body_rotation(r, x, p.q, SPACE2BODY);
  norm_r = sqrt( SQ(r[0]) + SQ(r[1]) + SQ(r[2]) );
  assert(positive_mp(norm_r));
  double dmy_norm = 1.0/norm_r;
  for(int d = 0; d < DIM; d++){
    r[d] *= dmy_norm;
  }

  //phi vector
  v_cross(phi, e3, r);  

  double dot_e3_r = v_inner_prod(e3, r);
  dmy_norm = sqrt(SQ(phi[0]) + SQ(phi[1]) + SQ(phi[2]));

  // No tangential velocity at the janus poles (set null tangent vectors) !
  if(non_zero_mp(dmy_norm) && less_than_mp(ABS(dot_e3_r), 1.0)){

    //phi normal vector
    dmy_norm = 1.0/dmy_norm;
    for(int d = 0; d < DIM; d++){
      phi[d] *= dmy_norm;
    }
    
    //theta normal vector
    theta_angle = acos(dot_e3_r);
    v_cross(theta, phi, r);
    dmy_norm = 1.0/sqrt(SQ(theta[0]) + SQ(theta[1]) + SQ(theta[2]));
    for(int d = 0; d < DIM; d++){
      theta[d] *= dmy_norm;
    }    
    
    //phi angle
    {
      double r12[DIM];
      for(int d = 0; d < DIM; d++){
        r12[d] = r[d] - e3[d] * dot_e3_r;
      }
      double dot_e1_r12 = v_inner_prod(e1, r12);
      double dot_e2_r12 = v_inner_prod(e2, r12);
      phi_angle = atan2(dot_e2_r12, dot_e1_r12);
    }
    
    rigid_body_rotation(r, p.q, BODY2SPACE);
    rigid_body_rotation(theta, p.q, BODY2SPACE);
    rigid_body_rotation(phi, p.q, BODY2SPACE);

  }else{ // r parallel to janus axis (tangent vectors not uniquely defined)
    norm_r = sqrt( SQ(x[0]) + SQ(x[1]) + SQ(x[2]) );
    dmy_norm = 1.0 / norm_r;
    for(int d = 0; d < DIM; d++){
      r[d] = x[d]*dmy_norm;
      phi[d] = theta[d] = 0.0;
    }

    if(dot_e3_r > 0.0){ // parallel
      theta_angle = 0.0;
      phi_angle = 0.0;
    }else{ // anti-parallel 
      theta_angle = M_PI;
      phi_angle = 0.0;
    }
  }
}


/*!
  \brief Compute linear velocity due to rigid body rotation (with \c ROTATION OFF)
  \details Trivial since rotation of the particle is ignored
  \f[
  \vec{v} = 0
  \f]
  \param[in] omega angular velocity 
  \param[in] r radial distance vector
  \param[out] v linear velocity due to rotation
 */
inline void Angular2v_rot_off(const double *omega, const double *r, double *v){
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
}
/*!
  \brief Compute linear velocity of a given point due to rigid body rotation (with \c ROTATION ON)
  \details
  \f[
  \vec{v} = \vec{\omega}\times\vec{r}
  \f]
  \param[in] omega angular velocity
  \param[in] r radial distance vector to given point
  \param[out] v linear velocity due to rotation
 */
inline void Angular2v_rot_on(const double *omega, const double *r, double *v){
    v[0] = omega[1] * r[2] - omega[2] * r[1];
    v[1] = omega[2] * r[0] - omega[0] * r[2];
    v[2] = omega[0] * r[1] - omega[1] * r[0];
}

/*!
  \brief Reset scalar field to a given value over specified domain
  \details
  \f[
  \phi(\vec{r}) = a, \qquad \vec{r}\in \mathcal{D}
  \f]
  \param[out] phi scalar field to reset
  \param[in] nx x domain \f$0\ldots N_x - 1\f$
  \param[in] ny y domain \f$0\ldots N_y - 1\f$
  \param[in] nz z domain \f$0\ldots N_z - 1\f$
  \param[in] value scalar value used to reset field
 */
inline void Reset_phi_primitive(double *phi
				,const int nx
				,const int ny
				,const int nz
				,const double &value
    ){
    int im;
#pragma omp parallel for schedule(dynamic, 1) private(im)
    for(int i=0;i<nx;i++){
	for(int j=0;j<ny;j++){
	    for(int k=0;k<nz;k++){
		im =(i*NY*NZ_)+(j*NZ_)+k; 
		phi[im] = value;
	    }
	}
    }
}

/*!
  \brief Reset scalar field to a given value over the whole domain
  \details \f[\phi(\vec{r}) = a\f]
  \param[out] phi scalar field to reset
  \param[in] value scalar value used to reset field
 */
inline void Reset_phi(double *phi, const double value = 0.0){
  Reset_phi_primitive(phi, NX, NY, NZ_, value);
}

/*!
  \brief Reset scalar and vector field to zero over the whole domain
  \details 
  \f{eqnarray*}{
  \phi(\vec{r}) &=& 0 \\
  \vec{u}(\vec{r}) &=& (0, 0, 0)
  \f}
  \param[out] phi scalar field to reset
  \param[out] up vector field to reset
 */
inline void Reset_phi_u(double *phi, double **up){
    int im;
//#pragma omp parallel for schedule(dynamic, 8) private(im)
#pragma omp parallel private(im)
    {
#pragma omp for nowait schedule(dynamic, 1)
	for(int i=0;i<NX;i++){
	    for(int j=0;j<NY;j++){
		for(int k=0;k<NZ_;k++){
		    im =(i*NY*NZ_)+(j*NZ_)+k; 
		    phi[im] = 0.; 
		}
	    }
	}
#pragma omp for nowait schedule(dynamic, 1)
	for(int i=0;i<NX;i++){
	    for(int j=0;j<NY;j++){
		for(int k=0;k<NZ_;k++){
		    im =(i*NY*NZ_)+(j*NZ_)+k; 
		    up[0][im] = 0.;
		}
	    }
	}
#pragma omp for nowait schedule(dynamic, 1)
	for(int i=0;i<NX;i++){
	    for(int j=0;j<NY;j++){
		for(int k=0;k<NZ_;k++){
		    im =(i*NY*NZ_)+(j*NZ_)+k; 
		    up[1][im] = 0.;
		}
	    }
	}
#pragma omp for schedule(dynamic, 1)
	for(int i=0;i<NX;i++){
	    for(int j=0;j<NY;j++){
		for(int k=0;k<NZ_;k++){
		    im =(i*NY*NZ_)+(j*NZ_)+k; 
		    up[2][im] = 0.;
		}
	    }
	}
    }
}

inline void Reset_u(double **up){
  int im;
#pragma omp parallel private(im)
  {
#pragma omp for nowait schedule(dynamic, 1)
    for(int i = 0; i < NX; i++){
      for(int j = 0; j < NY; j++){
	for(int k = 0; k < NZ_; k++){
	  im = (i * NY * NZ_) + (j * NZ_) + k;
	  up[0][im] = 0.0;
	}
      }
    }/* end omp for up[0] */
    
#pragma omp for nowait schedule(dynamic, 1)
    for(int i = 0; i < NX; i++){
      for(int j = 0; j < NY; j++){
	for(int k = 0; k < NZ_; k++){
	  im = (i * NY * NZ_) + (j * NZ_) + k;
	  up[1][im] = 0.0;
	}
      }
    }/* end omp for up[1] */
  
#pragma omp for nowait schedule(dynamic, 1)
    for(int i = 0; i < NX; i++){
      for(int j = 0; j < NY; j++){
	for(int k = 0; k < NZ_; k++){
	  im = (i * NY * NZ_) + (j * NZ_) + k;
	up[2][im] = 0.0;
	} 
      }
    }/* end omp for up[2] */
  }
}

#endif
