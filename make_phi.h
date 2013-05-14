//
// $Id: make_phi.h,v 1.3 2006/07/28 17:14:55 nakayama Exp $
//
#ifndef MAKE_PHI_H
#define MAKE_PHI_H

#include <assert.h> 
#include "input.h"
#include "variable.h"
#include "profile.h"

extern void (*Angular2v)(const double *omega, const double *r, double *v);
extern int NP_domain;
extern int NP_domain_interface;
extern int **Sekibun_cell;
extern int **Sekibun_cell_interface;

void Make_phi_u_advection(double *phi, double **up, Particle *p);
void Make_phi_particle(double *phi
		       ,Particle *p
		       ,const double radius = RADIUS
		       );
void Make_phi_u_particle(double *phi, double **up, Particle *p);
void Make_phi_particle_OBL(double *phi
			   ,Particle *p
			   ,const double radius = RADIUS
    );
void Make_phi_u_particle_OBL(double *phi, double **up, Particle *p);
void Make_surface_normal(double **surface_normal
			   ,const Particle *p);
void Make_rho_field(double *phi,Particle *p);

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

inline void Angular2v_rot_off(const double *omega, const double *r, double *v){
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
}
inline void Angular2v_rot_on(const double *omega, const double *r, double *v){
    v[0] = omega[1] * r[2] - omega[2] * r[1];
    v[1] = omega[2] * r[0] - omega[0] * r[2];
    v[2] = omega[0] * r[1] - omega[1] * r[0];
}
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
inline void Reset_phi(double *phi, const double value = 0.0){
  Reset_phi_primitive(phi, NX, NY, NZ_, value);
}
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

#endif
