//
// $Id: make_phi.h,v 1.16 2005/09/25 17:37:12 nakayama Exp $
//
#ifndef MAKE_PHI_H
#define MAKE_PHI_H

#include <assert.h> 
#include "input.h"
#include "variable.h"
#include "profile.h"

extern void (*Angular2v)(const double *omega, const double *r, double *v);
extern int NP_domain;
extern int NP_domain_extended;
extern int NP_domain_interface;
extern int NP_domain_interface_extended;
extern int **Sekibun_cell;
extern int **Sekibun_cell_extended;
extern int **Sekibun_cell_interface;
extern int **Sekibun_cell_interface_extended;

void Make_phi_u_advection(Value &phi, Value *up, Particle *p);
void Make_phi_particle(Value &phi
		       ,Particle *p
		       ,const double radius = RADIUS
		       );
void Make_phi_u_particle(Value &phi, Value *up, Particle *p);
void Make_phi_particle_extended(Value &phi,Particle *p);
void Make_phi_u_particle_extended(Value &phi,Value *up,Particle *p);
void Make_surface_normal(Value surface_normal[DIM]
			   ,const Particle *p);
void Make_phi_force(Value &phi, Value *fp, const Particle *p);
void Make_rho_field(Value &phi,const Particle *p);

void Make_phi_u_zwall(double u_wall[DIM]
		      ,Value &phi
		      ,Value *up
		      ,Particle *p
		      ,const int &SW_UP
		      ,const double &dx
		      ,const int &np_domain
		      ,int **sekibun_cell
		      ,const int Nlattice[DIM]
		      );

inline int Particle_cell(const double *xp
			 ,const double &dx
			 ,int *x_int
			 ,double *residue){
  const double idx = 1./dx;
    // r_i の左下のメッシュ座標を計算
    int sw_in_cell=0;
    for(int d=0; d<DIM; d++){
	double dmy = (xp[d] *idx);
	x_int[d] = (int)dmy;
	
	{
	    assert(xp[d] >= 0);
	    assert(xp[d] < L_particle[d]);
	}

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
inline void Reset_phi_primitive(Value &phi
				,const int nx
				,const int ny
				,const int nz
				,const double &value
				){
  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){
      for(int k=0;k<nz;k++){
	phi[i][j][k] = value;
      }
    }
  }
}
inline void Reset_phi(Value &phi_p, const double value = 0.0){
  Reset_phi_primitive(phi_p, NX, NY, NZ_, value);
}
inline void Reset_phi_extended(Value &phi
			       ,const double value = 0.0
			       ){
  Reset_phi_primitive(phi, N2X, N2Y, N2Z_, value);
}
inline void Reset_phi_u(Value &phi, Value *up){
    Reset_phi(phi);
    for(int d=0;d<DIM;d++){
	Reset_phi(up[d]);
    }
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

inline void Symmetrize_real(Value &f){
  for(int i=0;i<NX;i++){
    for(int j=1;j<NY_shear;j++){
      for(int k=0;k<NZ;k++){
	f[i][NY-j][k] = f[i][j][k];
      }
    }
  }
}
inline void Antisymmetrize_real(Value &f){
  for(int i=0;i<NX;i++){
    for(int k=0;k<NZ;k++){
      f[i][0][k] = 0.0;
      f[i][NY_shear][k] = 0.0;
    }
    for(int j=1;j<NY_shear;j++){
      for(int k=0;k<NZ;k++){
	f[i][NY-j][k] = -f[i][j][k];
      }
    }
  }
}
inline void Symmetrize_u(Value u[DIM]){
    Symmetrize_real(u[0]);
    Antisymmetrize_real(u[1]);
    Symmetrize_real(u[2]);
}
#endif
