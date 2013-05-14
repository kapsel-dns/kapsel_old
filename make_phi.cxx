//
// $Id: make_phi.cxx,v 1.14 2005/09/25 17:39:08 nakayama Exp $
//

#include "make_phi.h"

void (*Angular2v)(const double *omega, const double *r, double *v);

int NP_domain;
int NP_domain_extended;
int NP_domain_interface;
int NP_domain_interface_extended;
int **Sekibun_cell;
int **Sekibun_cell_extended;
int **Sekibun_cell_interface;
int **Sekibun_cell_interface_extended;

inline void Make_phi_u_primitive(Value &phi
				 ,Value *up
				 ,Particle *p
				 ,const int &SW_UP
				 ,const double &dx
				 ,const int &np_domain
				 ,int **sekibun_cell
				 ,const int Nlattice[DIM]
				 ,const double radius = RADIUS
				 ){
  static const double dv= POW3(dx);
  static double volume;
  
  for(int n=0; n < Particle_Number; n++){
    volume = 0.0;
    double xp[DIM],vp[DIM],omega_p[DIM];
    for(int d=0;d<DIM;d++){
      xp[d] = p[n].x[d];
      vp[d] = p[n].v[d];
      omega_p[d] = p[n].omega[d];
      
      if(0){
	fprintf(stderr, "ccc:%g %g %g\n"
		,p[n].x[0]/L[0]
		,p[n].x[1]/L[1]
		,p[n].x[2]/L[2]
		);
	assert(p[n].x[d] >= 0);
	assert(p[n].x[d] < L[d]);
	fprintf(stderr, "ddd:%g %g %g\n"
		,p[n].x[0]
		,p[n].x[1]
		,p[n].x[2]
		);
      }
    }
    
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell 
      = Particle_cell(xp, dx, x_int, residue);// {1,0} が返ってくる
    sw_in_cell = 1;
    int r_mesh[DIM];
    double r[DIM];
    for(int mesh=0; mesh < np_domain; mesh++){
      Relative_coord(sekibun_cell[mesh]
		     , x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
      double x[DIM];
      for(int d=0;d<DIM;d++){
	x[d] = r_mesh[d] * dx;
      }
      double dmy = Distance(x, xp);
      double dmy_phi= Phi(dmy,radius);
      phi[r_mesh[0]][r_mesh[1]][r_mesh[2]] += dmy_phi;
      volume += dmy_phi;
      
      if(SW_UP){
	double v_rot[DIM];
	Angular2v(omega_p, r, v_rot);
	for(int d=0;d<DIM;d++){
	  up[d][r_mesh[0]][r_mesh[1]][r_mesh[2]] += 
	    ( (vp[d]+v_rot[d]) * dmy_phi);
	}
      }
    }
    p[n].eff_mass_ratio 
      = 1.0;
    //= MASS[p[n].spec]/(RHO_particle[p[n].spec] * volume * dv);
  }
}


void Make_phi_u_advection(Value &phi, Value *up, Particle *p){
  static const double dv= POW3(DX);
  static double volume;
  
  int *nlattice;
  if(SW_EQ == Shear_Navier_Stokes){
    nlattice = Ns_shear;
  }else {
    nlattice = Ns;
  } 

  for(int n=0; n < Particle_Number; n++){
    volume = 0.0;
    double xp[DIM],vp[DIM];
    for(int d=0;d<DIM;d++){
      xp[d] = p[n].x[d];
      vp[d] = p[n].v[d];
      {
	assert(p[n].x[d] >= 0);
	assert(p[n].x[d] < L[d]);
      }
    }
    
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell 
      = Particle_cell(xp, DX, x_int, residue);// {1,0} が返ってくる
    sw_in_cell = 1;
    int r_mesh[DIM];
    double r[DIM];
    for(int mesh=0; mesh < NP_domain; mesh++){
      Relative_coord(Sekibun_cell[mesh], x_int, residue, sw_in_cell, nlattice, DX, r_mesh, r);
      double x[DIM];
      for(int d=0;d<DIM;d++){
	x[d] = r_mesh[d] * DX;
      }
      double dmy = Distance(x, xp);
      double dmy_phi= Phi(dmy);
      phi[r_mesh[0]][r_mesh[1]][r_mesh[2]] += dmy_phi;
      volume += dmy_phi;
      
      for(int d=0;d<DIM;d++){
	up[d][r_mesh[0]][r_mesh[1]][r_mesh[2]] += ( vp[d] * dmy_phi);
      }
    }
    p[n].eff_mass_ratio 
      = 1.0;
    //= MASS[p[n].spec]/(RHO_particle[p[n].spec] * volume * dv);
  }
  if(SW_EQ == Shear_Navier_Stokes){
    Symmetrize_real(phi);
    Symmetrize_u(up);
  }
}
void Make_phi_particle(Value &phi
		       ,Particle *p
		       ,const double radius
		       ){
  const int SW_UP = 0;
  Value *dmy_up;
  int *nlattice;
  if(SW_EQ == Shear_Navier_Stokes){
    nlattice = Ns_shear;
  }else {
    nlattice = Ns;
  } 
  Make_phi_u_primitive(phi, dmy_up, p, SW_UP,DX,NP_domain
		       ,Sekibun_cell
		       ,nlattice
		       ,radius);
  if(SW_EQ == Shear_Navier_Stokes){
    Symmetrize_real(phi);
  }
}
void Make_phi_u_particle(Value &phi
			 ,Value *up
			 ,Particle *p
			 ){
  const int SW_UP = 1;
  int *nlattice;
  if(SW_EQ == Shear_Navier_Stokes){
    nlattice = Ns_shear;
  }else {
    nlattice = Ns;
  } 
  Make_phi_u_primitive(phi, up, p, SW_UP,DX,NP_domain,Sekibun_cell,nlattice);
  if(SW_EQ == Shear_Navier_Stokes){
    Symmetrize_real(phi);
    Symmetrize_u(up);
  }
}
void Make_phi_particle_extended(Value &phi
				,Particle *p
				){
  const int SW_UP = 0;
  Value *dmy_up;
  Make_phi_u_primitive(phi, dmy_up, p
		       ,SW_UP,.5*DX,NP_domain_extended
		       ,Sekibun_cell_extended,N2s);
}
void Make_phi_u_particle_extended(Value &phi
			 ,Value *up
			 ,Particle *p
			 ){
  const int SW_UP = 1;
  Make_phi_u_primitive(phi, up, p
		       ,SW_UP,.5*DX,NP_domain_extended
		       ,Sekibun_cell_extended,N2s);
}

/////////////
void Make_surface_normal(Value surface_normal[DIM]
			   ,const Particle *p){
  
  for(int d=0;d<DIM;d++){
    Reset_phi(surface_normal[d]);
  }
  for(int n=0; n < Particle_Number; n++){
    double xp[DIM];
    for(int d=0;d<DIM;d++){
      xp[d] = p[n].x[d];
      
      assert(xp[d] >= 0);
      assert(xp[d] < L[d]);
    }
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell 
      = Particle_cell(xp, DX, x_int, residue);// {1,0} が返ってくる
    sw_in_cell = 1;
    int r_mesh[DIM];
    double r[DIM];
    for(int mesh=0; mesh < NP_domain; mesh++){
      Relative_coord(Sekibun_cell[mesh]
		     ,x_int, residue, sw_in_cell, Ns, DX, r_mesh, r);
      double x[DIM];
      for(int d=0;d<DIM;d++){
        x[d] = r_mesh[d] * DX;
      }
      double dmy_r = Distance(x, xp);
      double dmy = ABS(dmy_r - RADIUS);
      if(dmy < HXI){
	double ir = 1./dmy_r;
	for(int d=0;d<DIM;d++){
	  surface_normal[d][r_mesh[0]][r_mesh[1]][r_mesh[2]] 
	    += (x[d]-xp[d]) * ir;
	}
      }
    }
  }
}
/////////////
void Make_phi_force(Value &phi, Value *fp, const Particle *p){

  for(int n=0; n < Particle_Number; n++){
    double xp[DIM],f[DIM];
    double imass = IMASS[p[n].spec];
    for(int d=0;d<DIM;d++){
      xp[d] = p[n].x[d];
      f[d] = p[n].fr[d] + p[n].fv[d];
      
      {
	assert(p[n].x[d] >= 0);
	assert(p[n].x[d] < L[d]);
      }
    }
    
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell 
      = Particle_cell(xp, DX, x_int, residue);// {1,0} が返ってくる
    sw_in_cell = 1;
    int r_mesh[DIM];
    double r[DIM];
    for(int mesh=0; mesh < NP_domain; mesh++){
      Relative_coord(Sekibun_cell[mesh]
		     ,x_int, residue, sw_in_cell
		     ,Ns, DX, r_mesh, r);
      double x[DIM];
      for(int d=0;d<DIM;d++){
	x[d] = r_mesh[d] * DX;
      }
      double dmy = Distance(x, xp);
      double dmy_phi= Phi(dmy);
      phi[r_mesh[0]][r_mesh[1]][r_mesh[2]] += dmy_phi;
      
      for(int d=0;d<DIM;d++){
	fp[d][r_mesh[0]][r_mesh[1]][r_mesh[2]] += ( f[d] * dmy_phi * imass);
      }
    }
  }
}
/////////////
inline void Make_rho_field_primitive(Value &phi
				     ,const Particle *p
				     ,const double &dx
				     ,const int &np_domain
				     ,int **sekibun_cell
				     ,const int Nlattice[DIM]
				     ){
  static const double dv= POW3(dx);
  
  for(int n=0; n < Particle_Number; n++){
    double drho = RHO_particle[p[n].spec] - RHO;
    double xp[DIM];
    for(int d=0;d<DIM;d++){
      xp[d] = p[n].x[d];
      
    }
    
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell 
      = Particle_cell(xp, dx, x_int, residue);// {1,0} が返ってくる
    sw_in_cell = 1;
    int r_mesh[DIM];
    double r[DIM];
    for(int mesh=0; mesh < np_domain; mesh++){
      Relative_coord(sekibun_cell[mesh]
		     , x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
      double x[DIM];
      for(int d=0;d<DIM;d++){
	x[d] = r_mesh[d] * dx;
      }
      double dmy = Distance(x, xp);
      double dmy_phi= Phi(dmy) * drho;
      phi[r_mesh[0]][r_mesh[1]][r_mesh[2]] += dmy_phi;
    }
  }
  for(int i=0;i<Nlattice[0];i++){
    for(int j=0;j<Nlattice[1];j++){
      for(int k=0;k<Nlattice[2];k++){
	phi[i][j][k] += RHO; 
      }
    }
  }
}
void Make_rho_field(Value &phi
		    ,const Particle *p
		    ){
  int *nlattice;
  if(SW_EQ == Shear_Navier_Stokes){
    nlattice = Ns_shear;
  }else {
    nlattice = Ns;
  } 
  Make_rho_field_primitive(phi, p, DX, NP_domain,Sekibun_cell,nlattice);
  if(SW_EQ == Shear_Navier_Stokes){
    Symmetrize_real(phi);
  }
}

void Make_phi_u_zwall(double u_wall[DIM]
		      ,Value &phi
		      ,Value *up
		      ,Particle *p
		      ,const int &SW_UP
		      ,const double &dx
		      ,const int &np_domain
		      ,int **sekibun_cell
		      ,const int Nlattice[DIM]
		      ){
  const int k_radius_wall = (int)ceil(A_radius_wall + HXI/dx);

  for(int i=0;i<Nlattice[0];i++){
    for(int j=0;j<Nlattice[1];j++){
      {
	int k=0;
	double dmy_phi = Phi(0., A_radius_wall);
	phi[i][j][k] += dmy_phi;
      }
      for(int k=1;k<=k_radius_wall;k++){
	double dmy_r = (double)k*dx;
	double dmy_phi = Phi(dmy_r, A_radius_wall);
	phi[i][j][k] += dmy_phi;
	int k_mirror = Nlattice[2] - k;
	phi[i][j][k_mirror] += dmy_phi;
	if(SW_UP){
	  for(int d=0;d<DIM;d++){
	    double dmy_u = u_wall[d] * dmy_phi;
	    up[d][i][j][k] += dmy_u;
	    up[d][i][j][k_mirror] += dmy_u;
	  }
	}
      }
    }
  }
}
