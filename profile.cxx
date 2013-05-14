//
// $Id: profile.cxx,v 1.9 2005/12/01 09:33:29 nakayama Exp $
//
#include "profile.h"

inline void Copy_sekibun_cell(const int &np_domain
			 ,int** &src
			 ,int** &dest
			 ){
  {
    dest = alloc_2d_int(np_domain,DIM);
    for(int n = 0;n<np_domain;n++){
      for(int d=0;d<DIM;d++){
	dest[n][d] = src[n][d];
      }
    }
    free_2d_int(src);
  }
}
void Particle_domain(
		     double (*profile_func)(const double &x, const double radius)
		     ,int &np_domain
		     ,int** &sekibun_cell
		     ,int &np_domain_extended
		     ,int** &sekibun_cell_extended
		     ,int &np_domain_interface
		     ,int** &sekibun_cell_interface
		     ,int &np_domain_interface_extended
		     ,int** &sekibun_cell_interface_extended
		     ,int &np_domain_exponential
		     ,int** &sekibun_cell_exponential
		     ){
  // sekibun_cell の初期化関数
  
  double dmy_xp[DIM];
  double dmy_xp_extended[DIM];
  const double dx_extended = 0.5 * DX;

  for(int d= 0; d< DIM; d++){
    //dmy_xp[d] = 0.5* DX; 
    //dmy_xp_extended[d] = 0.5* dx_extended;
    dmy_xp[d] = 0.;
    dmy_xp_extended[d] = 0.;
  }


  int n_mesh = (NX+1)*(NY+1)*(NZ+1);
  int n_mesh2 = (N2X+1)*(N2Y+1)*(N2Z+1);
  int **dmy_sekibun_cell = alloc_2d_int(n_mesh,DIM);
  int **dmy_sekibun_cell_interface = alloc_2d_int(n_mesh,DIM);
  int **dmy_sekibun_cell_extended = alloc_2d_int(n_mesh2,DIM);
  int **dmy_sekibun_cell_interface_extended = alloc_2d_int(n_mesh2,DIM);
  
  int **dmy_sekibun_cell_exponential = alloc_2d_int(n_mesh,DIM);

  int mesh = 0;
  int mesh_interface = 0;
  int mesh_exponential = 0;
  for(int i= -HNX; i<=HNX; i++){
    for(int j= -HNY; j<=HNY; j++){
      for(int k= -HNZ; k<=HNZ; k++){
	double dmy_x[DIM];
	dmy_x[0] = (double)i*DX; 
	dmy_x[1] = (double)j*DX; 
	dmy_x[2] = (double)k*DX; 
	
	double dmy = sqrt(
			  SQ(dmy_x[0]-dmy_xp[0])
			  +SQ(dmy_x[1]-dmy_xp[1])
			  +SQ(dmy_x[2]-dmy_xp[2])
			  );
	double dmy_phi=profile_func(dmy,RADIUS+DX);
	double dmy_phi_exponential = Phi_exponential(dmy);
	if(dmy_phi_exponential > 1.e-10){
	  dmy_sekibun_cell_exponential[mesh_exponential][0] = i;
	  dmy_sekibun_cell_exponential[mesh_exponential][1] = j;
	  dmy_sekibun_cell_exponential[mesh_exponential][2] = k;
	  mesh_exponential++;
	}
	if(dmy_phi > 0.0){
	  //if(dmy_phi > 1.e-10){
	  dmy_sekibun_cell[mesh][0] = i;
	  dmy_sekibun_cell[mesh][1] = j;
	  dmy_sekibun_cell[mesh][2] = k;
	  mesh++;
	  if(dmy_phi < 1.0){
	    dmy_sekibun_cell_interface[mesh_interface][0] = i;
	    dmy_sekibun_cell_interface[mesh_interface][1] = j;
	    dmy_sekibun_cell_interface[mesh_interface][2] = k;
	    mesh_interface++;
	  }
	}
      }
    }
  }
  np_domain = mesh;
  np_domain_interface = mesh_interface;
  np_domain_exponential = mesh_exponential;

  Copy_sekibun_cell(np_domain, dmy_sekibun_cell, sekibun_cell);
  Copy_sekibun_cell(np_domain_interface, dmy_sekibun_cell_interface, sekibun_cell_interface);
  Copy_sekibun_cell(np_domain_exponential, dmy_sekibun_cell_exponential, sekibun_cell_exponential);

  int mesh_extended=0;
  int mesh_interface_extended=0;
  for(int i= -NX; i<=NX; i++){
    for(int j= -NY; j<=NY; j++){
      for(int k= -NZ; k<=NZ; k++){
	double dmy_x[DIM];
	dmy_x[0] = (double)i*dx_extended;
	dmy_x[1] = (double)j*dx_extended;
	dmy_x[2] = (double)k*dx_extended;
	
	double dmy = Distance(dmy_x, dmy_xp_extended);
	double dmy_phi=profile_func(dmy,RADIUS+dx_extended);
	if(dmy_phi > 0.0){
	  dmy_sekibun_cell_extended[mesh_extended][0] = i;
	  dmy_sekibun_cell_extended[mesh_extended][1] = j;
	  dmy_sekibun_cell_extended[mesh_extended][2] = k;
	  mesh_extended++;
	  if(dmy_phi < 1.0){
	    dmy_sekibun_cell_interface_extended[mesh_interface_extended][0] = i;
	    dmy_sekibun_cell_interface_extended[mesh_interface_extended][1] = j;
	    dmy_sekibun_cell_interface_extended[mesh_interface_extended][2] = k;
	    mesh_interface_extended++;
	  }
	}
      }
    }
  }
  np_domain_extended = mesh_extended;
  np_domain_interface_extended = mesh_interface_extended;
  Copy_sekibun_cell(np_domain_extended
		    ,dmy_sekibun_cell_extended
		    ,sekibun_cell_extended);
  Copy_sekibun_cell(np_domain_interface_extended
		    ,dmy_sekibun_cell_interface_extended
		    ,sekibun_cell_interface_extended);

}
