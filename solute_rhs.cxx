/*!
  \file solute_rhs.cxx
  \brief Routines to compute the terms appearing in the right hand side of the solute advection diffusion equation
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
 */

#include "solute_rhs.h"

double *Valency_e;

double *Total_solute;

double **Surface_normal;
double **Concentration;
double **Concentration_rhs0;
double **Concentration_rhs1;

//////

void Mem_alloc_solute(void){
  Valency_e = alloc_1d_double(N_spec);
  Total_solute = alloc_1d_double(N_spec);

  Surface_normal= (double **)malloc(sizeof(double *) * DIM);

  Concentration= (double **)malloc(sizeof(double *) * N_spec);
  Concentration_rhs0= (double **)malloc(sizeof(double *) * N_spec);
  Concentration_rhs1= (double **)malloc(sizeof(double *) * N_spec);

  for(int n=0;n<N_spec;n++){
    Concentration[n]= alloc_1d_double(NX*NY*NZ_);
    Concentration_rhs0[n]= alloc_1d_double(NX*NY*NZ_);
    Concentration_rhs1[n]= alloc_1d_double(NX*NY*NZ_);
  }
  for(int d=0;d<DIM;d++){
    Surface_normal[d] = alloc_1d_double(NX*NY*NZ_);
  }
}

void Solute_impermeability(Particle *p
			   ,double **solute_flux_x
			   ,double **surface_normal
			   ){
    double xp[DIM];
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell;
    int r_mesh[DIM];
    double r[DIM];
    double normal[DIM];
    double norm2;
	double dmy;
#pragma omp parallel for schedule(dynamic, 1) private(xp,x_int,residue,sw_in_cell,r_mesh,r,normal,norm2,dmy)
  for(int n=0; n < Particle_Number; n++){
  // double xp[DIM];
    for(int d=0;d<DIM;d++){
      xp[d] = p[n].x[d];
    }

    //int x_int[DIM];
    //double residue[DIM];
    //int sw_in_cell 
    sw_in_cell 
      = Particle_cell(xp, DX, x_int, residue);// {1,0} が返ってくる
    sw_in_cell = 1;
    //int r_mesh[DIM];
    //double r[DIM];
    for(int mesh=0; mesh < NP_domain; mesh++){
      Relative_coord(Sekibun_cell[mesh]
		     ,x_int, residue, sw_in_cell, Ns, DX, r_mesh, r);

    //  double normal[DIM];
      //double norm2 = 0.;
      norm2 = 0.;
      for(int d=0;d<DIM;d++){
	normal[d] = surface_normal[d][(r_mesh[0]*NY*NZ_)+(r_mesh[1]*NZ_)+r_mesh[2]];
	norm2 += SQ(normal[d]);
      } 
      if(norm2 > 0.){
	//double dmy = 0.;
	dmy = 0.;
	for(int d=0;d<DIM;d++){
	  dmy += 
	    solute_flux_x[d][(r_mesh[0]*NY*NZ_)+(r_mesh[1]*NZ_)+r_mesh[2]]
	    *normal[d];
	}
	dmy /= norm2;
	for(int d=0; d < DIM; d++ ){ 
	  solute_flux_x[d][(r_mesh[0]*NY*NZ_)+(r_mesh[1]*NZ_)+r_mesh[2]]
	    -= normal[d] * dmy;
	}
      }
    }
  }
}

void Rescale_solute(double *rescale_factor
		    ,double *total_solute
		    ,double **conc_k
		    ,Particle *p
		    ,double *phi // working memory
		    ,double *conc_x // working memory
		    ){

  Reset_phi(phi);
  Make_phi_particle(phi,p);

  int im;
  double dmy;

  for(int n=0;n<N_spec;n++){
    A_k2a_out(conc_k[n], conc_x);
    //double dmy = Count_single_solute(conc_x, phi);
    dmy = Count_single_solute(conc_x, phi);
    rescale_factor[n] = total_solute[n]/dmy;
#pragma omp parallel for schedule(dynamic, 1) private(im)
    for(int i=0;i<NX;i++){
      for(int j=0;j<NY;j++){
	    for(int k=0;k<NZ_;k++){
		  im=(i*NY*NZ_)+(j*NZ_)+k;
	  conc_k[n][im] *= rescale_factor[n];
	}
      }
    }
  }
}

double Count_single_solute(double *conc_x
			  ,double *phi_p
			  ){
  // set phi(x) before calling 
  
  double dmy = 0.;
  int im;
#pragma omp parallel for schedule(dynamic, 1) reduction(+:dmy) private(im)
  for(int i=0;i<NX;i++){
    for(int j=0;j<NY;j++){
      for(int k=0;k<NZ;k++){
		//int im=(i*NY*NZ_) + (j*NZ_) +k;
		im=(i*NY*NZ_)+(j*NZ_)+k;
	    dmy += (1.-phi_p[im]) * conc_x[im];
      }
    }
  }
  return dmy * DX3;
}

void Solute_solver_rhs_nonlinear_x_single(double **grad_potential
					  ,double *concentration_x
					  ,double **solute_flux
					  ,double &valency_e
					  ,double &onsager_coeff
					  ){
  int im;
  double dmy_u[DIM];
  double dmy_grad_pot[DIM];
  double dmy_conc;
  double dmy_conc_flux[DIM];
  double dmy_interaction;
#pragma omp parallel for schedule(dynamic, 1) private(im,dmy_u,dmy_grad_pot,dmy_conc,dmy_conc_flux,dmy_interaction)
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ; k++){
	im=(i*NY*NZ_)+(j*NZ_)+k;
	//double dmy_u[DIM];
	//double dmy_grad_pot[DIM];

	for(int d=0;d<DIM;d++){
	  //dmy_grad_pot[d] = grad_potential[d][i][j][k];
	  dmy_grad_pot[d] = grad_potential[d][im];
	}
	//double dmy_conc = concentration_x[im];
	dmy_conc = concentration_x[im];

	//double dmy_conc_flux[DIM];
	for(int d=0;d<DIM;d++){
	  dmy_interaction = dmy_grad_pot[d]*valency_e;
	  //double dmy_interaction = dmy_grad_pot[d]*valency_e;
	  dmy_conc_flux[d] = dmy_conc * (
					 onsager_coeff * dmy_interaction
					 );
	}
	for(int d=0;d<DIM;d++){
	  solute_flux[d][im] += dmy_conc_flux[d];
	}
      }
    }
  }
} 
void Add_advection_flux(double **solute_flux
			,double **u_solvent
			,double *concentration_x
			){
  int im;
  double dmy_conc;
#pragma omp parallel for schedule(dynamic, 1) private(im,dmy_conc)
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ; k++){
		//int im=(i*NY*NZ_)+(j*NZ_)+k;
		im=(i*NY*NZ_)+(j*NZ_)+k;
	    dmy_conc = concentration_x[im];
	//double dmy_conc = concentration_x[im];
	for(int d=0;d<DIM;d++){
	  solute_flux[d][im] += dmy_conc * -u_solvent[d][im];
	}
      }
    }
  }
}

void Diffusion_flux_single(double **diff_flux_x
			    ,double *conc_k
			    ,double &onsager_coeff
			    ,double *dmy_value // working memory
			    ){
  int im;
#pragma omp parallel for schedule(dynamic, 1) private(im)
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ_; k++){
		//int im=(i*NY*NZ_)+(j*NZ_)+k;
		im=(i*NY*NZ_)+(j*NZ_)+k;
	dmy_value[im] = 
	  kBT * onsager_coeff * conc_k[im];
      }
    }
  }	      

  A_k2da_k(dmy_value, diff_flux_x);
  U_k2u(diff_flux_x);
}

