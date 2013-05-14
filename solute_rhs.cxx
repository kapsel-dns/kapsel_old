//
// $Id: solute_rhs.cxx,v 1.8 2006/02/01 05:11:48 kin Exp $
//
#include "solute_rhs.h"

double *Valency_e;

int NP_domain_exponential;
int **Sekibun_cell_exponential;

double *Total_solute;
Value *Concentration;

Value *Concentration_rhs0;
Value *Concentration_rhs1;

Value Surface_normal[DIM];
//////

void Mem_alloc_solute(void){
  Valency_e = alloc_1d_double(N_spec);
  Onsager_coeff = alloc_1d_double(N_spec);
  Total_solute = alloc_1d_double(N_spec);
  Concentration = (Value *)malloc(sizeof(Value *) * N_spec);

  Concentration_rhs0= (Value *)malloc(sizeof(Value *) * N_spec);
  Concentration_rhs1= (Value *)malloc(sizeof(Value *) * N_spec);

  for(int n=0;n<N_spec;n++){
    Concentration[n] = alloc_3d_double(NX, NY, NZ_);
    Concentration_rhs0[n] = alloc_3d_double(NX, NY, NZ_);
    Concentration_rhs1[n] = alloc_3d_double(NX, NY, NZ_);
  }
  for(int d=0;d<DIM;d++){
    Surface_normal[d] = alloc_3d_double(NX, NY, NZ_);
  }
}


void Solute_impermeability(const Particle *p
			   ,Value solute_flux_x[DIM]
			   ,const Value surface_normal[DIM]
			   ){
  for(int n=0; n < Particle_Number; n++){
    double xp[DIM];
    for(int d=0;d<DIM;d++){
      xp[d] = p[n].x[d];
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

      double normal[DIM];
      double norm2 = 0.;
      for(int d=0;d<DIM;d++){
	normal[d] = surface_normal[d][r_mesh[0]][r_mesh[1]][r_mesh[2]];
	norm2 += SQ(normal[d]);
      } 
      if(norm2 > 0.){
	double dmy = 0.;
	for(int d=0;d<DIM;d++){
	  dmy += 
	    solute_flux_x[d][r_mesh[0]][r_mesh[1]][r_mesh[2]]
	    *normal[d];
	}
	dmy /= norm2;
	for(int d=0; d < DIM; d++ ){ 
	  solute_flux_x[d][r_mesh[0]][r_mesh[1]][r_mesh[2]]
	    -= normal[d] * dmy;
	}
      }
    }
  }
}
void Rescale_solute(double *rescale_factor
		    ,const double *total_solute
		    ,Value *conc_k
		    ,Particle *p
		    ,Value &phi_p // working memory
		    ,Value &conc_x // working memory
		    ){
  Reset_phi(phi_p);
  Make_phi_particle(phi_p,p);
  for(int n=0;n<N_spec;n++){
    A_k2a_out(conc_k[n], conc_x);
    double dmy = Count_single_solute(conc_x, phi_p);
    rescale_factor[n] = total_solute[n]/dmy;
    //fprintf(stderr, "bbb:%g %g\n", dmy, rescale);
    for(int i=0;i<NX;i++){
      for(int j=0;j<NY;j++){
	for(int k=0;k<NZ_;k++){
	  conc_k[n][i][j][k] *= rescale_factor[n];
	}
      }
    }
  }
}

double Count_single_solute(const Value &conc_x
			  ,const Value &phi_p 
			  ){
  // set phi(x) before calling 
  
  double dmy = 0.;
  for(int i=0;i<NX;i++){
    for(int j=0;j<NY;j++){
      for(int k=0;k<NZ;k++){
	dmy += (1.-phi_p[i][j][k]) * conc_x[i][j][k];
      }
    }
  }
  return dmy * DX3;
}

void Make_phi_exclude(Value &phi
		      ,Particle *p
		      ,const double &dx
		      ,const int &np_domain
		      ,int **sekibun_cell
		      ,const int Nlattice[DIM]
		      ){
  for(int n=0; n < Particle_Number; n++){
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
      double dmy = 0.;
      for(int d=0;d<DIM;d++){
	x[d] = r_mesh[d] * dx;
	dmy += SQ(x[d]-xp[d]);
      }
      double dmy_phi= Phi_exponential(dmy);
      phi[r_mesh[0]][r_mesh[1]][r_mesh[2]] += dmy_phi;
    }
  }
}

void Init_solute(Value &conc
		 ,Particle *p
		 ){
  {
    for(int n=0;n<N_spec;n++){
      Valency_e[n] = 1.;
      Onsager_coeff[n] = Onsager_solute_coeff;
    }
  }

  Reset_phi(conc);
  {
    for(int i=0;i<NX;i++){
      for(int j=0;j<NY;j++){
	for(int k=0;k<NZ;k++){
	  conc[i][j][k] = Mean_Bulk_concentration;
	}
      }
    }
  }
  {
    Reset_phi(phi);
    Make_phi_particle(phi,p);
    Total_solute[0] = Count_single_solute(conc, phi);
  }
  A2a_k(conc);
}

int NS_source = 0;
void Solute_solver_rhs_nonlinear_x_single(const Value grad_potential[DIM]
					  ,const Value &concentration_x
					  ,Value solute_flux[DIM]
					  ,const double &valency_e
					  ,const double &onsager_coeff
					  ,Value &omega_rhs
					  ){
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ; k++){
	double dmy_u[DIM];
	double dmy_grad_pot[DIM];
	double dmy_omega_rhs;
	for(int d=0;d<DIM;d++){
	  dmy_grad_pot[d] = grad_potential[d][i][j][k];
	}
	double dmy_conc = concentration_x[i][j][k];
	if(NS_source){
	  dmy_omega_rhs = -IRHO * dmy_conc * valency_e;
	}
	double dmy_conc_flux[DIM];
	for(int d=0;d<DIM;d++){
	  double dmy_interaction = dmy_grad_pot[d]*valency_e;
	  dmy_conc_flux[d] = dmy_conc * (
					 onsager_coeff * dmy_interaction
					 );
	}
	for(int d=0;d<DIM;d++){
	  solute_flux[d][i][j][k] += dmy_conc_flux[d];
	}
	if(NS_source){
	  omega_rhs[i][j][k] += dmy_omega_rhs;
	}
      }
    }
  }
} 
void Add_advection_flux(Value solute_flux[DIM]
			,const Value u_solvent[DIM]
			,const Value &concentration_x
			){
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ; k++){
	double dmy_conc = concentration_x[i][j][k];
	for(int d=0;d<DIM;d++){
	  solute_flux[d][i][j][k] += dmy_conc * -u_solvent[d][i][j][k];
	}
      }
    }
  }
}

void Diffusion_flux_single(Value diff_flux_x[DIM]
			    ,const Value &conc_k
			    ,const double &onsager_coeff
			    ,Value &dmy_value // working memory
			    ){
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ_; k++){
	dmy_value[i][j][k] = 
	  kBT * onsager_coeff * conc_k[i][j][k];
      }
    }
  }	      
  A_k2da_k(dmy_value, diff_flux_x);
  U_k2u(diff_flux_x);
}

