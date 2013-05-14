/*!
  \file fluid_solver.cxx
  \author Y. Nakayama
  \date 2006/08/15
  \version 1.2
  \brief Solver for Navier-Stokes equations 
 */
#include "fluid_solver.h"

double *Pressure;
double **Shear_force;
double **Shear_force_k;
double **f_ns0;
double **f_ns1;

////////// inline functions
// Shear_Navier_Stokes
const int Max_mean_shear_mode_PBC = 9;//HNY/9;
void Mean_shear_sustaining_force_PBC_OBL(double **u){
    for(int d = 0; d < DIM; d++){
	Reset_phi(Shear_force[d]);
	Reset_phi(Shear_force_k[d]);
    }


    const double dmy0= -2./SQ(M_PI) * LY* NY * NX * NZ * 0. * .5;
    for(int m = 1; m < 2*Max_mean_shear_mode_PBC + 2; m++){
	
	int j = m;

	Shear_force_k[0][j*NZ_] += -degree_oblique*u[1][j*NZ_] - u[0][j*NZ_];
	Shear_force_k[0][j*NZ_ + 1] += -degree_oblique*u[1][j*NZ_ + 1] - u[0][j*NZ_ + 1];
	u[0][j*NZ_] = -degree_oblique*u[1][j*NZ_];
	u[0][j*NZ_ + 1] = -degree_oblique*u[1][j*NZ_ + 1];

    }

    U_k2u(Shear_force_k);
    for(int i = 0; i < NX; i++){
	for(int j = 0; j < NY; j++){
	    for(int k = 0; k < NZ; k++){
		
		Shear_force[0][i*NY*NZ_ + j*NZ_ + k] += Shear_force_k[0][i*NY*NZ_ + j*NZ_ + k]; 

	    }
	}
    }
}
inline void Mean_shear_sustaining_yforce_PBC(double **zeta, double uk_dc[DIM], double **Shear_force,const CTime &jikan){
    
    for(int d=0;d<DIM;d++){
	Reset_phi(Shear_force[d]);
    }
    
    static const double qLY=LY/4.;
    static int j1=NY/4; 
    static int j2=3*NY/4; 
    static double dmy2=Shear_rate*qLY*DX;	
    
    Zeta_k2u(zeta, uk_dc, u);
    
    double dmy_1;    
    double dmy_2;    
    int im1;
    int im2;
#pragma omp parallel for schedule(dynamic, 1) private(dmy_1, dmy_2, im1, im2) 
    for(int i=0; i<NX; i++){
	for(int k=0; k<NZ; k++){
	    im1=(i*NY*NZ_)+(j1*NZ_)+k;
	    im2=(i*NY*NZ_)+(j2*NZ_)+k;
	    dmy_1=u[0][im1];
	    dmy_2=u[0][im2];
	    u[0][im1] = -dmy2; 		 
	    u[0][im2] = dmy2; 		 
	    Shear_force[0][im1] = u[0][im1] - dmy_1; 
	    Shear_force[0][im2] = u[0][im2] - dmy_2; 
	}
    }
    U2zeta_k(zeta, uk_dc, u);
}

inline void Mean_shear_sustaining_kforce_PBC(double **zeta, double uk_dc[DIM], double **Shear_force,const CTime &jikan){
 
  for(int d=0;d<DIM;d++){
    Reset_phi(Shear_force[d]);
  }
  
  static const double qLY=LY/4.;
  static const double HLY_shear=LY/2.;
  double time=jikan.ts*jikan.dt_fluid;
 
  Zeta_k2u(zeta, uk_dc, u);
  
  double dmy2;
  double dmy3;
  int im;
#pragma omp parallel for schedule(dynamic, 1) private(dmy2, dmy3, im)
  for(int i=0; i<NX; i++){
   for(int j=0; j<NY; j++){
	for(int k=0; k<NZ; k++){
		 im=(i*NY*NZ_)+(j*NZ_)+k;
         dmy3=u[0][im];
        
         if(int(j/qLY) == 0){
		dmy2 = -Shear_rate*j*DX;
          }else  if(int(j/qLY) == 3){
		dmy2 = Shear_rate*LY - Shear_rate*j*DX ;
	      }else{	
		dmy2 = -Shear_rate*HLY_shear + Shear_rate*j*DX;
          }

	    u[0][im] = dmy2*cos(Shear_frequency*time); 		 

	    Shear_force[0][im] = u[0][im] - dmy3; 
	 }
    }
   }
    U2zeta_k(zeta, uk_dc, u);
}

// ion ??
inline void Field_solver_Euler(const int &dim, double **zeta_k, const CTime &jikan, double **rhs, const Index_range &ijk_range){
int im;
  for(int d=0; d<dim; d++){
#pragma omp parallel for schedule(dynamic, 1) private(im)
    for(int i=ijk_range.istart; i<=ijk_range.iend; i++){
      for(int j=ijk_range.jstart; j<=ijk_range.jend; j++){
    	for(int k=ijk_range.kstart; k<=ijk_range.kend; k++){
	  im=(i*NY*NZ_)+(j*NZ_)+k;
	  zeta_k[d][im] += (jikan.dt_fluid * rhs[d][im]);
	}
      }
    }
  }
}
inline void Rhs_NS(double **zeta
		   ,double uk_dc[DIM]
		   ,double **rhs
		   ,const Index_range *ijk_range
		   ,const int &n_ijk_range
		   ,Particle *p
		   ){
  Zeta_k2advection_k(zeta, uk_dc, f_ns0);
  for(int n=0;n<n_ijk_range;n++){
    Add_zeta_viscous_term(zeta, rhs, ijk_range[n]);
  }
}

////////// inline functions end
void Mem_alloc_NS_solver(void){
    Pressure = alloc_1d_double(NX*NY*NZ_);
    f_ns0 = (double **) malloc(sizeof (double *)*DIM-1);
    f_ns1 = (double **) malloc(sizeof (double *)*DIM-1);

    Shear_force = (double **) malloc(sizeof (double *)*DIM);
    Shear_force_k = (double **) malloc(sizeof (double *)*DIM);
    for(int d=0;d<DIM-1;d++){
	f_ns0[d] = alloc_1d_double(NX*NY*NZ_);
	f_ns1[d] = alloc_1d_double(NX*NY*NZ_);
    }
    for(int d=0;d<DIM;d++){
	Shear_force[d] = alloc_1d_double(NX*NY*NZ_);
	Shear_force_k[d] = alloc_1d_double(NX*NY*NZ_);
    }
}

// Navior-Stokes 
////
inline void Rhs_solvent(double **u_solvent
			,double **rhs_solvent
			 ,double rhs_uk_dc[DIM]
			 ){
  

  U2advection_k(u_solvent, rhs_solvent); 
  for(int d=0; d<DIM; d++){
    rhs_uk_dc[d] = 0.0;
  }
}

inline void Rhs_NS_solute(Particle *p
			   ,double **zeta_k
			   ,double uk_dc[DIM]
			   ,double **u // working memory
			   ,double **concentration_k
			   ,double ** rhs_ns
			   ,double **rhs_solute
			   ,const Index_range *ijk_range
			   ,const int &n_ijk_range
			   ,double **solute_flux // working memory
			   ,double **grad_potential
			   ,double **surface_normal // working memory
			   ,double rhs_uk_dc[DIM] 
			   ){

  //Truncate_vector_two_third_rule(zeta_k, DIM-1);

  for(int d=0; d<(DIM-1); d++){
    Truncate_two_third_rule_ooura(zeta_k[d]);
  }

  Zeta_k2u(zeta_k, uk_dc, u);
  Make_surface_normal(surface_normal, p);

for(int n=0; n<N_spec; n++){
    //Truncate_two_third_rule(concentration_k[n]);
    Truncate_two_third_rule_ooura(concentration_k[n]);
    Diffusion_flux_single(solute_flux,concentration_k[n],Onsager_coeff[n], rhs_ns[0]);
    A_k2a_out(concentration_k[n], rhs_solute[n]);
    Solute_solver_rhs_nonlinear_x_single(grad_potential
					 ,rhs_solute[n]
					 ,solute_flux
					 ,Valency_e[n]
					 ,Onsager_coeff[n]
					 );
    Solute_impermeability(p, solute_flux, surface_normal);
    Add_advection_flux(solute_flux, u, rhs_solute[n]);
    U2u_k(solute_flux);
    U_k2divergence_k(solute_flux, rhs_solute[n]);

  }

  Rhs_solvent(u, rhs_ns, rhs_uk_dc);
  for(int n=0;n<n_ijk_range;n++){
    Add_zeta_viscous_term(zeta_k, rhs_ns, ijk_range[n]);
  }
}
inline void Add_constant_field_k(double **grad_potential_k
				 ,double e_ext[DIM]
				 ,const CTime &jikan
				 ){
  static const double nxnynz = (double)(NX*NY*NZ);
  double amp = nxnynz;
  if(AC){
    amp = sin( Angular_Frequency * jikan.time);
  }

  for(int d=0;d<DIM;d++){
    grad_potential_k[d][0] -= e_ext[d] * amp;
  }
}

inline void Rhs_NS_Nernst_Planck(Particle *p
				  ,double **zeta
				  ,double uk_dc[DIM]
				  ,double **u // working memory
				  ,double **concentration_k
				  ,double **rhs_ns
				  ,double **rhs_solute
				  ,const Index_range *ijk_range
				  ,const int &n_ijk_range
				  ,double **solute_flux // working memory
				  ,double **grad_potential // working memory
				  ,double **surface_normal // working memory
				  ,double rhs_uk_dc[DIM] 
				  ,const CTime &jikan
				  ){
  { // potential gradient in real space

    Conc_k2charge_field(p, concentration_k, u[0], u[1], u[2]);
    A2a_k(u[0]);
    Charge_field_k2Coulomb_potential_k_PBC(u[0]);

    //Truncate_two_third_rule(u[0]);
    Truncate_two_third_rule_ooura(u[0]);

    A_k2da_k(u[0], grad_potential);
    Add_constant_field_k(grad_potential, E_ext, jikan);
    U_k2u(grad_potential);
  }

  Rhs_NS_solute(p, zeta, uk_dc, u
		 ,concentration_k, rhs_ns, rhs_solute
		 ,ijk_range, n_ijk_range
		 ,solute_flux, grad_potential
		 ,surface_normal
		 , rhs_uk_dc
		 );

}
void Ion_diffusion_solver_Euler(double **zeta
				,const CTime &jikan
				,double uk_dc[DIM]
				,double **concentration_k
				,Particle *p
				,const Index_range *ijk_range
				,const int &n_ijk_range
				){
  double dc_rhs[DIM]; 
  {
    Rhs_NS_Nernst_Planck(p, zeta, uk_dc, u, concentration_k, f_ns0, Concentration_rhs0, ijk_range, n_ijk_range, up, f_particle, Surface_normal, dc_rhs, jikan);
  }
  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Euler(N_spec, concentration_k, jikan, Concentration_rhs0, ijk_range[n]);
  }
  {
    //double rescale_factor[N_spec];
    double *rescale_factor = new double[N_spec];
    Rescale_solute(rescale_factor
		   ,Total_solute
		   ,concentration_k, p, up[0], up[1]);
    delete [] rescale_factor;
  }
}

void NSsolute_solver_Euler(double **zeta
			   ,const CTime &jikan
			   ,double uk_dc[DIM]
			   ,double **concentration_k
			   ,Particle *p
			   ,const Index_range *ijk_range
			   ,const int &n_ijk_range
			   ){
    double dc_rhs[DIM]; 

    Rhs_NS_Nernst_Planck(p, zeta, uk_dc, u, concentration_k, f_ns1, Concentration_rhs0, ijk_range, n_ijk_range, up, f_particle, Surface_normal, dc_rhs, jikan);

  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Euler(DIM-1, zeta, jikan, f_ns1, ijk_range[n]);
    Field_solver_Euler(N_spec, concentration_k, jikan, Concentration_rhs0, ijk_range[n]);
  }

  for(int d=0;d<DIM;d++){
    uk_dc[d] += jikan.dt_fluid *  dc_rhs[d];
  }
}

void NS_solver_slavedEuler(double **zeta, const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p){
 
  //advection term on the rhs of NS equation (with minus sign)
  Zeta_k2advection_k(zeta, uk_dc, f_ns0);

  double dmy0 = -NU*jikan.dt_fluid; 
  const double dmy1 = 1./NU;
  int im;
  double dmy;
  for(int n=0;n<n_ijk_range;n++){
#pragma omp parallel for schedule(dynamic, 1) private(im, dmy) 
      for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	  for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	  for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
		im=(i*NY*NZ_)+(j*NZ_)+k;
	    dmy = exp(dmy0*K2[im])-1.; 
	        zeta[0][im] += 
	        dmy * (zeta[0][im] 
		     - dmy1 * IK2[im] * f_ns0[0][im]);
	        zeta[1][im] += 
	        dmy * (zeta[1][im] 
		     - dmy1 * IK2[im] * f_ns0[1][im]);
	  }
	}
      }
    }
}

// Shear_Navier_Stokes
void NS_solver_slavedEuler_Shear_PBC(double **zeta
				     ,const CTime &jikan
				     ,double uk_dc[DIM]
				     ,const Index_range *ijk_range
				     ,const int &n_ijk_range
				     ,Particle *p
				     ,double **force
    ){
    
    Zeta_k2advection_k(zeta, uk_dc, f_ns0);
    
    double dmy0 = -NU*jikan.dt_fluid; 
    const double dmy1 = 1./NU;
    int im;
    double dmy; 
    for(int n=0;n<n_ijk_range;n++){
#pragma omp parallel for schedule(dynamic, 1) private(dmy, im) 
	for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	    for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
		for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
		    im=(i*NY*NZ_)+(j*NZ_)+k;
		    dmy = exp(dmy0*K2[im])-1.; 
		    zeta[0][im] += 
			dmy * (zeta[0][im] 
			       - dmy1 * IK2[im] * f_ns0[0][im]);
		    zeta[1][im] += 
			dmy * (zeta[1][im] 
			       - dmy1 * IK2[im] * f_ns0[1][im]);
		}
	    }
	}
    }
    
    
    if(!Shear_AC){
	Mean_shear_sustaining_yforce_PBC(zeta, uk_dc, force, jikan);
    }else{ 
	Mean_shear_sustaining_kforce_PBC(zeta, uk_dc, force, jikan);
    }
}

void NS_solver_slavedEuler_Shear_OBL(double **zeta
				     ,const CTime &jikan
				     ,double uk_dc[DIM]
				     ,const Index_range *ijk_range
				     ,const int &n_ijk_range
				     ,Particle *p
				     ,double **force
    ){
    
    Zeta_k2advection_k_OBL(zeta, uk_dc, f_ns0);
    
    double dmy0 = -NU*jikan.dt_fluid; 
    const double dmy1 = 1./NU;
    int im;
    double dmy; 
    for(int n=0;n<n_ijk_range;n++){
#pragma omp parallel for schedule(dynamic, 1) private(dmy, im) 
	for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	    for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
		for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
		    im=(i*NY*NZ_)+(j*NZ_)+k;
		    dmy = exp(dmy0*K2[im])-1.; 
		    zeta[0][im] += 
			dmy * (zeta[0][im] 
			       - dmy1 * IK2[im] * f_ns0[0][im]);
		    zeta[1][im] += 
			dmy * (zeta[1][im] 
			       - dmy1 * IK2[im] * f_ns0[1][im]);
		}
	    }
	}
    }
}
