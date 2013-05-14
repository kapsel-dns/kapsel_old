//
// $Id: fluid_solver.cxx,v 1.44 2006/06/09 01:38:35 nakayama Exp $
//
#include "fluid_solver.h"

Value Pressure;
Value f_ns0[DIM-1];
Value f_ns1[DIM-1];

Value Shear_force[DIM];
////////// inline functions
const int Max_mean_shear_mode = 5;//HNY/9;
inline void Mean_shear_sustaining_force(Value zeta[DIM-1], double uk_dc[DIM], Value force[DIM]){
  for(int d=0;d<DIM;d++){
    Reset_phi(Shear_force[d]);
  }
  Zeta_k2u_k(zeta, uk_dc, u);
  const double dmy0= 4./SQ(M_PI) * LY_shear* HNY * NX * NZ * Shear_rate;
  for(int m=0;m<= Max_mean_shear_mode;m++){
    int j= 2*m+1;
    double dmy1 = dmy0/SQ(j);
    Shear_force[0][0][j][0] = dmy1 - u[0][0][j][0];
    Shear_force[0][0][j][1] =  - u[0][0][j][1];
    u[0][0][j][0] = dmy1;
    u[0][0][j][1] = 0.0;
  }
  U_k2zeta_k(u, zeta, uk_dc);
  Symmetrize_zetak(zeta);

  {
    //Solenoidal_uk(Shear_force);
  }
}
inline void Rhs_slippy_NS_viscosity(Value zeta[DIM-1], Value rhs[DIM-1], double force_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p){
  {
    for(int n=0;n<n_ijk_range;n++){
      for(int d=0;d<DIM-1;d++){
	for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	  for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	    for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	      up[d][i][j][k] = -(K2[i][j][k]* zeta[d][i][j][k]);
	    }
	  }
	}
      }
    }
    Truncate_vector_two_third_rule(up, DIM-1);
    double dmy[DIM]={0.,0.,0.};
    Zeta_k2u(up, dmy, u); 
  }

  {// making phi
    Reset_phi(phi);
    Make_phi_u_particle(phi, up, p);
  }
  double delta_nu = NU * (Nu_ratio -1.);
  for(int d=0;d<DIM;d++){
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ; k++){
	  double kinematic_viscosity = NU + delta_nu *phi[i][j][k];
	  u[d][i][j][k] *= kinematic_viscosity;
	}
      }
    }
  }
  double dmy[DIM];
  U2zeta_k(up, dmy, u);
  for(int d=0;d<DIM;d++){
    force_dc[d] += dmy[d];
  }
  { 
    for(int n=0;n<n_ijk_range;n++){
      for(int d=0;d<DIM-1;d++){
	for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	  for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	    for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	      rhs[d][i][j][k] += u[d][i][j][k];
	    }
	  }
	}
      }
    }
  }
}

inline void Field_solver_Euler(const int &dim, Value *zeta, const CTime &jikan, const Value *rhs, const Index_range &ijk_range){
  for(int d=0; d<dim; d++){
    for(int i=ijk_range.istart; i<=ijk_range.iend; i++){
      for(int j=ijk_range.jstart; j<=ijk_range.jend; j++){
    	for(int k=ijk_range.kstart; k<=ijk_range.kend; k++){
	  zeta[d][i][j][k] += (jikan.dt_fluid * rhs[d][i][j][k]);
	}
      }
    }
  }
}
inline void Field_solver_Heun(const int &dim, Value *zeta, const CTime &jikan, const Value *rhs0, const Value *rhs1, const Index_range &ijk_range){
  for(int d=0; d<dim; d++){
    for(int i=ijk_range.istart; i<=ijk_range.iend; i++){
      for(int j=ijk_range.jstart; j<=ijk_range.jend; j++){
    	for(int k=ijk_range.kstart; k<=ijk_range.kend; k++){
	  zeta[d][i][j][k] += (jikan.hdt_fluid 
			       * (rhs1[d][i][j][k] -rhs0[d][i][j][k]));
	}
      }
    }
  }
}
inline void Euler_solver_Euler(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p){
  Zeta_k2advection_k(zeta, uk_dc, f_ns0);
  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Euler(DIM-1, zeta, jikan, f_ns0, ijk_range[n]);
  }
}
inline void Euler_solver_Heun(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p){
  Zeta_k2advection_k(zeta, uk_dc, f_ns0);
  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Euler(DIM-1, zeta, jikan, f_ns0, ijk_range[n]);
  }

  Zeta_k2advection_k(zeta, uk_dc, f_particle);
  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Heun(DIM-1, zeta, jikan, f_ns0, f_particle, ijk_range[n]);
  }
}
inline void Rhs_NS(Value zeta[DIM-1]
		   ,const double uk_dc[DIM]
		   ,Value rhs[DIM-1]
		   ,const Index_range *ijk_range
		   ,const int &n_ijk_range
		   ,Particle *p
		   ){
  if(1){
    if(STOKES){
      SP_stokes_advection(zeta, uk_dc, rhs, ijk_range, n_ijk_range, p);
    }else {
      //SP_full_advection(zeta, uk_dc, rhs, ijk_range, n_ijk_range, p);
      Zeta_k2advection_k(zeta, uk_dc, rhs);
    }
  }else{
    for(int d=0;d<DIM-1;d++){
      Reset_phi(rhs[d]);
    }
  }
  for(int n=0;n<n_ijk_range;n++){
    Add_zeta_viscous_term(zeta, rhs, ijk_range[n]);
  }
}
////////// inline functions end
void Mem_alloc_NS_solver(void){
  Pressure = alloc_3d_double(NX, NY, NZ_);
  for(int d=0;d<DIM-1;d++){
    f_ns0[d] = alloc_3d_double(NX, NY, NZ_);
    f_ns1[d] = alloc_3d_double(NX, NY, NZ_);
  }
}


void Stokes_solver(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM]){
  static const double dmy = NU * jikan.dt_fluid;
  for(int d=0; d<DIM-1; d++){
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ_; k++){
	  zeta[d][i][j][k] *= exp(- dmy * K2[i][j][k] ); 
	}
      }
    }
  }
}

////////////////// SP Stokes advection
void SP_stokes_advection(Value zeta[DIM-1], const double uk_dc[DIM], Value rhs[DIM-1], const Index_range *ijk_range, const int &n_ijk_range, Particle *p){
  {// making phi
    Reset_phi_u(phi,up);
    Make_phi_u_advection(phi, up, p);
  }
  Solenoidal_u(up);
  {
    Truncate_vector_two_third_rule(zeta, DIM-1);
    
    Zeta_k2u(zeta, uk_dc, u); 
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ; k++){
	  double v[DIM];
	  double v_advection[DIM];
	  for(int d=0;d<DIM;d++){
	    v[d] = u[d][i][j][k];
	    v_advection[d] = up[d][i][j][k];
	  }
	  double sigma[DIM][DIM];
	  {///// advection of v
	    for(int m=0;m<DIM;m++){
	      for(int n=0;n<DIM;n++){
		sigma[m][n] = - v_advection[m] * v[n];
	      }
	    }
	  }
	  { 
	    up[0][i][j][k] = sigma[0][0] - sigma[2][2];
	    up[1][i][j][k] = sigma[0][1];
	    up[2][i][j][k] = sigma[0][2];
	    u[0][i][j][k] = sigma[1][0];
	    u[1][i][j][k] = sigma[1][1] - sigma[2][2];
	    u[2][i][j][k] = sigma[1][2];
	    rhs[0][i][j][k] = sigma[2][0];
	    rhs[1][i][j][k] = sigma[2][1];
	  }
	}
      }
    }
    {
      for(int d=0;d<DIM;d++){
	A2a_k(up[d]);
      }
      for(int d=0;d<DIM;d++){
	A2a_k(u[d]);
      }
      for(int d=0;d<DIM-1;d++){
	A2a_k(rhs[d]);
      }
    }
  }
  { 
    for(int n=0;n<n_ijk_range;n++){
      for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	  for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	    static double stress[DIM][DIM];
	    stress[0][0] = up[0][i][j][k];
	    stress[0][1] = up[1][i][j][k];
	    stress[0][2] = up[2][i][j][k];
	    stress[1][0] = u[0][i][j][k];
	    stress[1][1] = u[1][i][j][k];
	    stress[1][2] = u[2][i][j][k];
	    stress[2][0] = rhs[0][i][j][k];
	    stress[2][1] = rhs[1][i][j][k];
	    stress[2][2] = 0.0;
	    
	    {
	      static double ks[DIM];
	      ks[0] =KX_int[i][j][k] * WAVE_X;
	      ks[1] =KY_int[i][j][k] * WAVE_Y;
	      ks[2] =KZ_int[i][j][k] * WAVE_Z;
	      static double div_stress[DIM];
	      for(int d=0;d<DIM;d++){
		div_stress[d] = 0.0;
		for(int m=0;m<DIM;m++){
		  div_stress[d] += (ks[m] * stress[m][d]);
		}
	      }
	      up[0][i][j][k] =
		-(ks[1] * div_stress[2] - ks[2] * div_stress[1]);
	      up[1][i][j][k] =
		-(ks[2] * div_stress[0] - ks[0] * div_stress[2]);
	      up[2][i][j][k] =
		-(ks[0] * div_stress[1] - ks[1] * div_stress[0]);
	    }
	  }
	}
      }
    }
    Omega_k2zeta_k(up, rhs);
  }
}
void SP_full_advection(Value zeta[DIM-1], const double uk_dc[DIM], Value rhs[DIM-1], const Index_range *ijk_range, const int &n_ijk_range, Particle *p){
  {// making advection-velocity
    Truncate_vector_two_third_rule(zeta, DIM-1);
    Zeta_k2u(zeta, uk_dc, u);

    Reset_phi_u(phi,up);
    Make_phi_u_advection(phi, up, p);
    Make_f_particle_dt_sole(f_particle, u, up, phi);
    Add_f_particle(f_particle, u);
  }
  {
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ; k++){
	  double v[DIM];
	  double v_advection[DIM];
	  for(int d=0;d<DIM;d++){
	    v[d] = u[d][i][j][k];
	    v_advection[d] = f_particle[d][i][j][k];
	  }
	  double sigma[DIM][DIM];
	  {///// advection of v
	    for(int m=0;m<DIM;m++){
	      for(int n=0;n<DIM;n++){
		sigma[m][n] = - v_advection[m] * v[n];
	      }
	    }
	  }
	  { 
	    up[0][i][j][k] = sigma[0][0] - sigma[2][2];
	    up[1][i][j][k] = sigma[0][1];
	    up[2][i][j][k] = sigma[0][2];
	    u[0][i][j][k] = sigma[1][0];
	    u[1][i][j][k] = sigma[1][1] - sigma[2][2];
	    u[2][i][j][k] = sigma[1][2];
	    rhs[0][i][j][k] = sigma[2][0];
	    rhs[1][i][j][k] = sigma[2][1];
	  }
	}
      }
    }
    {
      for(int d=0;d<DIM;d++){
	A2a_k(up[d]);
      }
      for(int d=0;d<DIM;d++){
	A2a_k(u[d]);
      }
      for(int d=0;d<DIM-1;d++){
	A2a_k(rhs[d]);
      }
    }
  }
  { 
    for(int n=0;n<n_ijk_range;n++){
      for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	  for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	    static double stress[DIM][DIM];
	    stress[0][0] = up[0][i][j][k];
	    stress[0][1] = up[1][i][j][k];
	    stress[0][2] = up[2][i][j][k];
	    stress[1][0] = u[0][i][j][k];
	    stress[1][1] = u[1][i][j][k];
	    stress[1][2] = u[2][i][j][k];
	    stress[2][0] = rhs[0][i][j][k];
	    stress[2][1] = rhs[1][i][j][k];
	    stress[2][2] = 0.0;
	    
	    {
	      static double ks[DIM];
	      ks[0] =KX_int[i][j][k] * WAVE_X;
	      ks[1] =KY_int[i][j][k] * WAVE_Y;
	      ks[2] =KZ_int[i][j][k] * WAVE_Z;
	      static double div_stress[DIM];
	      for(int d=0;d<DIM;d++){
		div_stress[d] = 0.0;
		for(int m=0;m<DIM;m++){
		  div_stress[d] += (ks[m] * stress[m][d]);
		}
	      }
	      up[0][i][j][k] =
		-(ks[1] * div_stress[2] - ks[2] * div_stress[1]);
	      up[1][i][j][k] =
		-(ks[2] * div_stress[0] - ks[0] * div_stress[2]);
	      up[2][i][j][k] =
		-(ks[0] * div_stress[1] - ks[1] * div_stress[0]);
	    }
	  }
	}
      }
    }
    Omega_k2zeta_k(up, rhs);
  }
}
void NS_solver_Euler0(Value zeta[DIM-1]
		      ,const CTime &jikan
		      ,double uk_dc[DIM]
		      ,const Index_range *ijk_range
		      ,const int &n_ijk_range
		      ,Particle *p
		      ){
  Rhs_NS(zeta, uk_dc, f_ns0, ijk_range, n_ijk_range, p);
  if(kBT > 0.0){
    const double truncate_factor= (double)NX/(2.*TRN_X-1.)*(double)NY/(2.*TRN_Y-1.)*(double)NZ/(2.*TRN_Z-1.);
    Add_random_stress_term_NS(zeta, f_ns0, uk_dc, ijk_range, n_ijk_range, jikan, truncate_factor, p);
  }
  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Euler(DIM-1, zeta, jikan, f_ns0, ijk_range[n]);
  }
}
void NS_solver_Euler(Value zeta[DIM-1]
		     ,const CTime &jikan
		     ,double uk_dc[DIM]
		     ,const Index_range *ijk_range
		     ,const int &n_ijk_range
		     ,Particle *p
		     ){
  Rhs_NS(zeta, uk_dc, f_ns0, ijk_range, n_ijk_range, p);
  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Euler(DIM-1, zeta, jikan, f_ns0, ijk_range[n]);
  }
  int dmy=0;
  if(1){
    for(int n=0;n<n_ijk_range;n++){
      for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	  for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	    dmy++;
	  }
	}
      }
    }
    if(0){
      fprintf(stderr, "%d %d %g %g\n"
	      ,dmy, (2*TRN_X-1)*(2*TRN_Y-1)*(2*TRN_Z-1)
	      ,(double)dmy/( (2*TRN_X-1)*(2*TRN_Y-1)*(2*TRN_Z-1))
	      ,(double)dmy/( NX*NY*NZ));
    }
  }
  if(kBT > 0.0){
    //const double truncate_factor= 1.0;
    double truncate_factor= NX*NY*NZ/(double)dmy;
    //const double truncate_factor= (double)NX/(2.*TRN_X-1.)*(double)NY/(2.*TRN_Y-1.)*(double)NZ/(2.*TRN_Z-1.);
    //Add_random_stress_NS(zeta, f_ns0, uk_dc, ijk_range, n_ijk_range, jikan, truncate_factor, p);
    //Add_random_stress_NS0(zeta, f_ns0, uk_dc, ijk_range, n_ijk_range, jikan, truncate_factor, p);
    //Add_random_stress_NSm1(zeta, f_ns0, uk_dc, ijk_range, n_ijk_range, jikan, truncate_factor, p);
    Add_random_stress_NS1(zeta, f_ns0, uk_dc, ijk_range, n_ijk_range, jikan, truncate_factor, p);
    //Add_random_stress_NS4(zeta, f_ns0, uk_dc, ijk_range, n_ijk_range, jikan, truncate_factor, p);
    //Add_random_stress_NS2(zeta, f_ns0, uk_dc, ijk_range, n_ijk_range, jikan, truncate_factor, p);
    //Add_random_stress_NS3(zeta, f_ns0, uk_dc, ijk_range, n_ijk_range, jikan, truncate_factor, p);
  }
}
void NS_solver_Euler1(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p){
  int dmy=0;
  if(1){
    for(int n=0;n<n_ijk_range;n++){
      for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	  for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	    dmy++;
	  }
	}
      }
    }
    if(0){
      fprintf(stderr, "%d %d %g %g\n"
	      ,dmy, (2*TRN_X-1)*(2*TRN_Y-1)*(2*TRN_Z-1)
	      ,(double)dmy/( (2*TRN_X-1)*(2*TRN_Y-1)*(2*TRN_Z-1))
	      ,(double)dmy/( NX*NY*NZ));
    }
  }
  if(kBT > 0.0){
    double truncate_factor= NX*NY*NZ/(double)dmy;
    Add_random_stress_NS1(zeta, f_ns0, uk_dc, ijk_range, n_ijk_range, jikan, truncate_factor, p);
  }
  Rhs_NS(zeta, uk_dc, f_ns0, ijk_range, n_ijk_range, p);
  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Euler(DIM-1, zeta, jikan, f_ns0, ijk_range[n]);
  }
}
void NS_solver_Euler2(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p){
  Rhs_NS(zeta, uk_dc, f_ns0, ijk_range, n_ijk_range, p);
  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Euler(DIM-1, zeta, jikan, f_ns0, ijk_range[n]);
  }
}
void NS_solver_Euler3(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p){
  if(!STOKES){
    Euler_solver_Euler(zeta, jikan, uk_dc, ijk_range, n_ijk_range, p);
  }
  for(int d=0;d<DIM-1;d++){
    Reset_phi(f_ns0[d]);
  }
  for(int n=0;n<n_ijk_range;n++){
    Add_zeta_viscous_term(zeta, f_ns0, ijk_range[n]);
  }
  if(kBT > 0.0){
    const double truncate_factor= (double)NX/(2.*TRN_X-1.)*(double)NY/(2.*TRN_Y-1.)*(double)NZ/(2.*TRN_Z-1.);
    Add_random_stress_term_NS(zeta, f_ns0, uk_dc, ijk_range, n_ijk_range, jikan, truncate_factor, p);
  }
  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Euler(DIM-1, zeta, jikan, f_ns0, ijk_range[n]);
  }
}

void NS_solver_Heun(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p){
  Rhs_NS(zeta, uk_dc, f_ns0, ijk_range, n_ijk_range, p);
  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Euler(DIM-1, zeta, jikan, f_ns0, ijk_range[n]);
  }

  Rhs_NS(zeta, uk_dc, f_particle, ijk_range, n_ijk_range, p);
  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Heun(DIM-1, zeta, jikan, f_ns0, f_particle, ijk_range[n]);
  }
}
void NS_solver_slavedEuler(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p){
  if(STOKES){
    SP_stokes_advection(zeta, uk_dc, f_ns0, ijk_range, n_ijk_range, p);
  }else {
    SP_full_advection(zeta, uk_dc, f_ns0, ijk_range, n_ijk_range, p);
    //Zeta_k2advection_k(zeta, uk_dc, f_ns0);
  }

  double dmy0 = -NU*jikan.dt_fluid; 
  const double dmy1 = 1./NU;
  for(int n=0;n<n_ijk_range;n++){
    for(int d=0; d<DIM-1; d++){
      for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	  for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	    //double dmy = expm1(dmy0*K2[i][j][k]); 
	    double dmy = exp(dmy0*K2[i][j][k])-1.; 
	    zeta[d][i][j][k] += 
	      dmy * (zeta[d][i][j][k] 
		     - dmy1 * IK2[i][j][k] * f_ns0[d][i][j][k]);
	  }
	}
      }
    }
  }
}
void NS_solver_slavedHeun(Value zeta[DIM-1]
			  ,const CTime &jikan
			  ,double uk_dc[DIM]
			  ,const Index_range *ijk_range
			  ,const int &n_ijk_range
			  ,Particle *p
			  ){
  if(STOKES){
    SP_stokes_advection(zeta, uk_dc, f_ns0, ijk_range, n_ijk_range, p);
  }else {
    SP_full_advection(zeta, uk_dc, f_ns0, ijk_range, n_ijk_range, p);
    //Zeta_k2advection_k(zeta, uk_dc, rhs);
  }

  double dmy0 = -NU*jikan.dt_fluid; 
  double dmy1 = 1./NU;
  double dmy2 = 1./jikan.dt_fluid;
  for(int n=0;n<n_ijk_range;n++){
    for(int d=0; d<DIM-1; d++){
      for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	  for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	    double dmy00 = exp(dmy0*K2[i][j][k])-1.; //exp(-nu k^{2} h ) -1
	    double dmy01 = dmy1 * IK2[i][j][k]; // 1/( nu k^2)
	    double dmy02 = dmy01 * dmy2;// 1/( nu k^2 h)
	    double dmy03 = dmy01 * f_ns0[d][i][j][k];
	    zeta[d][i][j][k] += 
	      (
	       dmy00 * (zeta[d][i][j][k] - (1+dmy02) * dmy03)
	       -dmy03
	       );
	  }
	}
      }
    }
  }
  if(STOKES){
    SP_stokes_advection(zeta, uk_dc, f_ns0, ijk_range, n_ijk_range, p);
  }else {
    SP_full_advection(zeta, uk_dc, f_ns0, ijk_range, n_ijk_range, p);
    //Zeta_k2advection_k(zeta, uk_dc, rhs);
  }
  for(int n=0;n<n_ijk_range;n++){
    for(int d=0; d<DIM-1; d++){
      for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	  for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	    double dmy00 = exp(dmy0*K2[i][j][k])-1.; //exp(-nu k^{2} h ) -1
	    double dmy01 = dmy1 * IK2[i][j][k];// 1/( nu k^2)
	    double dmy02 = dmy01 * dmy2;// 1/( nu k^2 h)
	    zeta[d][i][j][k] += 
	      (dmy02 *(jikan.dt_fluid  + dmy00 * dmy01) * f_ns0[d][i][j][k]);
	  }
	}
      }
    }
  }
}

inline void FP_force_field(Value rhs[DIM-1], double force_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p){
  {/////////////////////// force field construction
    Reset_particle_force_current(p);
    Force(p);
    Reset_phi_u(phi, up);
    Make_phi_force(phi, up, p);
  }
  {
    // rotation of force field 
    double dmy[DIM]; 
    U2zeta_k(u, dmy, up);
    for(int n=0;n<n_ijk_range;n++){
      for(int d=0;d<DIM-1;d++){
	for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	  for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	    for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	      rhs[d][i][j][k] += u[d][i][j][k];
	    }
	  }
	}
      }
    }
    for(int d=0;d<DIM-1;d++){
      force_dc[d] += dmy[d];
    }
  }
}
inline void FPD_force_field(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p){

  for(int d=0;d<DIM-1;d++){
    Reset_phi(f_ns0[d]);
  }
  static double force_dc[DIM];
  for(int d=0;d<DIM;d++){ 
    force_dc[d] = 0.0;
  } 
  FP_force_field(f_ns0, force_dc, ijk_range, n_ijk_range, p);

  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Euler(DIM-1, zeta, jikan, f_ns0, ijk_range[n]);
  }
  if(1){
    for(int d=0;d<DIM;d++){
      uk_dc[d] += jikan.dt_fluid * force_dc[d];
    } 
  }
}
void Slippy_NS_solver_Euler(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p){
  
  {// advection
    for(int d=0;d<DIM-1;d++){
      Reset_phi(f_ns0[d]);
    }
    Advection_slippy_NS(zeta,uk_dc, f_ns0, ijk_range, n_ijk_range, p);
    for(int n=0;n<n_ijk_range;n++){
      Field_solver_Euler(DIM-1, zeta, jikan, f_ns0, ijk_range[n]);
    }
    if(jikan.ts == 0){
      MD_solver_position_Euler(p, jikan);
    }else{
      MD_solver_position_Euler(p, jikan);
      //MD_solver_position_AB2(p, jikan);
    }
  }
  double force_dc[DIM]={0.,0.,0.};
  {// non-advection
    for(int d=0;d<DIM-1;d++){
      Reset_phi(f_ns0[d]);
    }
    if(0){
      Rhs_slippy_NS_viscosity(zeta, f_ns0, force_dc, ijk_range, n_ijk_range, p);
      Rhs_slippy_NS_surface_normal(zeta, force_dc, f_ns0, ijk_range, n_ijk_range, p);
    }else{
      Rhs_slippy_NS_stress(zeta, uk_dc, f_ns0, ijk_range, n_ijk_range, p);
    }
    FP_force_field(f_ns0, force_dc, ijk_range, n_ijk_range, p);
    if(1){
      for(int d=0;d<DIM;d++){
	uk_dc[d] += jikan.dt_fluid * force_dc[d];
      } 
    }
    for(int n=0;n<n_ijk_range;n++){
      Field_solver_Euler(DIM-1, zeta, jikan, f_ns0, ijk_range[n]);
    }
  }
}

void QS_solver_Euler(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range){
  QS_rhs(zeta, uk_dc, Tensor_order, Stress, Q_rhs, ijk_range, n_ijk_range);
  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Euler(DIM-1, zeta, jikan, Stress, ijk_range[n]);
    Field_solver_Euler(QDIM, Tensor_order, jikan, Q_rhs, ijk_range[n]);
  }
}

inline double I_linear_gain(const Value q_k[QDIM], const Value q_rhs[QDIM]){
  static double h_control_coeff = .5; // must be <1
  static double h,dmy_rhs, dmy_q;
  h = DBL_MAX;
  for(int d=0;d<QDIM;d++){
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ_; k++){
	  dmy_q = q_k[d][i][j][k];
	  dmy_rhs = q_rhs[d][i][j][k];
	  if(dmy_rhs != 0.0 && dmy_q !=0.0){
	    dmy_rhs = ABS(dmy_q/dmy_rhs);
	    h = MIN(dmy_rhs, h);
	  }
	}
      }
    }
  }
  return h * h_control_coeff;
}
void AC_solver_Euler(CTime &jikan
		     ,Particle *p
		     ,const Index_range *ijk_range
		     ,const int &n_ijk_range
		     ){
  //AC_rhs_test(Tensor_order, Q_rhs);
  AC_rhs(Tensor_order, Q_rhs, p);
  if(0){
    double h = I_linear_gain(Tensor_order, Q_rhs);
    fprintf(stdout, "#time: %g %g\n",h ,jikan.dt_fluid);
    if(0){
      jikan.dt_md = jikan.dt_fluid = h;
      jikan.hdt_md = jikan.hdt_fluid = h * .5;
    }
  }
  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Euler(QDIM, Tensor_order, jikan, Q_rhs, ijk_range[n]);
  }
}
void OG_solver_Euler(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range){
  OG_rhs(zeta, uk_dc, Tensor_order, Stress, Q_rhs, ijk_range, n_ijk_range);
  for(int n=0;n<n_ijk_range;n++){
    //Field_solver_Euler(DIM-1, zeta, jikan, Stress, ijk_range[n]);
    Field_solver_Euler(QDIM, Tensor_order, jikan, Q_rhs, ijk_range[n]);
  }
}


inline void Rhs_solvent(Value u_solvent[DIM]
			,Value rhs_solvent[DIM-1]
			 ,Value grad_potential[DIM]
			 ,Value &source_scalar
			 ,double rhs_uk_dc[DIM]
			 ){
  U2advection_k(u_solvent, rhs_solvent); 
  if(NS_source){
    for(int d=0; d<DIM; d++){
      for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	  for(int k=0; k<NZ; k++){
	    grad_potential[d][i][j][k] *= source_scalar[i][j][k];
	  }
	}
      }	
    }
    U2zeta_k(u_solvent, rhs_uk_dc, grad_potential);
    for(int d=0; d<DIM-1; d++){
      for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	  for(int k=0; k<NZ_; k++){
	    rhs_solvent[d][i][j][k] += u_solvent[d][i][j][k];
	  }
	}
      }	
    }
  }else{
    for(int d=0; d<DIM; d++){
      rhs_uk_dc[d] = 0.0;
    }
  }
}
inline void Rhs_NS_solute(Particle *p
			   ,Value zeta_k[DIM-1]
			   ,const double uk_dc[DIM]
			   ,Value u[DIM] // working memory
			   ,Value *concentration_k
			   ,Value rhs_ns[DIM-1]
			   ,Value *rhs_solute
			   ,const Index_range *ijk_range
			   ,const int &n_ijk_range
			   ,Value solute_flux[DIM] // working memory
			   ,Value grad_potential[DIM] 
			   ,Value &omega_rhs_scalar // working memory
			   ,Value surface_normal[DIM] // working memory
			   ,double rhs_uk_dc[DIM] 
			   ){
  Truncate_vector_two_third_rule(zeta_k, DIM-1);
  Zeta_k2u(zeta_k, uk_dc, u);
  Make_surface_normal(surface_normal, p);
  Reset_phi(omega_rhs_scalar);
  for(int n=0; n<N_spec; n++){
    Truncate_two_third_rule(concentration_k[n]);
    Diffusion_flux_single(solute_flux,concentration_k[n],Onsager_coeff[n],rhs_ns[0]);
    A_k2a_out(concentration_k[n], rhs_solute[n]);
    Solute_solver_rhs_nonlinear_x_single(grad_potential
					 ,rhs_solute[n]
					 ,solute_flux
					 ,Valency_e[n]
					 ,Onsager_coeff[n]
					 ,omega_rhs_scalar
					 );
    Solute_impermeability(p, solute_flux, surface_normal);
    Add_advection_flux(solute_flux, u, rhs_solute[n]);
    U2u_k(solute_flux);
    U_k2divergence_k(solute_flux, rhs_solute[n]);
  }
  Rhs_solvent(u, rhs_ns, grad_potential, omega_rhs_scalar, rhs_uk_dc);
  for(int n=0;n<n_ijk_range;n++){
    Add_zeta_viscous_term(zeta_k, rhs_ns, ijk_range[n]);
  }
}
inline void Add_constant_field_k(Value grad_potential_k[DIM]
				 ,double e_ext[DIM]
				 ,const CTime &jikan
				 ){
  static const double nxnynz = (double)(NX*NY*NZ);
  static double amp = nxnynz;
  if(AC){
    amp *= sin( Angular_Frequency * jikan.time);
  }

  for(int d=0;d<DIM;d++){
    grad_potential_k[d][0][0][0] -= e_ext[d] * amp;
  }
}

inline void Rhs_NS_two_fluid(Particle *p
			      ,Value zeta_k[DIM-1]
			      ,const double uk_dc[DIM]
			      ,Value u[DIM] // working memory
			      ,Value *concentration_k
			      ,Value rhs_ns[DIM-1]
			      ,Value *rhs_solute
			      ,const Index_range *ijk_range
			      ,const int &n_ijk_range
			      ,Value solute_flux[DIM] // working memory
			      ,Value grad_potential[DIM] // working memory
			      ,Value surface_normal[DIM] // working memory
			      ,Value &omega_rhs_scalar // working memory
			      ,double rhs_uk_dc[DIM] 
			      ){
  
  { // potential gradient in x-space
    for(int d=0;d<DIM;d++){
      Reset_phi(grad_potential[d], -E_ext[d]);
    }
  }
  Rhs_NS_solute(p, zeta_k, uk_dc, u
		 ,concentration_k, rhs_ns, rhs_solute
		 ,ijk_range, n_ijk_range
		 ,solute_flux, grad_potential
		 ,omega_rhs_scalar,surface_normal
		 ,rhs_uk_dc
		 );
}
inline void Rhs_NS_Nernst_Planck(Particle *p
				  ,Value zeta_k[DIM-1]
				  ,const double uk_dc[DIM]
				  ,Value u[DIM] // working memory
				  ,Value *concentration_k
				  ,Value rhs_ns[DIM-1]
				  ,Value *rhs_solute
				  ,const Index_range *ijk_range
				  ,const int &n_ijk_range
				  ,Value solute_flux[DIM] // working memory
				  ,Value grad_potential[DIM] // working memory
				  ,Value surface_normal[DIM] // working memory
				  ,Value &omega_rhs_scalar // working memory
				  ,double rhs_uk_dc[DIM] 
				  ,const CTime &jikan
				  ){
  { // potential gradient in x-space
    Conc_k2charge_field(p, concentration_k, u[0], u[1], u[2]);
    A2a_k(u[0]);
    Charge_field_k2Coulomb_potential_k_PBC(u[0]);
    Truncate_two_third_rule(u[0]);
    A_k2da_k(u[0], grad_potential);
    Add_constant_field_k(grad_potential, E_ext, jikan);
    U_k2u(grad_potential);
  }
  Rhs_NS_solute(p, zeta_k, uk_dc, u
		 ,concentration_k, rhs_ns, rhs_solute
		 ,ijk_range, n_ijk_range
		 ,solute_flux, grad_potential
		 ,omega_rhs_scalar,surface_normal
		 , rhs_uk_dc
		 );
}

void Ion_diffusion_solver_Euler(Value zeta[DIM-1]
				,const CTime &jikan
				,double uk_dc[DIM]
				,Value *concentration_k
				,Particle *p
				,const Index_range *ijk_range
				,const int &n_ijk_range
				){
  double dc_rhs[DIM]; 
  {
    Rhs_NS_Nernst_Planck(p, zeta, uk_dc, u, concentration_k, f_ns0, Concentration_rhs0, ijk_range, n_ijk_range, up, f_particle, Surface_normal, phi, dc_rhs, jikan);
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
void NSsolute_solver_Euler(Value zeta[DIM-1]
			   ,const CTime &jikan
			   ,double uk_dc[DIM]
			   ,Value *concentration_k
			   ,Particle *p
			   ,const Index_range *ijk_range
			   ,const int &n_ijk_range
			   ){
  double dc_rhs[DIM]; 
  if(SW_EQ == Two_fluid){
    Rhs_NS_two_fluid(p, zeta, uk_dc, u, concentration_k, f_ns0, Concentration_rhs0, ijk_range, n_ijk_range, up, f_particle, Surface_normal, phi, dc_rhs);
  }else if(SW_EQ == Electrolyte){
    Rhs_NS_Nernst_Planck(p, zeta, uk_dc, u, concentration_k, f_ns0, Concentration_rhs0, ijk_range, n_ijk_range, up, f_particle, Surface_normal, phi, dc_rhs, jikan);
  }
  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Euler(DIM-1, zeta, jikan, f_ns0, ijk_range[n]);
    Field_solver_Euler(N_spec, concentration_k, jikan, Concentration_rhs0, ijk_range[n]);
  }
  for(int d=0;d<DIM;d++){
    uk_dc[d] += jikan.dt_fluid *  dc_rhs[d];
  }
}

void NSsolute_solver_Heun(Value zeta[DIM-1]
			  ,const CTime &jikan
			  ,double uk_dc[DIM]
			  ,Value *concentration_k
			  ,Particle *p
			  ,const Index_range *ijk_range
			  ,const int &n_ijk_range
			  ){
  double dc_rhs0[DIM]; 
  if(SW_EQ == Two_fluid){
    Rhs_NS_two_fluid(p, zeta, uk_dc, u, concentration_k, f_ns0, Concentration_rhs0, ijk_range, n_ijk_range, up, f_particle, Surface_normal, phi, dc_rhs0);
  }else if(SW_EQ == Electrolyte){
    Rhs_NS_Nernst_Planck(p, zeta, uk_dc, u, concentration_k, f_ns0, Concentration_rhs0, ijk_range, n_ijk_range, up, f_particle, Surface_normal, phi, dc_rhs0, jikan);
  }
  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Euler(DIM-1, zeta, jikan, f_ns0,ijk_range[n]);
    Field_solver_Euler(N_spec, concentration_k, jikan, Concentration_rhs0,ijk_range[n]);
  }
  for(int d=0;d<DIM;d++){
    uk_dc[d] += jikan.dt_fluid *  dc_rhs0[d];
  }

  double dc_rhs1[DIM]; 
  if(SW_EQ == Two_fluid){
    Rhs_NS_two_fluid(p, zeta, uk_dc, u, concentration_k, f_ns1, Concentration_rhs1,ijk_range, n_ijk_range, up, f_particle, Surface_normal, phi, dc_rhs1);
    //Rhs_NS_two_fluid(p, zeta, uk_dc, u, concentration_k, f_ns1, Concentration_rhs1,ijk_range, n_ijk_range, up, f_particle, Surface_normal, scalar, dc_rhs1);
  }else if(SW_EQ == Electrolyte){
    Rhs_NS_Nernst_Planck(p, zeta, uk_dc, u, concentration_k, f_ns1, Concentration_rhs1,ijk_range, n_ijk_range, up, f_particle, Surface_normal, phi, dc_rhs1, jikan);
    //Rhs_NS_Nernst_Planck(p, zeta, uk_dc, u, concentration_k, f_ns1, Concentration_rhs1,ijk_range, n_ijk_range, up, f_particle, Surface_normal, scalar, dc_rhs1, jikan);
  }
  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Heun(DIM-1, zeta, jikan, f_ns0, f_ns1, ijk_range[n]);
    Field_solver_Heun(N_spec, concentration_k, jikan, Concentration_rhs0, Concentration_rhs1, ijk_range[n]);
  }
  for(int d=0;d<DIM;d++){
    uk_dc[d] += jikan.hdt_fluid *  (dc_rhs1[d] - dc_rhs0[d]);
  }
}


void Shear_NS_solver_Euler(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p, Value force[DIM]){
  uk_dc[1] = 0.0;
  Symmetrize_zetak(zeta);
  Rhs_NS(zeta, uk_dc, f_ns0, ijk_range, n_ijk_range, p);
  Symmetrize_zetak(f_ns0);
  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Euler(DIM-1, zeta, jikan, f_ns0, ijk_range[n]);
  }
  Mean_shear_sustaining_force(zeta, uk_dc, force);
}
void Shear_NS_solver_Heun(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p, Value force[DIM]){
  uk_dc[1] = 0.0;
  Symmetrize_zetak(zeta);
  Rhs_NS(zeta, uk_dc, f_ns0, ijk_range, n_ijk_range, p);
  Symmetrize_zetak(f_ns0);
  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Euler(DIM-1, zeta, jikan, f_ns0, ijk_range[n]);
  }
  
  uk_dc[1] = 0.0;
  Symmetrize_zetak(zeta);
  Rhs_NS(zeta, uk_dc, f_particle, ijk_range, n_ijk_range, p);
  Symmetrize_zetak(f_particle);
  for(int n=0;n<n_ijk_range;n++){
    Field_solver_Heun(DIM-1, zeta, jikan, f_ns0, f_particle, ijk_range[n]);
  }
  Mean_shear_sustaining_force(zeta, uk_dc, force);
}
