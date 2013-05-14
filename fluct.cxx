//
// $Id: fluct.cxx,v 1.16 2005/10/08 13:10:49 nakayama Exp $
//
#include "fluct.h"


void Add_random_force(Particle *p, const CTime &jikan){
  static const double Zeta_drag = 6.* M_PI * ETA * RADIUS;
  static const double Zeta_drag_rot = 8.* M_PI * ETA * POW3(RADIUS);

  const double sdv_v = sqrt(Zeta_drag * jikan.dt_md * kBT);
  const double sdv_omega = sqrt(Zeta_drag_rot * jikan.dt_md * kBT);

  for(int n=0; n<Particle_Number; n++){
    double dmy[6];
    Gauss2(dmy);
    Gauss2(dmy+2);
    Gauss2(dmy+4);

    double imass = IMASS[p[n].spec];
    double imoi = IMOI[p[n].spec];
    for(int d=0; d<DIM; d++){  
      p[n].v[d] += sdv_v * dmy[d] * imass;
      p[n].omega[d] += sdv_omega * dmy[d+3] * imoi;
    }

  }
}
void BD_solver_position_Euler(Particle *p, const CTime &jikan){
  static const double iZeta_drag = 1./(6.* M_PI * ETA * RADIUS);
  static const double sdv_coeff = sqrt(2.* kBT*iZeta_drag);
  static double sdv, idt;
  idt = 1./jikan.dt_md;
  sdv = sdv_coeff * sqrt(jikan.dt_md);
  Force(p);
  for(int n=0; n<Particle_Number; n++){
    double dmy[4];
    Gauss2(dmy);
    Gauss2(dmy+2);
    for(int d=0; d<DIM; d++){  
      p[n].v[d] = (jikan.dt_md * iZeta_drag * 
		   (p[n].fr[d] + p[n].fv[d] )
		   + sdv * dmy[d]);
      p[n].x[d] += p[n].v[d];
      p[n].v[d] *= idt;
      p[n].x[d] = fmod(p[n].x[d]+L[d] , L[d]);

      p[n].fr_previous[d] = p[n].fr[d];
      p[n].fr[d] = 0.0;
      p[n].fv_previous[d] = p[n].fv[d];
      p[n].fv[d] = 0.0;

    }
  }
}

void BD_solver_momentum(Particle *p, const CTime &jikan){
  static const double Zeta_drag = 6.* M_PI * ETA * RADIUS;
  static const double iZeta_drag = 1./Zeta_drag;
  static double idt;
  idt = 1./jikan.dt_md;
  Force(p);
  for(int n=0; n<Particle_Number; n++){
    double dmy0 = IMASS[p[n].spec] * Zeta_drag * jikan.dt_md;
    //double dmy1= -expm1(-dmy0);
    double dmy1= 1.-exp(-dmy0);
    double dmy2 = MASS[p[n].spec] * iZeta_drag;
    double dmy3 =  dmy1 * dmy2;
    double dmy4 = 1.-dmy1;
    double dmy7 = iZeta_drag*(jikan.dt_md - dmy3);
    double dmy8 = iZeta_drag*(dmy1 -1. + dmy3*idt);
    
    double dmy5 = kBT * IMASS[p[n].spec];
    double sdv_r = dmy2 * sqrt(dmy5 * (2.*dmy0-dmy1*(2.+dmy1)));
    double sdv_v = sqrt(dmy5*dmy1*(2.-dmy1));
    double corr_rv = dmy5 * dmy2*SQ(dmy1)/(sdv_r*sdv_v);
    double dmy6  = sqrt(1.-SQ(corr_rv));
    for(int d=0; d<DIM; d++){  
      double dmy[2];
      Gauss2(dmy);
      p[n].x[d] += (
		    p[n].v[d] * dmy3 
		    + (p[n].fr[d] + p[n].fv[d]) * dmy7
		    + sdv_r * (dmy6 * dmy[0] + corr_rv * dmy[1])
		    );
      p[n].x[d] = fmod(p[n].x[d]+L[d] , L[d]);
      
      p[n].v[d] = p[n].v[d] * dmy4
	+ ( p[n].fr[d] + p[n].fv[d] ) * dmy8
	+sdv_v * dmy[1];
      
      p[n].fr_previous[d] = p[n].fr[d];
      p[n].fr[d] = 0.0;
      p[n].fv_previous[d] = p[n].fv[d];
      p[n].fv[d] = 0.0;
    } 
  }
  Force(p);
  for(int n=0; n<Particle_Number; n++){
    double dmy0 = IMASS[p[n].spec] * Zeta_drag * jikan.dt_md;
    //double dmy1= -expm1(-dmy0);
    double dmy1= 1.-exp(-dmy0);
    double dmy2 = MASS[p[n].spec] * iZeta_drag;
    double dmy3 =  iZeta_drag * (1. - dmy1 * dmy2 * idt) ;
    
    for(int d=0; d<DIM; d++){  
      p[n].v[d] += (( p[n].fr[d] + p[n].fv[d] )* dmy3);

      p[n].fr_previous[d] = p[n].fr[d];
      p[n].fr[d] = 0.0;
      p[n].fv_previous[d] = p[n].fv[d];
      p[n].fv[d] = 0.0;
    }
  }
}

void BD_solver_velocity_hydro(Particle *p, const CTime &jikan){
  static const double Zeta_drag = 6.* M_PI * ETA * RADIUS;
  static const double iZeta_drag = 1./Zeta_drag;
  static const double Zeta_drag_rot = 8.* M_PI * ETA * POW3(RADIUS);
  static const double iZeta_drag_rot = 1./Zeta_drag_rot;
  static double idt;
  idt = 1./jikan.dt_md;
  Force(p);
  for(int n=0; n<Particle_Number; n++){
    double dmy0 = IMASS[p[n].spec] * Zeta_drag * jikan.dt_md;
    double dmy1= 1.-exp(-dmy0);
    double dmy2 = MASS[p[n].spec] * SQ(iZeta_drag) * idt;

    double dmy_rot0 = IMOI[p[n].spec] * Zeta_drag_rot * jikan.dt_md;
    double dmy_rot1= 1.-exp(-dmy_rot0);
    double dmy_rot2 = MOI[p[n].spec] * SQ(iZeta_drag_rot) * idt;
    
    double dmy5 = kBT * IMASS[p[n].spec];
    double sdv_v = sqrt(dmy5*dmy1*(2.-dmy1));

    double dmy_rot5 = kBT * IMOI[p[n].spec];
    double sdv_omega = sqrt(dmy_rot5*dmy_rot1*(2.-dmy_rot1));

    double dmy_noise[DIM*2];
    Gauss2(dmy_noise);
    Gauss2(dmy_noise+2);
    Gauss2(dmy_noise+4);
    for(int d=0; d<DIM; d++){
      {
	p[n].v[d] += (
		      dmy1 * (-p[n].v[d]
			      + ( p[n].f_hydro[d]
				 + p[n].fr_previous[d] 
				 + p[n].fv[d] ) * iZeta_drag
			      - ( p[n].fr[d] -p[n].fr_previous[d]) * dmy2
			      )
		      + ( p[n].fr[d] -p[n].fr_previous[d]) * iZeta_drag
		      + sdv_v * dmy_noise[d]
		      );
	
	p[n].fr_previous[d] = p[n].fr[d];
	p[n].fr[d] = 0.0;
	p[n].fv_previous[d] = p[n].fv[d];
	p[n].fv[d] = 0.0;
	//p[n].f_hydro_previous[d] = p[n].f_hydro[d];
	p[n].f_hydro[d] = 0.0;
      } 
      {
	p[n].omega[d] += (
			  dmy_rot1 * (-p[n].omega[d]
				      + ( p[n].torque_hydro[d] 
					  + p[n].torquer_previous[d] 
					  + p[n].torquev[d] ) * iZeta_drag_rot
				      - ( p[n].torquer[d] 
					  - p[n].torquer_previous[d]) * dmy_rot2
				      )
			  + ( p[n].torquer[d] -p[n].torquer_previous[d]) * iZeta_drag_rot
			  + sdv_omega * dmy_noise[d]
			  );
	
	p[n].torquer_previous[d] = p[n].torquer[d];
	p[n].torquer[d] = 0.0;
	p[n].torquev_previous[d] = p[n].torquev[d];
	p[n].torquev[d] = 0.0;
	//p[n].torque_hydro_previous[d] = p[n].torque_hydro[d];
	p[n].torque_hydro[d] = 0.0;
      } 
    }
  }
}

void Add_random_stress_NSm1(Value zeta[DIM-1], Value rhs[DIM-1], double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, const CTime &jikan, const double &truncate_factor, Particle *p){
  
  Reset_phi(phi);
  Make_phi_u_particle(phi, u, p);
  {
    A2a_k(phi);
    Truncate_two_third_rule(phi);
    A_k2a(phi);
  }
  for(int d=0;d<DIM;d++){
    Reset_phi(up[d]);
  }
  for(int d=0;d<DIM;d++){
    Reset_phi(u[d]);
  }
  
  static const double sdv_coeff = sqrt( 2.*(NX*NY*NZ)/POW3(DX) * ETA * kBT);
  static double sdv;
  static double sdv_diag;
  sdv = sdv_coeff * sqrt(jikan.dt_fluid * truncate_factor);
  sdv_diag = sdv * Root_two;
  for(int n=0;n<n_ijk_range;n++){
    for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
      for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
    	for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k+=2){
	  static double random_stress[DIM][DIM][2];
	  {
	    Gauss2(random_stress[0][0]);
	    Gauss2(random_stress[0][1]);
	    Gauss2(random_stress[0][2]);
	    Gauss2(random_stress[1][1]);
	    Gauss2(random_stress[1][2]);
	    Gauss2(random_stress[2][2]);
	    
	    for(int d=0;d<2;d++){
	      random_stress[0][1][d] *= sdv;
	      random_stress[0][2][d] *= sdv;
	      random_stress[1][2][d] *= sdv;
	      
	      random_stress[0][0][d] *= sdv_diag;
	      random_stress[1][1][d] *= sdv_diag;
	      random_stress[2][2][d] *= sdv_diag;
	    }
	  }
	  {
	    for(int d=0;d<2;d++){
	      up[0][i][j][k+d]  = random_stress[0][0][d];
	      up[1][i][j][k+d]  = random_stress[0][1][d];
	      up[2][i][j][k+d]  = random_stress[0][2][d];
	      u[0][i][j][k+d]  = random_stress[1][1][d];
	      u[1][i][j][k+d]  = random_stress[1][2][d];
	      u[2][i][j][k+d]  = random_stress[2][2][d];
	    }
	  }
	}
      }
    }
  }
  {
    for(int d=0;d<DIM;d++){
      up[d][0][0][1] = 0.;
      A_k2a(up[d]);
    }
    for(int d=0;d<DIM;d++){
      u[d][0][0][1] = 0.;
      A_k2a(u[d]);
    }
  }
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ; k++){
	static double random_stress[DIM][DIM];
	static double dmy_phi;
	
	if(1){
	  dmy_phi = (1.- phi[i][j][k]);
	}else{
	  dmy_phi = 1.;
	}
	random_stress[0][0] = up[0][i][j][k] * dmy_phi; 
	random_stress[0][1] = up[1][i][j][k] * dmy_phi;
	random_stress[0][2] = up[2][i][j][k] * dmy_phi;
	random_stress[1][1] = u[0][i][j][k] * dmy_phi;
	random_stress[1][2] = u[1][i][j][k] * dmy_phi;
	random_stress[2][2] = u[2][i][j][k] * dmy_phi;
	  
	up[0][i][j][k] = random_stress[0][0]-random_stress[2][2];
	up[1][i][j][k] = random_stress[0][1];
	up[2][i][j][k] = random_stress[0][2];
	u[0][i][j][k]  = random_stress[1][1]-random_stress[2][2];
	u[1][i][j][k]  = random_stress[1][2];
      }
    }
  }
  {
    for(int d=0;d<DIM;d++){
      A2a_k(up[d]);
    }
    A2a_k(u[0]);
    A2a_k(u[1]);
  }
  
  for(int n=0;n<n_ijk_range;n++){
    for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
      for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
    	for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	  static double random_stress[DIM][DIM];
	  random_stress[0][0] = up[0][i][j][k];
	  random_stress[0][1] = up[1][i][j][k];
	  random_stress[0][2] = up[2][i][j][k];
	  random_stress[1][1] = u[0][i][j][k];
	  random_stress[1][2] = u[1][i][j][k];

	  random_stress[1][0] = random_stress[0][1];
	  random_stress[2][0] = random_stress[0][2];
	  random_stress[2][1] = random_stress[1][2];

	  random_stress[2][2] = 0.0;

	  {

	    double ks[DIM];
	    ks[0] =KX_int[i][j][k] * WAVE_X;
	    ks[1] =KY_int[i][j][k] * WAVE_Y;
	    ks[2] =KZ_int[i][j][k] * WAVE_Z;
	    double div_stress[DIM];
	    for(int d=0;d<DIM;d++){
	      div_stress[d] = 0.0;
	      for(int m=0;m<DIM;m++){
		div_stress[d] += (ks[m] * random_stress[m][d]);
	      }
	    }
	    up[0][i][j][k] =
	      -IRHO*(ks[1] * div_stress[2] - ks[2] * div_stress[1]);
	    up[1][i][j][k] =
	      -IRHO*(ks[2] * div_stress[0] - ks[0] * div_stress[2]);
	    up[2][i][j][k] =
	      -IRHO*(ks[0] * div_stress[1] - ks[1] * div_stress[0]);
	  }
	}
      }
    }
  }
  Omega_k2zeta_k(up, rhs);
  
  for(int n=0;n<n_ijk_range;n++){
    for(int d=0;d<DIM-1;d++){
      for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	  for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	    zeta[d][i][j][k] += rhs[d][i][j][k];
	  }
	}
      }
    }
  }
}

void Add_random_stress_NS(Value zeta[DIM-1], Value rhs[DIM-1], double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, const CTime &jikan, const double &truncate_factor, Particle *p){

  for(int d=0;d<DIM;d++){
    Reset_phi(up[d]);
  }

  static const double sdv_coeff = sqrt( 2.*(NX*NY*NZ)/POW3(DX) * ETA * kBT);
  static double sdv;
  static double sdv_diag;
  sdv = sdv_coeff * sqrt(jikan.dt_fluid * truncate_factor);
  sdv_diag = sdv * Root_two;
  for(int n=0;n<n_ijk_range;n++){
    for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
      for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
    	for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k+=2){
	  static double random_stress[DIM][DIM][2];
	  {
	    Gauss2(random_stress[0][0]);
	    Gauss2(random_stress[0][1]);
	    Gauss2(random_stress[0][2]);
	    Gauss2(random_stress[1][1]);
	    Gauss2(random_stress[1][2]);
	    Gauss2(random_stress[2][2]);
	    
	    for(int d=0;d<2;d++){
	      random_stress[0][1][d] *= sdv;
	      random_stress[0][2][d] *= sdv;
	      random_stress[1][2][d] *= sdv;
	      
	      random_stress[1][0][d] = random_stress[0][1][d];
	      random_stress[2][0][d] = random_stress[0][2][d];
	      random_stress[2][1][d] = random_stress[1][2][d];

	      random_stress[0][0][d] *= sdv_diag;
	      random_stress[1][1][d] *= sdv_diag;
	      random_stress[2][2][d] *= sdv_diag;
	    }
	  }
	  {
	    double ks[DIM];
	    ks[0] =KX_int[i][j][k] * WAVE_X;
	    ks[1] =KY_int[i][j][k] * WAVE_Y;
	    ks[2] =KZ_int[i][j][k] * WAVE_Z;
	    double div_stress[DIM][2];
	    for(int n=0;n<DIM;n++){
	      for(int d=0;d<2;d++){
		div_stress[n][d] = 0.0;
		for(int m=0;m<DIM;m++){
		  div_stress[n][d] 
		    += (ks[m] * random_stress[m][n][d]);
		}
	      }
	    }
	    for(int d=0;d<2;d++){
	      up[0][i][j][k+d] =
		-IRHO*(ks[1] * div_stress[2][d] - ks[2] * div_stress[1][d]);
	      up[1][i][j][k+d] =
		-IRHO*(ks[2] * div_stress[0][d] - ks[0] * div_stress[2][d]);
	      up[2][i][j][k+d] =
		-IRHO*(ks[0] * div_stress[1][d] - ks[1] * div_stress[0][d]);
	    }
	  }
	}
      }
    }
  }
  {
    for(int d=0;d<2;d++){
      up[d][0][0][1] = 0.;
    }
  }
  if(1){
    double dmy_dc[DIM]={0.,0.,0.};
    Omega_k2u_k(up, dmy_dc, u);
    U_k2u(u);
    Reset_phi(phi);
    Make_phi_u_particle(phi, up, p);
    for(int d=0; d<DIM; d++){
      for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	  for(int k=0; k<NZ; k++){
	    u[d][i][j][k] *= (1.-phi[i][j][k]);
	  }
	}
      }
    }
    U2zeta_k(rhs, dmy_dc, u);
    for(int d=0;d<DIM;d++){
      uk_dc[d] += dmy_dc[d];
    }
  }else{
    Omega_k2zeta_k(up, rhs);
  }
  for(int n=0;n<n_ijk_range;n++){
    for(int d=0;d<DIM-1;d++){
      for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	  for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	    zeta[d][i][j][k] += rhs[d][i][j][k];
	  }
	}
      }
    }
  }
}


void Add_random_stress_NS0(Value zeta[DIM-1], Value rhs[DIM-1], double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, const CTime &jikan, const double &truncate_factor, Particle *p){

  for(int d=0;d<DIM;d++){
    Reset_phi(up[d]);
  }

  static const double sdv_coeff = sqrt( 2.*(NX*NY*NZ)/POW3(DX) * ETA * kBT);
  static double sdv;
  static double sdv_diag;
  sdv = sdv_coeff * sqrt(jikan.dt_fluid * truncate_factor);
  sdv_diag = sdv * Root_two;
  for(int n=0;n<n_ijk_range;n++){
    for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
      for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
    	for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k+=2){
	  static double random_stress[DIM][DIM][2];
	  {
	    Gauss2(random_stress[0][0]);
	    Gauss2(random_stress[0][1]);
	    Gauss2(random_stress[0][2]);
	    Gauss2(random_stress[1][1]);
	    Gauss2(random_stress[1][2]);
	    Gauss2(random_stress[2][2]);
	    
	    for(int d=0;d<2;d++){
	      random_stress[0][1][d] *= sdv;
	      random_stress[0][2][d] *= sdv;
	      random_stress[1][2][d] *= sdv;
	      
	      random_stress[1][0][d] = random_stress[0][1][d];
	      random_stress[2][0][d] = random_stress[0][2][d];
	      random_stress[2][1][d] = random_stress[1][2][d];

	      random_stress[0][0][d] *= sdv_coeff;
	      random_stress[1][1][d] *= sdv_coeff;
	      random_stress[2][2][d] *= sdv_coeff;
	    }
	  }
	  {
	    double ks[DIM];
	    ks[0] =KX_int[i][j][k] * WAVE_X;
	    ks[1] =KY_int[i][j][k] * WAVE_Y;
	    ks[2] =KZ_int[i][j][k] * WAVE_Z;
	    double div_stress[DIM][2]={0.,0.
				       ,0.,0.
				       ,0.,0.};
	    for(int n=0;n<DIM;n++){
	      for(int m=0;m<DIM;m++){
		div_stress[n][0] += (ks[m] * random_stress[m][n][1]);
		div_stress[n][1] += (-ks[m] * random_stress[m][n][0]);
	      }
	    }
	    for(int n=0;n<DIM;n++){
	      for(int d=0;d<2;d++){
		up[n][i][j][k+d] = IRHO * div_stress[n][d];
	      }
	    }
	  }
	}
      }
    }
  }
  {
    for(int d=0;d<2;d++){
      up[d][0][0][1] = 0.;
    }
  }
  {
    Solenoidal_uk(up);
    U_k2u(up);
    Reset_phi(phi);
    Make_phi_u_particle(phi, u, p);
    for(int d=0; d<DIM; d++){
      for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	  for(int k=0; k<NZ; k++){
	    up[d][i][j][k] *= (1.-phi[i][j][k]);
	  }
	}
      }
    }
    double dmy_dc[DIM]={0.,0.,0.};
    U2zeta_k(rhs, dmy_dc, up);
    for(int d=0;d<DIM;d++){
      uk_dc[d] += dmy_dc[d];
    }
  }
  for(int n=0;n<n_ijk_range;n++){
    for(int d=0;d<DIM-1;d++){
      for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	  for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	    zeta[d][i][j][k] += rhs[d][i][j][k];
	  }
	}
      }
    }
  }
}
void Add_random_stress_NS1(Value zeta[DIM-1], Value rhs[DIM-1], double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, const CTime &jikan, const double &truncate_factor, Particle *p){
  
  Reset_phi(phi);
  Make_phi_u_particle(phi, u, p);
  if(0){
    A2a_k(phi);
    Truncate_two_third_rule(phi);
    A_k2a(phi);
  }
  
  static const double sdv_coeff = sqrt( 2./POW3(DX) * ETA * kBT);
  static double sdv;
  sdv = sdv_coeff * sqrt(jikan.dt_fluid * truncate_factor);
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ; k++){
	static double random_stress[DIM][DIM];
	static double dmy[2];
	static double dmy_phi;
	
	if(0){
	  dmy_phi = (1.- phi[i][j][k]) * sdv;
	}else{
	  dmy_phi = sdv;
	}
	Gauss2(dmy);
	random_stress[0][1] = dmy[0] * dmy_phi; 
	random_stress[0][2] = dmy[1] * dmy_phi; 
	Gauss2(dmy);
	random_stress[1][2] = dmy[0] * dmy_phi; 
	
	dmy_phi *= Root_two;
	random_stress[0][0] = dmy[1] * dmy_phi; 
	Gauss2(dmy);
	random_stress[1][1] = dmy[0] * dmy_phi; 
	random_stress[2][2] = dmy[1] * dmy_phi; 
	
	up[0][i][j][k] = random_stress[0][0]-random_stress[2][2];
	up[1][i][j][k] = random_stress[0][1];
	up[2][i][j][k] = random_stress[0][2];
	u[0][i][j][k]  = random_stress[1][1]-random_stress[2][2];
	u[1][i][j][k]  = random_stress[1][2];
      }
    }
  }
  {
    for(int d=0;d<DIM;d++){
      A2a_k(up[d]);
    }
    A2a_k(u[0]);
    A2a_k(u[1]);
  }
  Truncate_vector_two_third_rule(up, 3);
  Truncate_vector_two_third_rule(u, 2);
  
  for(int n=0;n<n_ijk_range;n++){
    for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
      for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
    	for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	  static double random_stress[DIM][DIM];
	  random_stress[0][0] = up[0][i][j][k];
	  random_stress[0][1] = up[1][i][j][k];
	  random_stress[0][2] = up[2][i][j][k];
	  random_stress[1][1] = u[0][i][j][k];
	  random_stress[1][2] = u[1][i][j][k];

	  random_stress[1][0] = random_stress[0][1];
	  random_stress[2][0] = random_stress[0][2];
	  random_stress[2][1] = random_stress[1][2];

	  random_stress[2][2] = 0.0;

	  {

	    double ks[DIM];
	    ks[0] =KX_int[i][j][k] * WAVE_X;
	    ks[1] =KY_int[i][j][k] * WAVE_Y;
	    ks[2] =KZ_int[i][j][k] * WAVE_Z;
	    double div_stress[DIM];
	    for(int d=0;d<DIM;d++){
	      div_stress[d] = 0.0;
	      for(int m=0;m<DIM;m++){
		div_stress[d] += (ks[m] * random_stress[m][d]);
	      }
	    }
	    up[0][i][j][k] =
	      -IRHO*(ks[1] * div_stress[2] - ks[2] * div_stress[1]);
	    up[1][i][j][k] =
	      -IRHO*(ks[2] * div_stress[0] - ks[0] * div_stress[2]);
	    up[2][i][j][k] =
	      -IRHO*(ks[0] * div_stress[1] - ks[1] * div_stress[0]);
	  }
	}
      }
    }
  }
  Omega_k2zeta_k(up, rhs);
  
  for(int n=0;n<n_ijk_range;n++){
    for(int d=0;d<DIM-1;d++){
      for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	  for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	    zeta[d][i][j][k] += rhs[d][i][j][k];
	  }
	}
      }
    }
  }
}

void Add_random_stress_NS2(Value zeta[DIM-1], Value rhs[DIM-1], double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, const CTime &jikan, const double &truncate_factor, Particle *p){

  static const double sdv_coeff = sqrt( 2./POW3(DX) * ETA * kBT);
  static double sdv;
  static double sdv_diag;
  sdv = sdv_coeff * sqrt(jikan.dt_fluid * truncate_factor);
  sdv_diag = sdv * Root_two;
  for(int d=0;d<DIM;d++){
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ; k++){
	  static double random_stress[DIM][DIM];
	  static double dmy[2];
	  
	  Gauss2(dmy);
	  random_stress[0][0] = dmy[0] * sdv_diag; 
	  random_stress[0][1] = dmy[1] * sdv; 
	  Gauss2(dmy);
	  random_stress[0][2] = dmy[0] * sdv; 
	  random_stress[1][1] = dmy[1] * sdv_diag; 
	  Gauss2(dmy);
	  random_stress[1][2] = dmy[0] * sdv; 
	  random_stress[2][2] = dmy[1] * sdv_diag; 

	  up[0][i][j][k] = random_stress[0][0];
	  up[1][i][j][k] = random_stress[0][1];
	  up[2][i][j][k] = random_stress[0][2];
	  u[0][i][j][k]  = random_stress[1][1];
	  u[1][i][j][k]  = random_stress[1][2];
	  u[2][i][j][k]  = random_stress[2][2];
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
  }
  Truncate_vector_two_third_rule(up, 3);
  Truncate_vector_two_third_rule(u, 3);
  for(int n=0;n<n_ijk_range;n++){
    for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
      for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
    	for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k+=2){
	  static double random_stress[DIM][DIM][2];
	  {
	    for(int d=0;d<2;d++){
	      random_stress[0][0][d] = up[0][i][j][k+d];
	      random_stress[0][1][d] = up[1][i][j][k+d];
	      random_stress[0][2][d] = up[2][i][j][k+d];
	      random_stress[1][1][d] = u[0][i][j][k+d];
	      random_stress[1][2][d] = u[1][i][j][k+d];
	      random_stress[2][2][d] = u[2][i][j][k+d];
	      
	      random_stress[1][0][d] = random_stress[0][1][d];
	      random_stress[2][0][d] = random_stress[0][2][d];
	      random_stress[2][1][d] = random_stress[1][2][d];
	    }
	  }
	  {
	    double ks[DIM];
	    ks[0] = KX_int[i][j][k] * WAVE_X;
	    ks[1] = KY_int[i][j][k] * WAVE_Y;
	    ks[2] = KZ_int[i][j][k] * WAVE_Z;
	    double div_stress[DIM][2];
	    for(int n=0;n<DIM;n++){
	      div_stress[n][0] = 0.0;
	      div_stress[n][1] = 0.0;
	      for(int m=0;m<DIM;m++){
		div_stress[n][0] += (ks[m] * random_stress[m][n][1]);
		div_stress[n][1] += (-ks[m] * random_stress[m][n][0]);
	      }
	    }
	    for(int n=0;n<DIM;n++){
	      for(int d=0;d<2;d++){
		up[n][i][j][k+d] = IRHO * div_stress[n][d];
	      }
	    }
	  }
	}
      }
    }
  }
  {
    for(int d=0;d<DIM;d++){
      up[d][0][0][1] = 0.;
    }
  }
  {
    Solenoidal_uk(up);
    U_k2u(up);
    Reset_phi(phi);
    Make_phi_u_particle(phi, u, p);
    for(int d=0; d<DIM; d++){
      for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	  for(int k=0; k<NZ; k++){
	    up[d][i][j][k] *= (1.-phi[i][j][k]);
	  }
	}
      }
    }
    double dmy_dc[DIM]={0.,0.,0.};
    U2zeta_k(rhs, dmy_dc, up);
    for(int d=0;d<DIM;d++){
      uk_dc[d] += dmy_dc[d];
    }
  }
  for(int n=0;n<n_ijk_range;n++){
    for(int d=0;d<DIM-1;d++){
      for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	  for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	    zeta[d][i][j][k] += rhs[d][i][j][k];
	  }
	}
      }
    }
  }
}


void Add_random_stress_NS3(Value zeta[DIM-1], Value rhs[DIM-1], double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, const CTime &jikan, const double &truncate_factor, Particle *p){

  Reset_phi(phi);
  Make_phi_u_particle(phi, u, p);
  
  static const double sdv_coeff = sqrt( 2./POW3(DX) * ETA * kBT);
  static double sdv;
  sdv = sdv_coeff * sqrt(jikan.dt_fluid * truncate_factor);
  for(int d=0;d<DIM;d++){
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ; k++){
	  static double random_stress[DIM][DIM];
	  static double dmy[2];
	  static double dmy_phi;
	  
	  if(1){
	    dmy_phi = (1.- phi[i][j][k]) * sdv;
	  }else{
	    dmy_phi = sdv;
	  }
	  Gauss2(dmy);
	  random_stress[0][1] = dmy[0] * dmy_phi; 
	  random_stress[0][2] = dmy[1] * dmy_phi; 
	  Gauss2(dmy);
	  random_stress[1][0] = dmy[0] * dmy_phi; 
	  random_stress[1][2] = dmy[1] * dmy_phi; 
	  Gauss2(dmy);
	  random_stress[2][0] = dmy[0] * dmy_phi; 
	  random_stress[2][1] = dmy[1] * dmy_phi; 

	  Gauss2(dmy);
	  random_stress[0][0] = dmy[0] * dmy_phi; 
	  random_stress[1][1] = dmy[1] * dmy_phi; 
	  Gauss2(dmy);
	  random_stress[2][2] = dmy[0] * dmy_phi; 
	  
	  up[0][i][j][k] = random_stress[0][0]-random_stress[2][2];
	  up[1][i][j][k] = random_stress[0][1];
	  up[2][i][j][k] = random_stress[0][2];
	  u[0][i][j][k]  = random_stress[1][0];
	  u[1][i][j][k]  = random_stress[1][1]-random_stress[2][2];
	  u[2][i][j][k]  = random_stress[1][2];
	  rhs[0][i][j][k]  = random_stress[2][0];
	  rhs[1][i][j][k]  = random_stress[2][1];
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
  Truncate_vector_two_third_rule(up, DIM);
  Truncate_vector_two_third_rule(u, DIM);
  Truncate_vector_two_third_rule(rhs, DIM-1);
  
  for(int n=0;n<n_ijk_range;n++){
    for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
      for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
    	for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	  static double random_stress[DIM][DIM];
	  random_stress[0][0] = up[0][i][j][k];
	  random_stress[0][1] = up[1][i][j][k];
	  random_stress[0][2] = up[2][i][j][k];
	  random_stress[1][0] = u[0][i][j][k];
	  random_stress[1][1] = u[1][i][j][k];
	  random_stress[1][2] = u[2][i][j][k];
	  random_stress[2][0] = rhs[0][i][j][k];
	  random_stress[2][1] = rhs[1][i][j][k];
	  random_stress[2][2] = 0.0;

	  {
	    double ks[DIM];
	    ks[0] =KX_int[i][j][k] * WAVE_X;
	    ks[1] =KY_int[i][j][k] * WAVE_Y;
	    ks[2] =KZ_int[i][j][k] * WAVE_Z;
	    double div_stress[DIM];
	    for(int d=0;d<DIM;d++){
	      div_stress[d] = 0.0;
	      for(int m=0;m<DIM;m++){
		div_stress[d] += (ks[m] * random_stress[m][d]);
	      }
	    }
	    up[0][i][j][k] =
	      -IRHO*(ks[1] * div_stress[2] - ks[2] * div_stress[1]);
	    up[1][i][j][k] =
	      -IRHO*(ks[2] * div_stress[0] - ks[0] * div_stress[2]);
	    up[2][i][j][k] =
	      -IRHO*(ks[0] * div_stress[1] - ks[1] * div_stress[0]);
	  }
	}
      }
    }
  }
  Omega_k2zeta_k(up, rhs);
  
  for(int n=0;n<n_ijk_range;n++){
    for(int d=0;d<DIM-1;d++){
      for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	  for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	    zeta[d][i][j][k] += rhs[d][i][j][k];
	  }
	}
      }
    }
  }
}
void Add_random_stress_NS4(Value zeta[DIM-1], Value rhs[DIM-1], double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, const CTime &jikan, const double &truncate_factor, Particle *p){

  Reset_phi(phi);
  Make_phi_u_particle(phi, u, p);
  
  static const double sdv_coeff = sqrt( 2./POW3(DX) * ETA * kBT);
  static double sdv;
  //sdv = sdv_coeff * sqrt(jikan.dt_fluid * truncate_factor);
  sdv = sdv_coeff * sqrt(jikan.dt_fluid);
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ; k++){
	static double random_stress[DIM][DIM];
	static double dmy[2];
	static double dmy_phi;
	
	if(1){
	  dmy_phi = (1.- phi[i][j][k]) * sdv;
	}else{
	  dmy_phi = sdv;
	}
	Gauss2(dmy);
	random_stress[0][1] = dmy[0] * dmy_phi; 
	random_stress[0][2] = dmy[1] * dmy_phi; 
	Gauss2(dmy);
	random_stress[1][2] = dmy[0] * dmy_phi; 
	  
	dmy_phi *= Root_two;
	random_stress[0][0] = dmy[1] * dmy_phi; 
	Gauss2(dmy);
	random_stress[1][1] = dmy[0] * dmy_phi; 
	random_stress[2][2] = dmy[1] * dmy_phi; 
	  
	up[0][i][j][k] = random_stress[0][0]-random_stress[2][2];
	up[1][i][j][k] = random_stress[0][1];
	up[2][i][j][k] = random_stress[0][2];
	u[0][i][j][k]  = random_stress[1][1]-random_stress[2][2];
	u[1][i][j][k]  = random_stress[1][2];
      }
    }
  }
  {
    for(int d=0;d<DIM;d++){
      A2a_k(up[d]);
    }
    A2a_k(u[0]);
    A2a_k(u[1]);
  }
  
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ_; k+=2){
	static double random_stress[DIM][DIM][2];
	{
	  for(int d=0;d<2;d++){
	    random_stress[0][0][d] = up[0][i][j][k+d];
	    random_stress[0][1][d] = up[1][i][j][k+d];
	    random_stress[0][2][d] = up[2][i][j][k+d];
	    random_stress[1][1][d] = u[0][i][j][k+d];
	    random_stress[1][2][d] = u[1][i][j][k+d];
	    random_stress[2][2][d] = 0.0;
	    
	    random_stress[1][0][d] = random_stress[0][1][d];
	    random_stress[2][0][d] = random_stress[0][2][d];
	    random_stress[2][1][d] = random_stress[1][2][d];
	  }
	}
	{
	  double ks[DIM];
	  ks[0] = KX_int[i][j][k] * WAVE_X;
	  ks[1] = KY_int[i][j][k] * WAVE_Y;
	  ks[2] = KZ_int[i][j][k] * WAVE_Z;
	  double div_stress[DIM][2];
	  for(int n=0;n<DIM;n++){
	    div_stress[n][0] = 0.0;
	    div_stress[n][1] = 0.0;
	    for(int m=0;m<DIM;m++){
	      div_stress[n][0] += (ks[m] * random_stress[m][n][1]);
	      div_stress[n][1] += (-ks[m] * random_stress[m][n][0]);
	    }
	  }
	  for(int n=0;n<DIM;n++){
	    for(int d=0;d<2;d++){
	      up[n][i][j][k+d] = IRHO * div_stress[n][d];
	    }
	  }
	}
      }
    }
  }
  Solenoidal_uk(up);
  U_k2u(up);

  Zeta_k2u(zeta, uk_dc, u);
  for(int d=0;d<DIM;d++){
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ; k++){
	  u[d][i][j][k] += up[d][i][j][k];
	}
      }
    }
  }
  U2zeta_k(zeta, uk_dc, u);
}
void Add_random_stress_term_NS(Value zeta[DIM-1], Value rhs[DIM-1], double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, const CTime &jikan, const double &truncate_factor, Particle *p){

  Reset_phi(phi);
  Make_phi_particle(phi, p);
  
  static const double sdv_coeff = sqrt( 2./POW3(DX) * ETA * kBT);
  static double sdv;
  sdv = sdv_coeff * sqrt(jikan.dt_fluid * truncate_factor);
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ; k++){
	static double random_stress[DIM][DIM];
	static double dmy[2];
	static double dmy_phi;
	
	if(1){
	  dmy_phi = (1.- phi[i][j][k]) * sdv;
	}else{
	  dmy_phi = sdv;
	}
	Gauss2(dmy);
	random_stress[0][1] = dmy[0] * dmy_phi; 
	random_stress[0][2] = dmy[1] * dmy_phi; 
	Gauss2(dmy);
	random_stress[1][2] = dmy[0] * dmy_phi; 
	  
	dmy_phi *= Root_two;
	random_stress[0][0] = dmy[1] * dmy_phi; 
	Gauss2(dmy);
	random_stress[1][1] = dmy[0] * dmy_phi; 
	random_stress[2][2] = dmy[1] * dmy_phi; 
	  
	up[0][i][j][k] = random_stress[0][0]-random_stress[2][2];
	up[1][i][j][k] = random_stress[0][1];
	up[2][i][j][k] = random_stress[0][2];
	u[0][i][j][k]  = random_stress[1][1]-random_stress[2][2];
	u[1][i][j][k]  = random_stress[1][2];
      }
    }
  }
  {
    for(int d=0;d<DIM;d++){
      A2a_k(up[d]);
    }
    A2a_k(u[0]);
    A2a_k(u[1]);
  }
  Truncate_vector_two_third_rule(up, 3);
  Truncate_vector_two_third_rule(u, 2);
  
  for(int n=0;n<n_ijk_range;n++){
    for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
      for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
    	for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	  static double random_stress[DIM][DIM];
	  random_stress[0][0] = up[0][i][j][k];
	  random_stress[0][1] = up[1][i][j][k];
	  random_stress[0][2] = up[2][i][j][k];
	  random_stress[1][1] = u[0][i][j][k];
	  random_stress[1][2] = u[1][i][j][k];

	  random_stress[1][0] = random_stress[0][1];
	  random_stress[2][0] = random_stress[0][2];
	  random_stress[2][1] = random_stress[1][2];

	  random_stress[2][2] = 0.0;

	  {

	    double ks[DIM];
	    ks[0] =KX_int[i][j][k] * WAVE_X;
	    ks[1] =KY_int[i][j][k] * WAVE_Y;
	    ks[2] =KZ_int[i][j][k] * WAVE_Z;
	    double div_stress[DIM];
	    for(int d=0;d<DIM;d++){
	      div_stress[d] = 0.0;
	      for(int m=0;m<DIM;m++){
		div_stress[d] += (ks[m] * random_stress[m][d]);
	      }
	    }
	    up[0][i][j][k] =
	      -IRHO*(ks[1] * div_stress[2] - ks[2] * div_stress[1]);
	    up[1][i][j][k] =
	      -IRHO*(ks[2] * div_stress[0] - ks[0] * div_stress[2]);
	    up[2][i][j][k] =
	      -IRHO*(ks[0] * div_stress[1] - ks[1] * div_stress[0]);
	  }
	}
      }
    }
  }
  Omega_k2zeta_k(up, u);
  for(int d=0; d<DIM-1; d++){
    for(int n=0;n<n_ijk_range;n++){
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
void Add_random_stress_x_NS1(Value u[DIM]
			     ,Particle *p
			     ,const Index_range *ijk_range
			     ,const int &n_ijk_range
			     ,const CTime &jikan
			     ,const double &truncate_factor
			     ,Value f_particle[DIM] // working memory
			     ,Value up[DIM] // working memory
			     ){
  Reset_phi(phi);
  Make_phi_u_particle(phi, up, p);
  
  static const double sdv_coeff = sqrt( 2./POW3(DX) * ETA * kBT);
  static double sdv;
  sdv = sdv_coeff * sqrt(jikan.dt_fluid * truncate_factor);
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ; k++){
	static double random_stress[DIM][DIM];
	static double dmy[2];
	static double dmy_phi;
	
	if(1){
	  dmy_phi = (1.- phi[i][j][k]) * sdv;
	}else{
	  dmy_phi = sdv;
	}
	Gauss2(dmy);
	random_stress[0][1] = dmy[0] * dmy_phi; 
	random_stress[0][2] = dmy[1] * dmy_phi; 
	Gauss2(dmy);
	random_stress[1][2] = dmy[0] * dmy_phi; 
	  
	dmy_phi *= Root_two;
	random_stress[0][0] = dmy[1] * dmy_phi; 
	Gauss2(dmy);
	random_stress[1][1] = dmy[0] * dmy_phi; 
	random_stress[2][2] = dmy[1] * dmy_phi; 
	  
	up[0][i][j][k] = random_stress[0][0]-random_stress[2][2];
	up[1][i][j][k] = random_stress[0][1];
	up[2][i][j][k] = random_stress[0][2];
	f_particle[0][i][j][k]  = random_stress[1][1]-random_stress[2][2];
	f_particle[1][i][j][k]  = random_stress[1][2];
      }
    }
  }
  {
    for(int d=0;d<DIM;d++){
      A2a_k(up[d]);
    }
    A2a_k(f_particle[0]);
    A2a_k(f_particle[1]);
  }
  Truncate_vector_two_third_rule(up, 3);
  Truncate_vector_two_third_rule(f_particle, 2);
  
  for(int n=0;n<n_ijk_range;n++){
    for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
      for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
    	for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	  static double random_stress[DIM][DIM];
	  random_stress[0][0] = up[0][i][j][k];
	  random_stress[0][1] = up[1][i][j][k];
	  random_stress[0][2] = up[2][i][j][k];
	  random_stress[1][1] = f_particle[0][i][j][k];
	  random_stress[1][2] = f_particle[1][i][j][k];
	  
	  random_stress[1][0] = random_stress[0][1];
	  random_stress[2][0] = random_stress[0][2];
	  random_stress[2][1] = random_stress[1][2];
	  
	  random_stress[2][2] = 0.0;
	  
	  {
	    
	    double ks[DIM];
	    ks[0] =KX_int[i][j][k] * WAVE_X;
	    ks[1] =KY_int[i][j][k] * WAVE_Y;
	    ks[2] =KZ_int[i][j][k] * WAVE_Z;
	    double div_stress[DIM];
	    for(int d=0;d<DIM;d++){
	      div_stress[d] = 0.0;
	      for(int m=0;m<DIM;m++){
		div_stress[d] += (ks[m] * random_stress[m][d]);
	      }
	    }
	    up[0][i][j][k] =
	      -IRHO*(ks[1] * div_stress[2] - ks[2] * div_stress[1]);
	    up[1][i][j][k] =
	      -IRHO*(ks[2] * div_stress[0] - ks[0] * div_stress[2]);
	    up[2][i][j][k] =
	      -IRHO*(ks[0] * div_stress[1] - ks[1] * div_stress[0]);
	  }
	}
      }
    }
  }
  //  Omega_k2zeta_k(up, rhs);
  double dmy_dc[DIM] = {0.,0.,0.};
  Omega_k2u_k(up, dmy_dc, f_particle);
  U_k2u(f_particle);
  
  for(int d=0; d<DIM; d++){
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ; k++){
	  u[d][i][j][k] += f_particle[d][i][j][k];
	}
      }
    }
  }
}
void Add_random_stress_x_NS2(Value u[DIM]
			     ,Particle *p
			     ,const Index_range *ijk_range
			     ,const int &n_ijk_range
			     ,const CTime &jikan
			     ,const double &truncate_factor
			     ,Value f_particle[DIM] // working memory
			     ,Value up[DIM] // working memory
			     ){
  Reset_phi(phi);
  Make_phi_u_particle(phi, up, p);
  
  static const double sdv_coeff = sqrt( 2./POW3(DX) * ETA * kBT);
  static double sdv;
  sdv = sdv_coeff * sqrt(jikan.dt_fluid * truncate_factor);
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ; k++){
	static double random_stress[DIM][DIM];
	static double dmy[2];
	static double dmy_phi;
	
	if(1){
	  dmy_phi = (1.- phi[i][j][k]) * sdv;
	}else{
	  dmy_phi = sdv;
	}
	Gauss2(dmy);
	random_stress[0][1] = dmy[0] * dmy_phi; 
	random_stress[0][2] = dmy[1] * dmy_phi; 
	Gauss2(dmy);
	random_stress[1][2] = dmy[0] * dmy_phi; 
	  
	dmy_phi *= Root_two;
	random_stress[0][0] = dmy[1] * dmy_phi; 
	Gauss2(dmy);
	random_stress[1][1] = dmy[0] * dmy_phi; 
	random_stress[2][2] = dmy[1] * dmy_phi; 
	  
	up[0][i][j][k] = random_stress[0][0];
	up[1][i][j][k] = random_stress[0][1];
	up[2][i][j][k] = random_stress[0][2];
	f_particle[0][i][j][k]  = random_stress[1][1];
	f_particle[1][i][j][k]  = random_stress[1][2];
	f_particle[2][i][j][k]  = random_stress[2][2];
      }
    }
  }
  {
    for(int d=0;d<DIM;d++){
      A2a_k(up[d]);
    }
    for(int d=0;d<DIM;d++){
      A2a_k(f_particle[d]);
    }
  }
  Truncate_vector_two_third_rule(up, 3);
  Truncate_vector_two_third_rule(f_particle, 3);
  
  for(int n=0;n<n_ijk_range;n++){
    for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
      for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
    	for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k+=2){
	  static double random_stress[DIM][DIM][2];
	  for(int d=0; d<2; d++){
	    random_stress[0][0][d] = up[0][i][j][k+d];
	    random_stress[0][1][d] = up[1][i][j][k+d];
	    random_stress[0][2][d] = up[2][i][j][k+d];
	    random_stress[1][1][d] = f_particle[0][i][j][k+d];
	    random_stress[1][2][d] = f_particle[1][i][j][k+d];
	    random_stress[2][2][d] = f_particle[2][i][j][k+d];
	    
	    random_stress[1][0][d] = random_stress[0][1][d];
	    random_stress[2][0][d] = random_stress[0][2][d];
	    random_stress[2][1][d] = random_stress[1][2][d];
	  }
	  
	  {
	    double ks[DIM];
	    ks[0] = KX_int[i][j][k] * WAVE_X;
	    ks[1] = KY_int[i][j][k] * WAVE_Y;
	    ks[2] = KZ_int[i][j][k] * WAVE_Z;
	    double div_stress[DIM][2];
	    for(int n=0;n<DIM;n++){
	      div_stress[n][0] = 0.0;
	      div_stress[n][1] = 0.0;
	      for(int m=0;m<DIM;m++){
		div_stress[n][0] += (ks[m] * random_stress[m][n][1]);
		div_stress[n][1] += (-ks[m] * random_stress[m][n][0]);
	      }
	    }
	    for(int n=0;n<DIM;n++){
	      for(int d=0;d<2;d++){
		up[n][i][j][k+d] = IRHO * div_stress[n][d];
	      }
	    }
	  }
	}
      }
    }
  }
  Solenoidal_uk(up);
  U_k2u(up);
  
  for(int d=0; d<DIM; d++){
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ; k++){
	  u[d][i][j][k] += up[d][i][j][k];
	}
      }
    }
  }
}


void Add_random_stress_x_NS3(Value u[DIM]
			     ,Particle *p
			     ,const Index_range *ijk_range
			     ,const int &n_ijk_range
			     ,const CTime &jikan
			     ,const double &truncate_factor
			     ,Value f_particle[DIM] // working memory
			     ,Value up[DIM] // working memory
			     ){
  static const double sdv_coeff = sqrt( 2./POW3(DX) * ETA * kBT);
  static double sdv;
  sdv = sdv_coeff * sqrt(jikan.dt_fluid * truncate_factor);
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ; k++){
	static double random_stress[DIM][DIM];
	static double dmy[2];
	static double dmy_phi;

	dmy_phi = sdv;

	Gauss2(dmy);
	random_stress[0][1] = dmy[0] * dmy_phi; 
	random_stress[0][2] = dmy[1] * dmy_phi; 
	Gauss2(dmy);
	random_stress[1][2] = dmy[0] * dmy_phi; 
	  
	dmy_phi *= Root_two;
	random_stress[0][0] = dmy[1] * dmy_phi; 
	Gauss2(dmy);
	random_stress[1][1] = dmy[0] * dmy_phi; 
	random_stress[2][2] = dmy[1] * dmy_phi; 
	  
	up[0][i][j][k] = random_stress[0][0];
	up[1][i][j][k] = random_stress[0][1];
	up[2][i][j][k] = random_stress[0][2];
	f_particle[0][i][j][k]  = random_stress[1][1];
	f_particle[1][i][j][k]  = random_stress[1][2];
	f_particle[2][i][j][k]  = random_stress[2][2];
      }
    }
  }
  {
    for(int d=0;d<DIM;d++){
      A2a_k(up[d]);
    }
    for(int d=0;d<DIM;d++){
      A2a_k(f_particle[d]);
    }
  }
  Truncate_vector_two_third_rule(up, 3);
  Truncate_vector_two_third_rule(f_particle, 3);
  
  for(int n=0;n<n_ijk_range;n++){
    for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
      for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
    	for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k+=2){
	  static double random_stress[DIM][DIM][2];
	  for(int d=0; d<2; d++){
	    random_stress[0][0][d] = up[0][i][j][k+d];
	    random_stress[0][1][d] = up[1][i][j][k+d];
	    random_stress[0][2][d] = up[2][i][j][k+d];
	    random_stress[1][1][d] = f_particle[0][i][j][k+d];
	    random_stress[1][2][d] = f_particle[1][i][j][k+d];
	    random_stress[2][2][d] = f_particle[2][i][j][k+d];
	    
	    random_stress[1][0][d] = random_stress[0][1][d];
	    random_stress[2][0][d] = random_stress[0][2][d];
	    random_stress[2][1][d] = random_stress[1][2][d];
	  }
	  
	  {
	    double ks[DIM];
	    ks[0] = KX_int[i][j][k] * WAVE_X;
	    ks[1] = KY_int[i][j][k] * WAVE_Y;
	    ks[2] = KZ_int[i][j][k] * WAVE_Z;
	    double div_stress[DIM][2];
	    for(int n=0;n<DIM;n++){
	      div_stress[n][0] = 0.0;
	      div_stress[n][1] = 0.0;
	      for(int m=0;m<DIM;m++){
		div_stress[n][0] += (ks[m] * random_stress[m][n][1]);
		div_stress[n][1] += (-ks[m] * random_stress[m][n][0]);
	      }
	    }
	    for(int n=0;n<DIM;n++){
	      for(int d=0;d<2;d++){
		up[n][i][j][k+d] = IRHO * div_stress[n][d];
	      }
	    }
	  }
	}
      }
    }
  }
  U_k2u(up);
  if(1){
    Reset_phi(phi);
    Make_phi_u_particle(phi, f_particle, p);
    for(int d=0; d<DIM; d++){
      for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	  for(int k=0; k<NZ; k++){
	    up[d][i][j][k] *= (1.-phi[i][j][k]);
	  }
	}
      }
    }
  }
  {
    U2u_k(up);
    Solenoidal_uk(up);
    for(int d=0;d<DIM;d++){
      up[d][0][0][0] = 0.0;
    }
    U_k2u(up);
  }
  for(int d=0; d<DIM; d++){
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ; k++){
	  u[d][i][j][k] += up[d][i][j][k];
	}
      }
    }
  }
}

void Add_random_stress_x_slavedNS(Value u[DIM]
			     ,Particle *p
			     ,const Index_range *ijk_range
			     ,const int &n_ijk_range
			     ,const CTime &jikan
			     ,const double &truncate_factor
			     ,Value f_particle[DIM] // working memory
			     ,Value up[DIM] // working memory
			     ){
  Reset_phi(phi);
  assert(RADIUS > XI);
  //Make_phi_particle(phi, p);
  //Make_phi_particle(phi, p, RADIUS-XI);
  Make_phi_particle(phi, p, RADIUS+XI);

  static const double sdv_coeff = sqrt( 2./POW3(DX) * ETA * kBT);
  static double sdv;
  sdv = sdv_coeff * sqrt(truncate_factor);
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ; k++){
	static double random_stress[DIM][DIM];
	static double dmy[2];
	static double dmy_phi;

	if(1){
	  dmy_phi = (1.- phi[i][j][k]) * sdv;
	}else{
	  dmy_phi = sdv;
	}
	Gauss2(dmy);
	random_stress[0][1] = dmy[0] * dmy_phi; 
	random_stress[0][2] = dmy[1] * dmy_phi; 
	Gauss2(dmy);
	random_stress[1][2] = dmy[0] * dmy_phi; 
	  
	dmy_phi *= Root_two;
	random_stress[0][0] = dmy[1] * dmy_phi; 
	Gauss2(dmy);
	random_stress[1][1] = dmy[0] * dmy_phi; 
	random_stress[2][2] = dmy[1] * dmy_phi; 
	  
	up[0][i][j][k] = random_stress[0][0];
	up[1][i][j][k] = random_stress[0][1];
	up[2][i][j][k] = random_stress[0][2];
	f_particle[0][i][j][k]  = random_stress[1][1];
	f_particle[1][i][j][k]  = random_stress[1][2];
	f_particle[2][i][j][k]  = random_stress[2][2];
      }
    }
  }
  {
    for(int d=0;d<DIM;d++){
      A2a_k(up[d]);
    }
    for(int d=0;d<DIM;d++){
      A2a_k(f_particle[d]);
    }
  }
  Truncate_vector_two_third_rule(up, 3);
  Truncate_vector_two_third_rule(f_particle, 3);

  double dmy0 = -2.*NU*jikan.dt_fluid; 
  const double dmy1 = 1/(2.*NU);
  for(int n=0;n<n_ijk_range;n++){
    for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
      for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
    	for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k+=2){
	  //double dmy = sqrt(-expm1(dmy0*K2[i][j][k])*dmy1 * IK2[i][j][k]);
	  double dmy = IRHO * sqrt((1.-exp(dmy0*K2[i][j][k]))*dmy1 * IK2[i][j][k]);
	  static double random_stress[DIM][DIM][2];
	  for(int d=0; d<2; d++){
	    random_stress[0][0][d] = up[0][i][j][k+d];
	    random_stress[0][1][d] = up[1][i][j][k+d];
	    random_stress[0][2][d] = up[2][i][j][k+d];
	    random_stress[1][1][d] = f_particle[0][i][j][k+d];
	    random_stress[1][2][d] = f_particle[1][i][j][k+d];
	    random_stress[2][2][d] = f_particle[2][i][j][k+d];
	    
	    random_stress[1][0][d] = random_stress[0][1][d];
	    random_stress[2][0][d] = random_stress[0][2][d];
	    random_stress[2][1][d] = random_stress[1][2][d];
	  }
	  
	  {
	    double ks[DIM];
	    ks[0] = KX_int[i][j][k] * WAVE_X;
	    ks[1] = KY_int[i][j][k] * WAVE_Y;
	    ks[2] = KZ_int[i][j][k] * WAVE_Z;
	    double div_stress[DIM][2];
	    for(int n=0;n<DIM;n++){
	      div_stress[n][0] = 0.0;
	      div_stress[n][1] = 0.0;
	      for(int m=0;m<DIM;m++){
		div_stress[n][0] += (ks[m] * random_stress[m][n][1]);
		div_stress[n][1] += (-ks[m] * random_stress[m][n][0]);
	      }
	    }
	    for(int n=0;n<DIM;n++){
	      for(int d=0;d<2;d++){
		up[n][i][j][k+d] = div_stress[n][d] * dmy;
	      }
	    }
	  }
	}
      }
    }
  }
  Solenoidal_uk(up);
  U_k2u(up);
  
  for(int d=0; d<DIM; d++){
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ; k++){
	  u[d][i][j][k] += up[d][i][j][k];
	}
      }
    }
  }
}
void Add_random_stress_x_up_slavedNS(Value up[DIM]
			     ,Particle *p
			     ,const Index_range *ijk_range
			     ,const int &n_ijk_range
			     ,const CTime &jikan
			     ,const double &truncate_factor
			     ,Value f_particle[DIM] // working memory
			     ){
  Reset_phi(phi);
  assert(RADIUS > XI);
  //Make_phi_particle(phi, p);
  //Make_phi_particle(phi, p, RADIUS-XI);
  Make_phi_particle(phi, p, RADIUS+XI);
  
  static const double sdv_coeff = sqrt( 2./POW3(DX) * ETA * kBT);
  static double sdv;
  sdv = sdv_coeff * sqrt(truncate_factor);
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ; k++){
	static double random_stress[DIM][DIM];
	static double dmy[2];
	static double dmy_phi;
	
	if(1){
	  dmy_phi = (1.- phi[i][j][k]) * sdv;
	}else{
	  dmy_phi = sdv;
	}
	Gauss2(dmy);
	random_stress[0][1] = dmy[0] * dmy_phi; 
	random_stress[0][2] = dmy[1] * dmy_phi; 
	Gauss2(dmy);
	random_stress[1][2] = dmy[0] * dmy_phi; 
	  
	dmy_phi *= Root_two;
	random_stress[0][0] = dmy[1] * dmy_phi; 
	Gauss2(dmy);
	random_stress[1][1] = dmy[0] * dmy_phi; 
	random_stress[2][2] = dmy[1] * dmy_phi; 
	  
	up[0][i][j][k] = random_stress[0][0];
	up[1][i][j][k] = random_stress[0][1];
	up[2][i][j][k] = random_stress[0][2];
	f_particle[0][i][j][k]  = random_stress[1][1];
	f_particle[1][i][j][k]  = random_stress[1][2];
	f_particle[2][i][j][k]  = random_stress[2][2];
      }
    }
  }
  {
    for(int d=0;d<DIM;d++){
      A2a_k(up[d]);
    }
    for(int d=0;d<DIM;d++){
      A2a_k(f_particle[d]);
    }
  }
  Truncate_vector_two_third_rule(up, 3);
  Truncate_vector_two_third_rule(f_particle, 3);

  double dmy0 = -2.*NU*jikan.dt_fluid; 
  const double dmy1 = 1/(2.*NU);
  for(int n=0;n<n_ijk_range;n++){
    for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
      for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
    	for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k+=2){
	  //double dmy = sqrt(-expm1(dmy0*K2[i][j][k])*dmy1 * IK2[i][j][k]);
	  double dmy = sqrt((1.-exp(dmy0*K2[i][j][k]))*dmy1 * IK2[i][j][k]);
	  static double random_stress[DIM][DIM][2];
	  for(int d=0; d<2; d++){
	    random_stress[0][0][d] = up[0][i][j][k+d];
	    random_stress[0][1][d] = up[1][i][j][k+d];
	    random_stress[0][2][d] = up[2][i][j][k+d];
	    random_stress[1][1][d] = f_particle[0][i][j][k+d];
	    random_stress[1][2][d] = f_particle[1][i][j][k+d];
	    random_stress[2][2][d] = f_particle[2][i][j][k+d];
	    
	    random_stress[1][0][d] = random_stress[0][1][d];
	    random_stress[2][0][d] = random_stress[0][2][d];
	    random_stress[2][1][d] = random_stress[1][2][d];
	  }
	  
	  {
	    double ks[DIM];
	    ks[0] = KX_int[i][j][k] * WAVE_X;
	    ks[1] = KY_int[i][j][k] * WAVE_Y;
	    ks[2] = KZ_int[i][j][k] * WAVE_Z;
	    double div_stress[DIM][2];
	    for(int n=0;n<DIM;n++){
	      div_stress[n][0] = 0.0;
	      div_stress[n][1] = 0.0;
	      for(int m=0;m<DIM;m++){
		div_stress[n][0] += (ks[m] * random_stress[m][n][1]);
		div_stress[n][1] += (-ks[m] * random_stress[m][n][0]);
	      }
	    }
	    for(int n=0;n<DIM;n++){
	      for(int d=0;d<2;d++){
		up[n][i][j][k+d] = IRHO * div_stress[n][d] * dmy;
	      }
	    }
	  }
	}
      }
    }
  }
  U_k2u(up);
  Reset_phi(phi);
  Make_phi_u_particle(phi, up, p);
}
