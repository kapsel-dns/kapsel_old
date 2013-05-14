///
// $Id: particle_solver.cxx,v 1.16 2005/09/25 17:46:03 nakayama Exp $
//
#include "particle_solver.h"

void MD_solver_position_Euler(Particle *p, const CTime &jikan){
    for(int n=0; n<Particle_Number; n++){
	for(int d=0; d<DIM; d++){  
	  p[n].x[d] += jikan.dt_md * p[n].v[d];
	  p[n].x[d] = fmod(p[n].x[d]+L_particle[d] , L_particle[d]);

	  assert(p[n].x[d] >= 0);
	  assert(p[n].x[d] < L[d]);

	}
    }
}
void MD_solver_position_CN(Particle *p, const CTime &jikan){
    for(int n=0; n<Particle_Number; n++){
	for(int d=0; d<DIM; d++){  
	  p[n].x[d] += jikan.hdt_md * (p[n].v[d]+p[n].v_old[d]);
	  p[n].x[d] = fmod(p[n].x[d]+L_particle[d] , L_particle[d]);

	  assert(p[n].x[d] >= 0);
	  assert(p[n].x[d] < L[d]);

	}
    }
}
void MD_solver_position_CN_shear(Particle *p, const CTime &jikan){
  const double Shear_max = Shear_rate * LY_shear;
  //fprintf(stdout, "%g %g %g\n", L_particle[0] , L_particle[1] , L_particle[2]); 
  for(int n=0; n<Particle_Number; n++){
    int d;
    {
      d=0;
      p[n].x[d] += jikan.hdt_md * (p[n].v[d]+p[n].v_old[d]);
      p[n].x[d] -= floor(p[n].x[d]*iL_particle[d])*L_particle[d];
      assert(jikan.hdt_md * (p[n].v[0]+p[n].v_old[0])< L_particle[0]);
    }
    {
      d=2;
      p[n].x[d] += jikan.hdt_md * (p[n].v[d]+p[n].v_old[d]);
      p[n].x[d] -= floor(p[n].x[d]*iL_particle[d])*L_particle[d];
      assert(jikan.hdt_md * (p[n].v[2]+p[n].v_old[2])< L_particle[2]);
    }
    {
      d=1;
      p[n].x[d] += jikan.hdt_md * (p[n].v[d]+p[n].v_old[d]);
      assert(jikan.hdt_md * (p[n].v[1]+p[n].v_old[1])< L_particle[1]);
      if(p[n].x[d] > L_particle[d]){
	p[n].x[0] -= Shear_strain;
	p[n].x[0] -= floor(p[n].x[0]*iL_particle[0])*L_particle[0];
	
	p[n].v[0] += Shear_max;
      }else if(p[n].x[d] < 0.0){
	p[n].x[0] += Shear_strain;
	p[n].x[0] -= floor(p[n].x[0]*iL_particle[0])*L_particle[0];

	p[n].v[0] -= Shear_max;
      }
      p[n].x[d] -= floor(p[n].x[d]*iL_particle[d])*L_particle[d];
      //fprintf(stdout, "%g\n", floor(p[n].x[d]*iL_particle[d]));
    }
  }
}
void MD_solver_position_AB2(Particle *p, const CTime &jikan){
  for(int n=0; n<Particle_Number; n++){
    for(int d=0; d<DIM; d++){  
      p[n].x[d] += jikan.hdt_md * (3.*p[n].v[d] - p[n].v_old[d]);
      p[n].x[d] = fmod(p[n].x[d]+L_particle[d] , L_particle[d]);

      assert(p[n].x[d] >= 0);
      assert(p[n].x[d] < L[d]);

    }
  }
}
void MD_solver_velocity_Euler(Particle *p, const CTime &jikan){
  Force(p);

  for(int n=0; n< Particle_Number; n++){
    double dmy = jikan.dt_md * IMASS[p[n].spec];
    double dmy_rot = jikan.dt_md * IMOI[p[n].spec];
    for(int d=0; d<DIM; d++){  
      {
	p[n].v_old[d] = p[n].v[d];
	p[n].v[d] += ( dmy * 
		       ( p[n].f_hydro[d] + p[n].fr[d] + p[n].fv[d] )
		       );
	
	p[n].fr_previous[d] = p[n].fr[d];
	p[n].fr[d] = 0.0;
	p[n].fv_previous[d] = p[n].fv[d];
	p[n].fv[d] = 0.0;
	//p[n].f_hydro_previous[d] = p[n].f_hydro[d];
	p[n].f_hydro[d] = 0.0;
      }
      {
	p[n].omega_old[d] = p[n].omega[d];
	p[n].omega[d] += ( dmy_rot 
			   * ( p[n].torque_hydro[d] 
			       + p[n].torquer[d] 
			       + p[n].torquev[d] )
			   );
	
	p[n].torquer_previous[d] = p[n].torquer[d];
	p[n].torquer[d] = 0.0;
	p[n].torquev_previous[d] = p[n].torquev[d];
	p[n].torquev[d] = 0.0;
	//p[n].torque_hydro_previous[d] = p[n].torque_hydro[d];
	p[n].torque_hydro[d] = 0.0;
      }
      //assert(p[n].v[d]*jikan.dt_md < SIGMA);
      //assert(p[n].v[d]*jikan.dt_md < L_particle[d]);
    }
  }
}
void MD_solver_velocity_Euler_hydro(Particle *p, const CTime &jikan){
  Force(p);

  for(int n=0; n< Particle_Number; n++){
    double dmy = jikan.dt_md * IMASS[p[n].spec];
    double dmy_rot = jikan.dt_md * IMOI[p[n].spec];
    for(int d=0; d<DIM; d++){  
      {
	p[n].v_old[d] = p[n].v[d];
	p[n].v[d] += ( dmy * 
		       ( p[n].f_hydro[d] + p[n].fv[d] 
			 + 0.5*(p[n].fr[d] + p[n].fr_previous[d] )
			 )
		       );
	
	p[n].fr_previous[d] = p[n].fr[d];
	p[n].fr[d] = 0.0;
	p[n].fv_previous[d] = p[n].fv[d];
	p[n].fv[d] = 0.0;
	//p[n].f_hydro_previous[d] = p[n].f_hydro[d];
	p[n].f_hydro[d] = 0.0;
      }
      {
	p[n].omega_old[d] = p[n].omega[d];
	p[n].omega[d] += ( dmy_rot 
			   * ( p[n].torque_hydro[d] 
			       + 0.5 * (p[n].torquer[d] + p[n].torquer[d] )
			       + p[n].torquev[d] )
			   );
	
	p[n].torquer_previous[d] = p[n].torquer[d];
	p[n].torquer[d] = 0.0;
	p[n].torquev_previous[d] = p[n].torquev[d];
	p[n].torquev[d] = 0.0;
	//p[n].torque_hydro_previous[d] = p[n].torque_hydro[d];
	p[n].torque_hydro[d] = 0.0;
      }
      //assert(p[n].v[d]*jikan.dt_md < SIGMA);
      //assert(p[n].v[d]*jikan.dt_md < L_particle[d]);
    }
  }
}
void MD_solver_velocity_AB2_hydro(Particle *p, const CTime &jikan){

  Force(p);

  for(int n=0; n< Particle_Number; n++){
    double dmy = jikan.hdt_md * IMASS[p[n].spec];
    double dmy_rot = jikan.hdt_md * IMOI[p[n].spec];
    for(int d=0; d<DIM; d++){  
      {
	p[n].v_old[d] = p[n].v[d];
	p[n].v[d] += 
	  dmy * (2.*p[n].f_hydro[d]
		 + 3.* p[n].fv[d] - p[n].fv_previous[d] // AB2
		 + p[n].fr[d] + p[n].fr_previous[d] // CN
		 );
      
	p[n].fr_previous[d] = p[n].fr[d];
	p[n].fr[d] = 0.0;
	p[n].fv_previous[d] = p[n].fv[d];
	p[n].fv[d] = 0.0;
	//p[n].f_hydro_previous[d] = p[n].f_hydro[d];
	p[n].f_hydro[d] = 0.0;
      }
      {
	p[n].omega_old[d] = p[n].omega[d];
	p[n].omega[d] += dmy_rot 
	  *( 2.* p[n].torque_hydro[d]
	     + 3.* p[n].torquev[d] -  p[n].torquev_previous[d] // AB2
	     + p[n].torquer[d] + p[n].torquer_previous[d] // CN
	     );
	
	p[n].torquer_previous[d] = p[n].torquer[d];
	p[n].torquer[d] = 0.0;
	p[n].torquev_previous[d] = p[n].torquev[d];
	p[n].torquev[d] = 0.0;
	//p[n].torque_hydro_previous[d] = p[n].torque_hydro[d];
	p[n].torque_hydro[d] = 0.0;
      }
      //assert(p[n].v[d]*jikan.dt_md < SIGMA);
      //assert(p[n].v[d]*jikan.dt_md < L_particle[d]);
    }
  }
}
void MD_solver_velocity_AB2(Particle *p, const CTime &jikan){

  Force(p);

  for(int n=0; n< Particle_Number; n++){
    double dmy = jikan.hdt_md * IMASS[p[n].spec];
    double dmy_rot = jikan.hdt_md * IMOI[p[n].spec];
    for(int d=0; d<DIM; d++){  
      {
	p[n].v_old[d] = p[n].v[d];
	p[n].v[d] += 
	  dmy * (2.*p[n].f_hydro[d]
		 + 3.* (p[n].fv[d] + p[n].fr[d] )
		 - (p[n].fv_previous[d] + p[n].fr_previous[d] )
		 );
      
	p[n].fr_previous[d] = p[n].fr[d];
	p[n].fr[d] = 0.0;
	p[n].fv_previous[d] = p[n].fv[d];
	p[n].fv[d] = 0.0;
	//p[n].f_hydro_previous[d] = p[n].f_hydro[d];
	p[n].f_hydro[d] = 0.0;
      }
      {
	p[n].omega_old[d] = p[n].omega[d];
	p[n].omega[d] += dmy_rot 
	  *( 2.* p[n].torque_hydro[d]
	     + 3.* ( p[n].torquev[d] + p[n].torquer[d] )
	     -  ( p[n].torquev_previous[d] + p[n].torquer_previous[d] )
	     );
	
	p[n].torquer_previous[d] = p[n].torquer[d];
	p[n].torquer[d] = 0.0;
	p[n].torquev_previous[d] = p[n].torquev[d];
	p[n].torquev[d] = 0.0;
	//p[n].torque_hydro_previous[d] = p[n].torque_hydro[d];
	p[n].torque_hydro[d] = 0.0;
      }
      //assert(p[n].v[d]*jikan.dt_md < SIGMA);
      //assert(p[n].v[d]*jikan.dt_md < L_particle[d]);
    }
  }
}
void MD_solver_velocity_sAB2_nohydro(Particle *p, const CTime &jikan){
  static const double Zeta_drag = 6.* M_PI * ETA * RADIUS;
  static const double iZeta_drag = 1./Zeta_drag;
  
  Force(p);
  
  for(int n=0; n< Particle_Number; n++){
    //double dmy = expm1(-Zeta_drag * IMASS[p[n].spec] * jikan.dt_md);
    double dmy = exp(-Zeta_drag * IMASS[p[n].spec] * jikan.dt_md)-1.;
    double dmy1 = -dmy * iZeta_drag *.5;
    for(int d=0; d<DIM; d++){  
      {
	p[n].v_old[d] = p[n].v[d];
	p[n].v[d] = (dmy+1.)*p[n].v[d] 
	  + dmy1 * (3.* (p[n].fr[d] + p[n].fv[d])
		    - (p[n].fr_previous[d] + p[n].fv_previous[d])
		    );
	p[n].fr_previous[d] = p[n].fr[d];
	p[n].fr[d] = 0.0;
	p[n].fv_previous[d] = p[n].fv[d];
	p[n].fv[d] = 0.0;
      }
    }
  }
}
