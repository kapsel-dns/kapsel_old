//
// $Id: fluid_solver.h,v 1.24 2005/09/14 08:41:30 nakayama Exp $
//
#ifndef FLUID_SOLVER_H
#define FLUID_SOLVER_H

#include <math.h>
#include "variable.h"
#include "input.h"
#include "f_particle.h"
#include "operate_omega.h"
#include "operate_Qian_Sheng.h"
#include "solute_rhs.h"
#include "fluct.h"

extern Value Pressure;
extern Value Shear_force[];
extern Value f_ns0[];
extern Value f_ns1[];

void Stokes_solver(Value *zeta, const CTime &jikan, double *uk_dc);
void SP_stokes_advection(Value zeta[DIM-1], const double uk_dc[DIM], Value rhs[DIM-1], const Index_range *ijk_range, const int &n_ijk_range, Particle *p);
void SP_full_advection(Value zeta[DIM-1], const double uk_dc[DIM], Value rhs[DIM-1], const Index_range *ijk_range, const int &n_ijk_range, Particle *p);
void Mem_alloc_NS_solver(void);
void NS_solver_Euler(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p);
void NS_solver_Euler1(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p);
void NS_solver_Euler2(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p);
void NS_solver_Euler3(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p);
void NS_solver_Euler0(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p);
void NS_solver_Heun(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p);
void NS_solver_slavedEuler(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p);
void NS_solver_slavedHeun(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p);
void QS_solver_Euler(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range);
void OG_solver_Euler(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range);
void AC_solver_Euler(CTime &jikan
		     ,Particle *p
		     ,const Index_range *ijk_range
		     ,const int &n_ijk_range);

void Slippy_NS_solver_Euler(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p);
void Slippy_NS_solver_Euler0(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p);
void Shear_NS_solver_Euler(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p, Value force[DIM]);
void Shear_NS_solver_Heun(Value zeta[DIM-1], const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p, Value force[DIM]);
void Ion_diffusion_solver_Euler(Value zeta[DIM-1]
				,const CTime &jikan
				,double uk_dc[DIM]
				,Value *concentration_k
				,Particle *p
				,const Index_range *ijk_range
				,const int &n_ijk_range
				);
void NSsolute_solver_Euler(Value zeta[DIM-1]
			   ,const CTime &jikan
			   ,double uk_dc[DIM]
			   ,Value *concentration_k
			   ,Particle *p
			   ,const Index_range *ijk_range
			   ,const int &n_ijk_range
			   );
void NSsolute_solver_Heun(Value zeta[DIM-1]
			  ,const CTime &jikan
			  ,double uk_dc[DIM]
			  ,Value *concentration_k
			  ,Particle *p
			  ,const Index_range *ijk_range
			  ,const int &n_ijk_range
			  );

inline void Calc_shear_stress(const CTime &jikan
			      ,const Particle *p
			      ,Value &phi
			      ,Value force[DIM]
			      ,double stress[DIM][DIM]
			      ){
  for(int m=0;m<DIM;m++){
    for(int n=0;n<DIM;n++){
      stress[m][n] = 0.0;
    }
  }

  {
    Reset_phi(phi);
    Make_rho_field(phi, p);
  }
  U_k2u(force);
  double mean_force[DIM] = {0.,0.,0.};
  for(int d=0;d<DIM;d++){
    for(int i=0; i<NX; i++){
      for(int j=0; j<HNY; j++){
	for(int k=0; k<NZ; k++){
	  force[d][i][j][k] *= phi[i][j][k];
	  mean_force[d] += force[d][i][j][k];
	}
      }
    }
  }
  //static const double ivolume = 1./((double)(HNY*NX*NZ));
  static const double ivolume = Ivolume * POW3(DX);
  double x[DIM];
  for(int df=0;df<DIM;df++){
    mean_force[df] *= ivolume;
    for(int i=0; i<NX; i++){
      for(int j=0; j<HNY; j++){
	for(int k=0; k<NZ; k++){
	  x[0] = (double)i*DX;
	  x[1] = (double)j*DX;
	  x[2] = (double)k*DX;
	  double f = force[df][i][j][k] - mean_force[df];
	  for(int d=0;d<DIM;d++){
	    stress[d][df] += ( x[d] * f );
	  }
	}
      }
    }
  }
  const double dmy = ivolume/jikan.dt_fluid *.5;
  for(int m=0;m<DIM;m++){
    for(int n=0;n<DIM;n++){
      stress[m][n] *= dmy;
    }
  }
}
#endif

