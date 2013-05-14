//
// $Id: fluid_solver.h,v 1.1 2006/06/27 18:41:28 nakayama Exp $
//
#ifndef FLUID_SOLVER_H
#define FLUID_SOLVER_H

#include <math.h>
#include "variable.h"
#include "input.h"
#include "f_particle.h"
#include "operate_omega.h"
#include "solute_rhs.h"
#include "fluct.h"

extern double *Pressure;
extern double **Shear_force;
extern double **Shear_force_k;
extern double **f_ns0;
extern double **f_ns1;
extern double **f_ns2;
extern double **f_ns3;
extern double **f_ns4;


// Memory 
void Mem_alloc_NS_solver(void);

// Navier_Stokes 
void NS_solver_slavedEuler(double **zeta, const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p);

// Shear_Navier_Stokes
void NS_solver_slavedEuler_Shear_PBC(double **zeta, const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p, double **force);
void NS_solver_slavedEuler_Shear_OBL(double **zeta, const CTime &jikan, double uk_dc[DIM], const Index_range *ijk_range, const int &n_ijk_range, Particle *p, double **force);
//Type of Rheology 
void Mean_shear_sustaining_yforce_PBC(double **u, double **force,const CTime &jikan);
void Mean_shear_sustaining_kforce_PBC(double **u, double **force,const CTime &jikan);
//void Mean_shear_sustaining_kforce_PBC(Value u[DIM], Value force[DIM],const CTime &jikan);
void Mean_shear_sustaining_force_PBC_OBL(double **u);

// ion 
void Ion_diffusion_solver_Euler(double **zeta
				,const CTime &jikan
				,double uk_dc[DIM]
				,double **concentration_k
				,Particle *p
				,const Index_range *ijk_range
				,const int &n_ijk_range
				);

void NSsolute_solver_Euler(double **zeta
			   ,const CTime &jikan
			   ,double uk_dc[DIM]
			   ,double **concentration_k
			   ,Particle *p
			   ,const Index_range *ijk_range
			   ,const int &n_ijk_range
			   );

inline void Calc_shear_stress(const CTime &jikan
			      ,Particle *p
			      ,double *phi
			      //,double **force
			      ,double **Shear_force
			      ,double stress[DIM][DIM]
			      ){

    
    {
	Reset_phi(phi);
	Make_rho_field(phi, p);
	//Make_phi(phi, p);
    }

   int im;
   if(Shear_AC){
   static double factor=Ivolume/jikan.dt_fluid;
   int jm;
   Inertia_stress=0.;
#pragma omp parallel for schedule(dynamic, 1) reduction(+:Inertia_stress) private(im,jm)
   for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
       for(int k=0; k<NZ; k++){
        jm=j;   
		im=(i*NY*NZ_)+(j*NZ_)+k;
        Inertia_stress +=  (double)jm*DX*(phi[im]*ucp[0][im]-rhop[im]); 
        rhop[im] = phi[im]*ucp[0][im]; 
       }
     }
   }
   Inertia_stress *= factor;
  }
    
    double mean_force_x = 0.;
#pragma omp parallel for schedule(dynamic, 1) reduction(+:mean_force_x) private(im)
	  for(int i=0; i<NX; i++){
		for(int j=0; j<NY; j++){
	      for(int k=0; k<NZ; k++){
			im=(i*NY*NZ_)+(j*NZ_)+k;
			//force[0][im] *= phi[im];
			Shear_force[0][im] *= phi[im];
			//mean_force_x += force[0][im];
			mean_force_x += Shear_force[0][im];
		   }
	     }
	   }
	static const double ivolume = Ivolume * POW3(DX);
	mean_force_x *= ivolume;

	 double stress_yx = 0.0;

	double y;
	double fx;
#pragma omp parallel for schedule(dynamic,1) reduction(+:stress_yx) private(y,fx,im)
	    for(int i=0; i<NX; i++){
		for(int j=0; j<NY; j++){
		for(int k=0; k<NZ; k++){
			im=(i*NY*NZ_)+(j*NZ_)+k;
			y = (double)j*DX;
			//fx = force[0][im] - mean_force_x;
			fx = Shear_force[0][im] - mean_force_x;
		    stress_yx += y * fx;
		}
		}
	    }
	 const double dmy = ivolume/jikan.dt_fluid; //Select this 
	 stress[1][0]=stress_yx*dmy;
}

inline void Calc_hydro_stress(const CTime &jikan
			      ,const Particle *p
			      ,double *phi
			      ,double *force
			      ,double stress[DIM][DIM]
    ){
    
    double stress_yx = 0.0;
    
    static const double ivolume = Ivolume * POW3(DX);
#pragma omp parallel for reduction(+:stress_yx) 
    for(int i = 0; i < NX; i++){
	for(int j = 0; j < NY; j++){
	    for(int k = 0; k < NZ; k++){
		/*
		int jj = (int)fmod((double)targetJ + j + NY, NY); 

		double y = (double)j*DX;
		double fx = force[i][jj][k] - mean_force_x;
		stress_yx += y * fx;
		*/
		stress_yx += force[i*NY*NZ_ + j*NZ_ + k];
	    }
	}	
    }
    
    static const double dmy = (double)ivolume/jikan.dt_fluid; //Select this 
    stress[1][0]=stress_yx*dmy;
    
}
#endif


