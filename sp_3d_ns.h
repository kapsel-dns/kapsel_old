//
// $Id: sp_3d_ns.h,v 1.9 2006/11/30 14:25:44 nakayama Exp $
//
#ifndef SP_3D_NS_H
#define SP_3D_NS_H


#define  NDEBUG
#include <assert.h> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "macro.h"
#include "input.h"
#include "variable.h"
#include "make_phi.h"
#include "avs_output.h"
#include "avs_output_p.h"
#include "resume.h"
#include "md_force.h"
#include "fluid_solver.h"
#include "particle_solver.h"
#include "f_particle.h"
#include "init_fluid.h"
#include "init_particle.h"
#include "operate_omega.h"
#include "operate_electrolyte.h"

enum Count_SW {INIT, ADD, MEAN, SNAP_MEAN, SHOW};

inline void Electrolyte_free_energy(const Count_SW &OPERATION
				    ,FILE *fout
				    ,Particle *p
				    ,double **Concentration_rhs1
				    ,const CTime &jikan
				    ){
  static const char *labels[]={""
			       ,"nu"
			       ,"radius"
			       ,"xi"
			       ,"M"
			       ,"I"
			       ,"kBT"
			       ,"kBT/M"
  };
  static char line_label[1<<10];
  static const char *labels_free_energy_nosalt[]={""
					   ,"total"
					   ,"ideal_gas"
					   ,"electrostatic"
					   ,"liquid_charge"
					   ,"total_counterion"
  };
  static const char *labels_free_energy_salt[]={""
					   ,"total"
					   ,"ideal_gas"
					   ,"electrostatic"
					   ,"liquid_charge"
					   ,"total_positive_ion"
					   ,"total_negative_ion"
  };
  //static char line_label_free_energy[1<<10];

  if(OPERATION == INIT){
    {
      sprintf(line_label,"#");
      for(int d=1;d<sizeof(labels)/sizeof(char *);d++){
	sprintf(line_label,"%s%d:%s ",line_label, d, labels[d]);
      }
    }
  }else if(OPERATION == SHOW
	   || OPERATION == MEAN
	   ){
    if(OPERATION == MEAN){
      fprintf(fout,"%s\n",line_label);
      fprintf(fout, "%g %g %g %g %g %g %g\n"
	      ,NU, A*DX, XI*DX
	      ,MASS[p[0].spec]
	      ,MOI[p[0].spec]
	      ,kBT
	      ,kBT*IMASS[p[0].spec]
	      );
    }
    {
      double free_energy[3];
      Calc_free_energy_PB(Concentration_rhs1, p, free_energy, up[0], up[1], up[2], jikan);
      double ion_density = 0.; 
      double *n_solute = new double[N_spec];
      Count_solute_each(n_solute, Concentration_rhs1, p, phi, up[0]);
      for(int n=0;n < N_spec;n++){
	ion_density += n_solute[n] * Valency_e[n];
      }
      char line0[1<<10];
      char line1[1<<10];
      int d;
      int dstart;
      {
	if(OPERATION == SHOW){
	  dstart = 1;
	  sprintf(line0,"#%d:%s ",dstart,"time");
	  sprintf(line1, "%d ",jikan.ts);
	}else {
	  dstart = 0;
	  sprintf(line0,"#");
	  sprintf(line1, "");
	}
      }
      if(N_spec == 1){
	for(d=1;d<sizeof(labels_free_energy_nosalt)/sizeof(char *);d++){
	  sprintf(line0,"%s%d:%s "
		  ,line0, dstart+d, labels_free_energy_nosalt[d]);
	}
	sprintf(line1, "%s%.15g %.15g %.15g %.15g %.15g"
		,line1
		,free_energy[0], free_energy[1], free_energy[2]
		,ion_density, n_solute[0]);
	sprintf(line0,"%s\n",line0);
	sprintf(line1,"%s\n",line1);
      }else{
	for(d=1;d<sizeof(labels_free_energy_salt)/sizeof(char *);d++){
	  sprintf(line0,"%s%d:%s "
		  ,line0, dstart+d, labels_free_energy_salt[d]);
	}
	sprintf(line1, "%s%.15g %.15g %.15g %.15g %.15g %.15g"
		,line1
		,free_energy[0], free_energy[1], free_energy[2]
		,ion_density, n_solute[0], n_solute[1]);
	sprintf(line0,"%s\n",line0);
	sprintf(line1,"%s\n",line1);
      }
      fprintf(fout, "%s%s", line0, line1);
      delete [] n_solute;
    }
  }else {
    fprintf(stderr, "invalid OPERATION in Electrolyte_free_energy().\n"); 
    exit_job(EXIT_FAILURE);
  }

}
inline double Calc_instantaneous_shear_rate(double **zeta
					    ,double uk_dc[DIM]
					    ,double **u //working memory
					    ){
  static const double hivolume = Ivolume * POW3(DX) * 2.;
  static int ny0=NY/4;
  static int ny1=3*NY/4;
  double srate_eff = 0.0;
  
  Zeta_k2u_k(zeta, uk_dc, u);
  A_k2dya_k(u[0], u[1]);
  A_k2a(u[1]);
 {
#pragma omp parallel for reduction(+:srate_eff) 
  for(int i=0; i<NX; i++){
    for(int j=ny0; j<ny1; j++){
      for(int k=0; k<NZ; k++){
	srate_eff += u[1][(i*NY*NZ_)+(j*NZ_)+k];
      }
    }
  }
 }
 return (srate_eff * hivolume);
}

inline void Calc_shear_rate_eff() {
    static const double ivolume = Ivolume * POW3(DX);
    double s_rate_eff = 0.0;
    
    {
#pragma omp parallel for reduction(+:s_rate_eff) 
	
	for(int i = 0; i < NX; i++){
	    for(int j = 0; j < NY - 1; j++){
		for(int k = 0; k < NZ; k++){
		    
		    s_rate_eff +=
			(ucp[0][(i*NY*NZ_)+((j+1)*NZ_)+k] - ucp[0][(i*NY*NZ_)+(j*NZ_)+k])/DX;
		    
		}
	    }
	}
    }

    s_rate_eff *= ivolume*NY/(NY - 1);
    
    Shear_rate_eff = s_rate_eff;
}

inline double Update_strain(double &shear_strain_realized
			    ,const CTime &jikan
			    ,double **zeta
			    ,double uk_dc[DIM]
			    ,double **u //working memory
			    ){
  double srate_eff = -Calc_instantaneous_shear_rate(zeta,uk_dc,u);
  shear_strain_realized += srate_eff * jikan.dt_fluid;
  return srate_eff;
}

inline void Mean_shear_stress(const Count_SW &OPERATION
			      ,FILE *fout
			      ,double *virial
			      ,Particle *p
			      ,const CTime &jikan
			      ,const double &srate_eff
			      ){
    static const int d_virial = 3;
    static const char *labels[]={""
				 ,"time"
				 ,"shear_rate_temporal"
				 ,"shear_strain_temporal"
				 ,"shear_stress_temporal"
				 ,"viscosity"
    };
    static const char *label1s[]={""
				 ,"time"
				 ,"shear_rate_temporal"
				 ,"shear_strain_temporal"
				 ,"shear_stress_temporal"
				 ,"shear_inertia_stress_temporal"
				 ,"apparent_shear_stress"
    };

    static char line_label[1<<10];

    if(OPERATION == INIT){
      if(!Shear_AC){
      {
	  sprintf(line_label,"#");
	  for(int d=1;d<sizeof(labels)/sizeof(char *);d++){
	  sprintf(line_label,"%s%d:%s ", line_label, d, labels[d]);
	  }
      }
      {
	  fprintf(fout,"%s\n",line_label);
      }
	  }else{
      {
	  sprintf(line_label,"#");
	  for(int d=1;d<sizeof(label1s)/sizeof(char *);d++){
	  sprintf(line_label,"%s%d:%s ", line_label, d, label1s[d]);
	  }
      }
      {
	  fprintf(fout,"%s\n",line_label);
      }
	  }

    }else if(OPERATION == SHOW){

	double stress[DIM][DIM] = {
	  {0.,0.,0.},
	  {0.,0.,0.},
	  {0.,0.,0.}
	};

	double hydro_stress[DIM][DIM] = {
	  {0.,0.,0.},
	  {0.,0.,0.},
	  {0.,0.,0.}
	};
	
	Calc_shear_stress(jikan, p, phi, Shear_force, stress);
	Calc_hydro_stress(jikan, p, phi, Hydro_force, hydro_stress);
	
	double strain_output = Shear_strain_realized;
    
	if (SW_EQ == Shear_Navier_Stokes_Lees_Edwards) {
	    fprintf(fout, "%g %g %g %g %g\n"
		    ,jikan.time
		    ,srate_eff
		    ,strain_output
		    ,hydro_stress[1][0]
		    ,(hydro_stress[1][0] + dev_shear_stress_lj)/srate_eff
		);
	} else if(!Shear_AC){
	    fprintf(fout, "%g %g %g %g %g\n"
		    ,jikan.time
		    ,srate_eff
		    ,strain_output
		    ,-stress[1][0]
		    ,-stress[1][0]/srate_eff
		);
	}else{
	    fprintf(fout, "%g %g %g %g %g %g\n"
		    ,jikan.time
		    ,srate_eff
		    ,strain_output
		    ,-stress[1][0]
		    ,Inertia_stress
		    ,-stress[1][0]+Inertia_stress
		);
	}
    }else {
	fprintf(stderr, "invalid OPERATION in Mean_shear_stress().\n"); 
	exit_job(EXIT_FAILURE);
    }
}
#endif

