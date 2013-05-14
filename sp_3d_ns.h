//
// $Id: sp_3d_ns.h,v 1.23 2006/06/11 16:14:30 nakayama Exp $
//
#ifndef SP_3D_NS_H
#define SP_3D_NS_H

//#define  NDEBUG
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
#include "operate_Qian_Sheng.h"
#include "operate_electrolyte.h"

enum Count_SW {INIT, ADD, MEAN, SNAP_MEAN, SHOW};

inline void Check_Poisson_Boltzmann(const Count_SW &OPERATION
				    ,FILE *fout
				    ,Particle *p
				    ,const Value *conc
				    ,Value &e_potential// working memory
				    ,Value &dmy_value0// working memory
				    ,Value &dmy_value1// working memory
				    ){
  if(Particle_Number != 1 
     || SW_EQ  != Electrolyte
     || N_spec != 2
     || Component_Number != 1
     ){
    fprintf(stderr,"###############################\n");
    fprintf(stderr,"# Check_Poisson_Boltzmann mode\n");
    fprintf(stderr,"###############################\n");
    fprintf(stderr, "invalid input parameters.\n"); 
    fprintf(stderr, "set Particle_spec[0].Particle_Number to 1\n");
    fprintf(stderr, "select Electrolyte in constitutive_eq\n");
    fprintf(stderr, "select salt in Add_salt\n");
    exit_job(EXIT_FAILURE);
  }
  static const char *labels[]={""
			       ,"xi"
			       ,"radius"
			       ,"volume_fraction"
			       ,"kappa_a"
			       ,"kBT"
			       ,"permittivity"
			       ,"Bjerrum_length"
			       ,"r"
			       ,"Z"
			       ,"Z+"
			       ,"Z-"
			       ,"Phi(r)"
			       ,"C+(r)"
			       ,"C-(r)"
  };
  static char line_label[1<<10];

  if(OPERATION == INIT){
    {
      sprintf(line_label,"#");
      for(int d=1;d<sizeof(labels)/sizeof(char *);d++){
	sprintf(line_label,"%s%d:%s ",line_label, d, labels[d]);
      }
    }
    {
      p[0].x[0] = HLX;
      p[0].x[1] = HLY;
      p[0].x[2] = HLZ;
    }
    Reset_particle_force_current(p);
    Reset_particle_force_previous(p);
    Reset_particle_velocity(p);
    {
      fprintf(fout,"###############################\n");
      fprintf(fout,"# Check_Poisson_Boltzmann mode\n");
      fprintf(fout,"###############################\n");
      fprintf(fout,"%s\n",line_label);
    }
  }else if(OPERATION == SHOW){
    {
      Conc_k2charge_field(p, conc, e_potential, dmy_value0, dmy_value1);
      A2a_k(e_potential);
      Charge_field_k2Coulomb_potential_k_PBC(e_potential);
      A_k2a(e_potential);
    }

    fprintf(fout,"%s\n",line_label);
    char line_show[1<<10];
    for(int i=HNX;i<NX;i++){
      int j=HNY;
      int k=HNZ;
      sprintf(line_show,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g"
	      ,XI,RADIUS,VF
	      ,RADIUS/Debye_length
	      ,kBT
	      ,Dielectric_cst
	      ,Bjerrum_length
	      ,ABS((double)i*DX-p[0].x[0])
	      ,Surface_charge[0]
	      ,Valency_positive_ion
	      ,Valency_negative_ion
	      ,e_potential[i][j][k]
	      ,conc[0][i][j][k]
	      ,conc[1][i][j][k]
	      );
      fprintf(fout,"%s\n",line_show);
    }
    exit_job(EXIT_SUCCESS);
  }else {
    fprintf(stderr, "invalid OPERATION in Check_Poisson_Boltzmann().\n"); 
    exit_job(EXIT_FAILURE);
  }
}
inline void Electrolyte_free_energy(const Count_SW &OPERATION
				    ,FILE *fout
				    ,Particle *p
				    ,Value *conc
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
      Calc_free_energy_PB(conc, p, free_energy, up[0], up[1], up[2], jikan);
      double ion_density = 0.; 
      //double n_solute[N_spec];
      double *n_solute = new double[N_spec];
      Count_solute_each(n_solute, conc, p, phi, up[0]);
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
inline void DKT(const Count_SW &OPERATION
		,FILE *fout
		,Particle *p
		,const CTime &jikan
		){
    if(Particle_Number != 2 
       || G_direction != 1
       ){
	fprintf(stderr,"###############################\n");
	fprintf(stderr,"# DKT mode\n");
	fprintf(stderr,"###############################\n");
	fprintf(stderr, "invalid input parameters.\n"); 
	fprintf(stderr, "set Particle_spec[0].Particle_Number to 2\n");
	fprintf(stderr, "set G_direction to \"-Y\"\n");
	exit_job(EXIT_FAILURE);
    }

    char *labels[]={""
		    ,"time_step"
		    ,"time"
		    ,"R1x","R1y","R1z"
		    ,"R2x","R2y","R2z"
		    ,"V1x","V1y","V1z"
		    ,"V2x","V2y","V2z"
		    ,"Omega1x","Omega1y","Omega1z"
		    ,"Omega2x","Omega2y","Omega2z"
		    ,"xi"
		    ,"radius"
		    ,"eta"
		    ,"g"
		    ,"rhop-rhof"
		    ,"rhof"
		    ,"Re"
		    ,"volume_fraction"
    };
    static double **x_show = alloc_2d_double(Particle_Number,DIM);
    static int **x_show_int = alloc_2d_int(Particle_Number,DIM);

    if(OPERATION == INIT){
	{
	    p[0].x[0] = HLX;
	    p[0].x[1] = LY - RADIUS;
	    p[0].x[2] = HLZ;

	    p[1].x[0] = HLX;
	    //p[1].x[1] = p[0].x[1] - 1.5*SIGMA;
	    p[1].x[1] = p[0].x[1] - HLY;
	    p[1].x[2] = HLZ;
	}
	Reset_particle_velocity(p);
	Reset_particle_force_current(p);
	Reset_particle_force_previous(p);

	for(int n=0;n< Particle_Number;n++){
	    for(int d=0;d< DIM; d++){
		x_show[n][d] = p[n].x[d];
		x_show_int[n][d] = 0;
	    }
	}
	{
	    fprintf(fout,"###############################\n");
	    fprintf(fout,"# DKT mode\n");
	    fprintf(fout,"###############################\n");
	    fprintf(fout,"#");
	    for(int d=1;d<sizeof(labels)/sizeof(char *);d++){
		fprintf(fout,"%d:%s ", d, labels[d]);
	    }
	    fprintf(fout,"\n");
	}
    }else if(OPERATION == SHOW){
      //double x_old[Particle_Number][DIM];
      double **x_old = alloc_2d_double(Particle_Number,DIM);
      for(int n=0;n< Particle_Number;n++){
	for(int d=0;d< DIM; d++){
	  x_old[n][d] = fmod(x_show[n][d],L_particle[d]);
	  if(x_old[n][d] < 0.){
	    x_old[n][d] += L_particle[d];
	  }
	  if(x_old[n][d]-p[n].x[d] > HL_particle[d]){
	    x_show_int[n][d]--;
	  }else if(p[n].x[d] - x_old[n][d] > HL_particle[d]){
	    x_show_int[n][d]++;
	  }
	  x_show[n][d] = p[n].x[d]-x_show_int[n][d]*L_particle[d];
	}
      }
      //double dmy_v[Particle_Number];
      double *dmy_v = alloc_1d_double(Particle_Number);
      double max_v = 0.;
      for(int n=0;n< Particle_Number;n++){
	double dmy = 0.;
	for(int d=0;d< DIM; d++){
	  dmy += SQ(p[n].v[d]);
	}
	dmy_v[n] = sqrt(dmy);
	max_v = MAX(max_v, dmy_v[n]);
      }
      char line[1<<10];
      sprintf(line, "%d %g",jikan.ts,jikan.ts*jikan.dt_fluid);
      for(int n=0;n< Particle_Number;n++){
	for(int d=0;d< DIM; d++){
	  sprintf(line, "%s %g",line,x_show[n][d]);
	}
      }
      for(int n=0;n< Particle_Number;n++){
	for(int d=0;d< DIM; d++){
	  sprintf(line, "%s %g",line,p[n].v[d]);
	}
      }
      for(int n=0;n< Particle_Number;n++){
	for(int d=0;d< DIM; d++){
	  sprintf(line, "%s %g",line,p[n].omega[d]);
	}
      }
      sprintf(line, "%s %g %g %g %g %g %g %g %g",line
	      ,XI,RADIUS,ETA,G
	      ,(MASS_RATIOS[p[0].spec] - 1.)*RHO
	      ,RHO
	      ,RADIUS/NU*max_v
	      ,VF
	      );
	
      fprintf(fout,"%s\n",line);
      free_2d_double(x_old);
      free_1d_double(dmy_v);
    }else {
	fprintf(stderr, "invalid OPERATION in DKT().\n"); 
	exit_job(EXIT_FAILURE);
    }
}

inline void Terminal_velocity(const Count_SW &OPERATION
			      ,FILE *fout
			      ,Particle *p
			      ,const CTime &jikan
			      ){

  if(Particle_Number != 1 
     ){
      fprintf(stderr,"###############################\n");
      fprintf(stderr,"# Terminal_velocity_measurement mode\n");
      fprintf(stderr,"###############################\n");
      fprintf(stderr, "invalid input parameters.\n"); 
      fprintf(stderr, "set Particle_spec[0].Particle_Number to 1\n");
      exit_job(EXIT_FAILURE);
  }
  static int cnt = 0;
  static double *meanv;
  static double *meanomega;
  static const char *labels[]={""
			       ,"time_step"
			       ,"time"
			       ,"#sample"
			       ,"Vx","Vy","Vz"
			       ,"Omegax","Omegay","Omegaz"
			       ,"xi"
			       ,"radius"
			       ,"eta"
			       ,"g"
			       ,"rhop-rhof"
			       ,"rhof"
			       ,"Re"
			       ,"vterm_dilute"
			       ,"volume_fraction"
  };
  static const char *labels_electrolyte[]={""
					   ,"kappa_a"
					   ,"kBT"
					   ,"Bjerrum_length"
					   ,"counterion_diffusivity"
  };
  static char line_label[1<<10];

  if(OPERATION == INIT){
    {
      sprintf(line_label,"#");
      for(int d=1;d<sizeof(labels)/sizeof(char *);d++){
	sprintf(line_label,"%s%d:%s ",line_label, d, labels[d]);
      }
      if(SW_EQ == Electrolyte){
	int dstart = sizeof(labels)/sizeof(char *)-1;
	for(int d=1;d<sizeof(labels_electrolyte)/sizeof(char *);d++){
	  sprintf(line_label,"%s%d:%s ",line_label, d+dstart, labels_electrolyte[d]);
	}
      }
    }
    meanv = alloc_1d_double(DIM);
    meanomega = alloc_1d_double(DIM);
    for(int d=0;d<DIM;d++){
      meanv[d] = 0.;
      meanomega[d] = 0.;
    }
    {
      //p[0].x[0] = HLX;
      p[0].x[0] = HLX+.5*DX;
      p[0].x[1] = HLY;
      p[0].x[2] = HLZ;
    }
    Reset_particle_velocity(p);
    Reset_particle_force_current(p);
    Reset_particle_force_previous(p);
    {

      fprintf(fout,"###############################\n");
      fprintf(fout,"# Terminal_velocity_measurement mode\n");
      fprintf(fout,"###############################\n");
      fprintf(fout,"%s\n",line_label);
    }
  }else if(OPERATION == SHOW
	   || OPERATION == ADD
	   || OPERATION == MEAN
	   ){
    double v_show[DIM];
    double omega_show[DIM];
    if(OPERATION == ADD){
      for(int d=0;d<DIM;d++){
	meanv[d] += p[0].v[d];
	meanomega[d] += p[0].omega[d];

	v_show[d] = p[0].v[d]; 
	omega_show[d] = p[0].omega[d]; 
      }
      cnt++;
    }else if(OPERATION == MEAN){
      {
	fprintf(fout,"%s\n",line_label);
      }

      double icnt = 1./(double)cnt;
      for(int d=0;d<DIM;d++){
	meanv[d] *= icnt;
	meanomega[d] *= icnt;

	v_show[d] = meanv[d]; 
	omega_show[d] = meanomega[d]; 
      }
    }else if(OPERATION == SHOW){
      for(int d=0;d<DIM;d++){
	v_show[d] = p[0].v[d]; 
	omega_show[d] = p[0].omega[d]; 
      }
    }
    double dmy_v = sqrt(
			SQ(v_show[0])
			+SQ(v_show[1])
			+SQ(v_show[2])
			);
    static const double Stokes_drag_coeff = 6.*M_PI*ETA*RADIUS;
    static const double buoyant_force = 4./3.*M_PI*POW3(RADIUS)*G*(MASS_RATIOS[p[0].spec] - 1.)*RHO;
    static const double vterm_dilute = -buoyant_force/Stokes_drag_coeff;
    char line_show[1<<10];
    sprintf(line_show
	    ,"%d %g %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g"
	    ,jikan.ts
	    ,jikan.ts*jikan.dt_fluid
	    ,cnt
	    ,v_show[0],v_show[1],v_show[2]
	    ,omega_show[0],omega_show[1],omega_show[2]
	    ,XI,RADIUS,ETA,G
	    ,(MASS_RATIOS[p[0].spec] - 1.)*RHO
	    ,RHO
	    ,RADIUS/NU*dmy_v
	    ,vterm_dilute
	    ,VF
	    );
    if(SW_EQ == Electrolyte){
      if(N_spec==1){
	sprintf(line_show
		,"%s %g %g %g %g"
		,line_show
		,RADIUS/ Debye_length
		,kBT
		,Bjerrum_length
		,kBT*Onsager_coeff_counterion
		);
      }else if(N_spec==2){
	sprintf(line_show
		,"%s %g %g %g %g %g"
		,line_show
		,RADIUS/ Debye_length
		,kBT
		,Bjerrum_length
		,kBT*Onsager_coeff_positive_ion
		,kBT*Onsager_coeff_negative_ion
		);
      }
    }
    fprintf(fout,"%s\n",line_show);

  }else {
    fprintf(stderr, "invalid OPERATION in Terminal_velocity().\n"); 
    exit_job(EXIT_FAILURE);
  }
}
inline void Lubrication(const Count_SW &OPERATION
			,FILE *fout
			,Particle *p
			,const CTime &jikan
			){
  if(Component_Number != 2 
     || Particle_Number != 2 
     || G_direction != 0
     ){
      fprintf(stderr,"###############################\n");
      fprintf(stderr,"# Lubrication_measurement mode\n");
      fprintf(stderr,"###############################\n");
      fprintf(stderr, "invalid input parameters.\n"); 
      fprintf(stderr, "set G_direction to -X\n"); 
      fprintf(stderr, "set Particle_spec[0].Particle_Number to 1\n"); 
      fprintf(stderr, "set Particle_spec[1].Particle_Number to 1\n"); 
    exit_job(EXIT_FAILURE);
  }
  char *labels[]={""
		  ,"time_step"
		  ,"time"
		  ,"(x2-x1)-r_eq"
		  ,"(y2-y1)-r_eq"
		  ,"(z2-z1)-r_eq"
		  ,"v2x-v1x"
		  ,"v2y-v1y"
		  ,"v2z-v1z"
		  ,"particle_volume"
		  ,"rho1-rhof"
		  ,"rho2-rhof"
		  ,"rhof"
		  ,"g_x"
		  ,"Re"
  };
  if(OPERATION == INIT){
    {
      double separation_half;
      if(MASS_RATIOS[p[0].spec] < MASS_RATIOS[p[1].spec]){ 
	separation_half = HLX*.5;
      }else {
	  separation_half = LJ_dia*pow(2.,1./6.)*.5;
      }
      p[0].x[0] = HLX - separation_half;
      p[0].x[1] = HLY;
      p[0].x[2] = HLZ;

      p[1].x[0] = HLX + separation_half;
      p[1].x[1] = HLY;
      p[1].x[2] = HLZ;
    }
    Reset_particle_velocity(p);
    Reset_particle_force_current(p);
    Reset_particle_force_previous(p);
    {
      fprintf(fout,"###############################\n");
      fprintf(fout,"# Lubrication_measurement mode");
      if(MASS_RATIOS[p[0].spec] < MASS_RATIOS[p[1].spec]){ 
	fprintf(fout," :approaching\n");
      }else {
	fprintf(fout," :retracting\n");
      }
      fprintf(fout,"###############################\n");
      fprintf(fout,"#");
      for(int d=1;d<sizeof(labels)/sizeof(char *);d++){
	fprintf(fout,"%d:%s ", d, labels[d]);
      }
      fprintf(fout,"\n");
    }
  }else if(OPERATION == SHOW){
    double r_eq = LJ_dia * pow(2.,1./6.);
    double dmy_v = sqrt(
			SQ(p[1].v[0]-p[0].v[0])
			+SQ(p[1].v[1]-p[0].v[1])
			+SQ(p[1].v[2]-p[0].v[2])
			);
    fprintf(fout,"%d %g %g %g %g %g %g %g %g %g %g %g %g %g\n"
	    ,jikan.ts
	    ,jikan.ts*jikan.dt_fluid
	    ,(p[1].x[0]-p[0].x[0])-r_eq
	    ,(p[1].x[1]-p[0].x[1])-r_eq
	    ,(p[1].x[2]-p[0].x[2])-r_eq
	    ,(p[1].v[0]-p[0].v[0])
	    ,(p[1].v[1]-p[0].v[1])
	    ,(p[1].v[2]-p[0].v[2])
	    ,4./3.*M_PI*POW3(RADIUS)
	    ,(MASS_RATIOS[p[0].spec] - 1.)*RHO
	    ,(MASS_RATIOS[p[1].spec] - 1.)*RHO
	    ,RHO
	    ,G
	    ,RADIUS/NU*dmy_v
	    );

  }else {
    fprintf(stderr, "invalid OPERATION in Lubrication().\n"); 
    exit_job(EXIT_FAILURE);
  }
}
inline void Mean_shear_stress(const Count_SW &OPERATION
			      ,FILE *fout
			      ,const double *virial
			      ,const Particle *p
			      ,const CTime &jikan
			      ,const Value zeta[DIM-1]
			      ,const double uk_dc[DIM]
			      ,Value u[DIM] //working memory
			      ){
    static double *sum_shear_stress;
    static const int d_virial = 3;
    static double **sum_shear_stress_hydro;
    static double sum_shear_rate = 0.0;
    static int cnt = 0;

    static const char *labels[]={""
				 ,"dev_sigma_yx"
				 ,"shear_rate_obs"
				 ,"shear_rate"
				 ,"ETA"
				 ,"VF"
				 ,"stress_Einstein"
				 ,"sigma_hydro_xx"
				 ,"sigma_hydro_yx"
				 ,"sigma_hydro_zx"
				 ,"sigma_hydro_xy"
				 ,"sigma_hydro_yy"
				 ,"sigma_hydro_zy"
				 ,"sigma_hydro_xz"
				 ,"sigma_hydro_yz"
				 ,"sigma_hydro_zz"
				 ,"cnt"
				 ,"strain"
				 ,"shear_rate_obs*LJ_time"
				 ,"shear_Re"
				 ,"dev_sigma_lub"
				 ,"dev_sigma_LJ"
    };
    static char line_label[1<<10];

    if(OPERATION == INIT){
      sum_shear_stress = alloc_1d_double(d_virial);
      for(int m=0;m<d_virial;m++){
	sum_shear_stress[m] = 0.;
      }
      sum_shear_stress_hydro = alloc_2d_double(DIM,DIM);
      for(int m=0;m<DIM;m++){
	for(int n=0;n<DIM;n++){
	  sum_shear_stress_hydro[m][n] = 0.0;
	}
      }
      {
	sprintf(line_label,"#");
	for(int d=1;d<sizeof(labels)/sizeof(char *);d++){
	  sprintf(line_label,"%s%d:%s ", line_label, d, labels[d]);
	}
      }
      {
	fprintf(fout,"%s\n",line_label);
      }
    }else if(OPERATION == ADD || OPERATION == SHOW){// conuting
	double srate_eff = 0.0;
	{
	    Zeta_k2u_k(zeta, uk_dc, u);
	    A_k2dya_k(u[0], u[1]);
	    A_k2a(u[1]);

	    for(int i=0; i<NX; i++){
	      for(int j=0; j<HNY; j++){
		for(int k=0; k<NZ; k++){
		  srate_eff += u[1][i][j][k];
		}
	      }
	    }
	    //static const double ivolume = 1./((double)(HNY*NX*NZ));	
	    static const double ivolume = Ivolume * POW3(DX);
	    srate_eff *= ivolume;
	    if(OPERATION == ADD){
	      sum_shear_rate += srate_eff;
	    }
	}
	double stress[DIM][DIM] = {
	  {0.,0.,0.},
	  {0.,0.,0.},
	  {0.,0.,0.}
	};
	if(HYDRO_int > 0){
	  Calc_shear_stress(jikan, p, phi, Shear_force, stress);
	}
	fprintf(fout, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d %g %g %g %g %g\n"
		,virial[0]
		,srate_eff
		,Shear_rate, ETA, VF
		,ETA*Shear_rate*(1.+2.5*VF)
		,stress[0][0], stress[1][0], stress[2][0]
		,stress[0][1], stress[1][1], stress[2][1]
		,stress[0][2], stress[1][2], stress[2][2]
		,cnt
		,Shear_strain + Shear_strain_int*L[0]
		,srate_eff*T_LJ
		,srate_eff*SQ(SIGMA)/NU
		,virial[1]
		,virial[2]
		);

	if(OPERATION == ADD){
	  for(int d=0;d<d_virial;d++){
	    sum_shear_stress[d] += virial[d];
	  }
	  for(int m=0;m<DIM;m++){
	    for(int n=0;n<DIM;n++){
	      sum_shear_stress_hydro[m][n] += stress[m][n];
	    }
	  }
	  cnt++; 
	}
    }else if(OPERATION == MEAN){ // returning averages
	double icnt = 1./(double)cnt;
	double mean[d_virial];
	for(int d=0;d<d_virial;d++){
	  mean[d] = sum_shear_stress[d] * icnt;
	}
	double mean_shear_rate = sum_shear_rate * icnt;
	double mean_hydro[DIM][DIM];
	for(int m=0;m<DIM;m++){
	  for(int n=0;n<DIM;n++){
	    mean_hydro[m][n] = sum_shear_stress_hydro[m][n] * icnt;
	  }
	}
	{
	  fprintf(fout,"%s\n",line_label);
	}

	{
	  fprintf(fout, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d %g %g %g %g %g\n"
		    ,mean[0]
		    ,mean_shear_rate
		    ,Shear_rate
		    ,ETA
		    ,VF
		    ,ETA*Shear_rate*(1.+2.5*VF)
		    ,mean_hydro[0][0],mean_hydro[1][0],mean_hydro[2][0]
		    ,mean_hydro[0][1],mean_hydro[1][1],mean_hydro[2][1]
		    ,mean_hydro[0][2],mean_hydro[1][2],mean_hydro[2][2]
		    ,cnt
		    ,mean_shear_rate*jikan.time
		    ,mean_shear_rate*T_LJ
		    ,mean_shear_rate*SQ(SIGMA)/NU
		    ,mean[1]
		    ,mean[2]
		    );
	}
    }else {
	fprintf(stderr, "invalid OPERATION in Mean_shear_stress().\n"); 
	exit_job(EXIT_FAILURE);
    }
}

inline void Count_noisy_NS(const Count_SW &OPERATION
			   ,FILE *fout
			   ,Particle *p
			   ,const double uk_dc[DIM]
			   ,const Value zeta[DIM-1]
			   ,Value &phi // working memory
			   ,Value u[DIM]// working memory
			   ,double return_v[DIM]
			   ){
  static int cnt_f = 0;
  static double *fluid_v_mean;
  static double *fluid_vsq_mean;
  
  static const char *labels[]={""
			       ,"nu"
			       ,"RHO"
			       ,"kBT"
			       ,"DX"
			       ,"Lx"
			       ,"Ly"
			       ,"Lz"
			       ,"#sample"
			       ,"vx","vx^2","dvx^2"
			       ,"vy","vy^2","dvy^2"
			       ,"vz","vz^2","dvz^2"
  };
  static char line_label[1<<10];

  if(OPERATION == INIT){
    {
      fluid_v_mean = alloc_1d_double(DIM);
      fluid_vsq_mean = alloc_1d_double(DIM);
      for(int d=0;d<DIM;d++){
	fluid_v_mean[d] = 0.0;
	fluid_vsq_mean[d] = 0.0;
      }
    }
    {
      sprintf(line_label,"#");
      for(int d=1;d<sizeof(labels)/sizeof(char *);d++){
	sprintf(line_label,"%s%d:%s ",line_label, d, labels[d]);
      }
    }
  }else if(OPERATION == ADD){
    double base_v[DIM];
    //const double ivolume = 1./(double)(NX*NY*NZ);
    const double ivolume = Ivolume * POW3(DX);
    for(int d=0;d<DIM;d++){
      base_v[d] = uk_dc[d] * ivolume;
      //  base_v[d] = 0.;
    }
    Zeta_k2u(zeta, uk_dc, u); 
    Reset_phi(phi);
    Make_phi_particle(phi, p);
    
    double v_dmy[DIM];
    double vsq_dmy[DIM];
    double fluid_volume=0.0;
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ; k++){
	  double dmy_phi = phi[i][j][k];
	  fluid_volume += (1.-dmy_phi);
	}
      }
    }
    for(int d=0;d<DIM;d++){
      v_dmy[d] = vsq_dmy[d] = 0.;
      for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	  for(int k=0; k<NZ; k++){
	    double dmy_phi = phi[i][j][k];
	    double dmy_u = u[d][i][j][k] - base_v[d];
	    v_dmy[d] += dmy_u * ( 1- dmy_phi);
	    vsq_dmy[d] += SQ(dmy_u) * ( 1- dmy_phi);
	  }
	}
      }
      double ifluid_volume = 1./fluid_volume;
      fluid_v_mean[d] += (v_dmy[d] * ifluid_volume);
      fluid_vsq_mean[d] += (vsq_dmy[d] * ifluid_volume);
      return_v[d] = fluid_v_mean[d];
    }
    cnt_f++; 
  }else if(OPERATION == MEAN){
    {
      fprintf(fout,"%s\n",line_label);
    }
    {
      fprintf(fout , "%g %g %g %g %g %g %g %d"
	      ,NU,RHO,kBT,DX,LX,LY,LZ,cnt_f);
      double icnt_f = 1./(double)cnt_f;
      for(int d=0;d<DIM;d++){
	fluid_v_mean[d] *= icnt_f;
	fluid_vsq_mean[d] *= icnt_f;
	fprintf(fout, " %g %g %g"
		,fluid_v_mean[d],fluid_vsq_mean[d]
		,fluid_vsq_mean[d]-SQ(fluid_v_mean[d]));
      }
      fprintf(fout, "\n");
    }
  }else {
    fprintf(stderr, "invalid OPERATION in Count_noisy_NS().\n"); 
    exit_job(EXIT_FAILURE);
  }
}
inline void Count_noisy_particle(const Count_SW &OPERATION
				 ,FILE *fout
				 ,const double base_v[DIM]
				 ,const Particle *p
				 ,const CTime &jikan
				 ){
  static double *particle_v_mean;
  static double *particle_vsq_mean;
  static double *particle_omega_mean;
  static double *particle_omegasq_mean;
  static double *particle_potential_energy;
  static int cnt = 0; 

  static int n_relax;
  static double *vp_mean_n;
  static double *op_mean_n;
  static double *vf_mean_n;
  static double *vpsq_mean_n;
  static double *opsq_mean_n;
  static double *vfsq_mean_n;

  static int temp_current=0;
  
  static const char *labels[]={""
			       ,"nu"
			       ,"radius"
			       ,"xi"
			       ,"M"
			       ,"I"
			       ,"kBT"
			       ,"kBT/M"
			       ,"#sample"
			       ,"Vx","Vx^2","dVx^2"
			       ,"Omegax","Omegax^2","dOmegax^2"
			       ,"Vy","Vy^2","dVy^2"
			       ,"Omegay","Omegay^2","dOmegay^2"
			       ,"Vz","Vz^2","dVz^2"
			       ,"Omegaz","Omegaz^2","dOmegaz^2"
			       ,"M_eff"
			       ,"Ux"
			       ,"Uy"
			       ,"Uz"
			       ,"U"
  };
  static char line_label[1<<10];
  if(OPERATION == INIT){
    {
      particle_v_mean = alloc_1d_double(DIM);
      particle_vsq_mean = alloc_1d_double(DIM);
      particle_omega_mean = alloc_1d_double(DIM);
      particle_omegasq_mean = alloc_1d_double(DIM);
      particle_potential_energy = alloc_1d_double(DIM);
      for(int d=0;d<DIM;d++){
	particle_v_mean[d] = 0.0;
	particle_vsq_mean[d] = 0.0;
	particle_omega_mean[d] = 0.0;
	particle_omegasq_mean[d] = 0.0;
	particle_potential_energy[d] = 0.0;
      }
    }
    
    {
      double mass_max = 0.0;
      for(int i=0; i<Component_Number; i++){
	mass_max = MAX(mass_max, MASS[i]);
      }
      const double Zeta_drag = 6.* M_PI * ETA * RADIUS;
      n_relax =(int)(mass_max/Zeta_drag/jikan.dt_fluid);
      int min_n_relax = 100;
      n_relax = n_relax>min_n_relax?n_relax:min_n_relax;
    }
    {
      vp_mean_n = alloc_1d_double(n_relax);
      op_mean_n = alloc_1d_double(n_relax);
      vf_mean_n = alloc_1d_double(n_relax);
      vpsq_mean_n = alloc_1d_double(n_relax);
      opsq_mean_n = alloc_1d_double(n_relax);
      vfsq_mean_n = alloc_1d_double(n_relax);
    }
    {
      for(int n=0;n<n_relax;n++){
	vp_mean_n[n] = 0.0;
	op_mean_n[n] = 0.0;
	vf_mean_n[n] = 0.0;
	vpsq_mean_n[n] = 0.0;
	opsq_mean_n[n] = 0.0;
	vfsq_mean_n[n] = 0.0;
      }
    }
    {
      sprintf(line_label,"#");
      for(int d=1;d<sizeof(labels)/sizeof(char *);d++){
	sprintf(line_label,"%s%d:%s ",line_label, d, labels[d]);
      }
    }
  }else if(OPERATION == ADD){
    double potential_energy[DIM];
    Calc_harmonic_energy(p, potential_energy);
    for(int d=0;d<DIM;d++){
      particle_potential_energy[d] += potential_energy[d];
    }
    double dmy_v_mean[DIM] = {0.,0.,0.};
    double dmy_vsq_mean[DIM] = {0.,0.,0.};
    double dmy_omega_mean[DIM] = {0.,0.,0.};
    double dmy_omegasq_mean[DIM] = {0.,0.,0.};

    for(int n=0; n<Particle_Number;n++){
      for(int d=0;d<DIM;d++){
	double dmy_v = p[n].v[d] - base_v[d];
	dmy_v_mean[d] += dmy_v;
	dmy_vsq_mean[d] += SQ(dmy_v);
	dmy_omega_mean[d] += p[n].omega[d];
	dmy_omegasq_mean[d] += SQ(p[n].omega[d]);
      }
    }
    if(0){
      fprintf(stdout, "%g %g %g #ccc:\n"
	      ,base_v[0]
	      ,base_v[1]
	      ,base_v[2]
	      );
    }
    vp_mean_n[temp_current] = 0.;
    op_mean_n[temp_current] = 0.;
    vpsq_mean_n[temp_current] = 0.;
    opsq_mean_n[temp_current] = 0.;
    double iNp = 1./(double)Particle_Number;  
    for(int d=0;d<DIM;d++){
      dmy_v_mean[d] *= iNp;
      dmy_vsq_mean[d] *= iNp;
      dmy_omega_mean[d] *= iNp;
      dmy_omegasq_mean[d] *= iNp;
      
      vp_mean_n[temp_current] +=  dmy_v_mean[d];
      vpsq_mean_n[temp_current] +=  dmy_vsq_mean[d];
      op_mean_n[temp_current] +=  dmy_omega_mean[d];
      opsq_mean_n[temp_current] +=  dmy_omegasq_mean[d];
      
      particle_v_mean[d] += dmy_v_mean[d];
      particle_vsq_mean[d] += dmy_vsq_mean[d];
      particle_omega_mean[d] += dmy_omega_mean[d];
      particle_omegasq_mean[d] += dmy_omegasq_mean[d];
    }
    cnt++;
    
    {
      vp_mean_n[temp_current] *= One_third;
      vpsq_mean_n[temp_current] *= One_third;
      op_mean_n[temp_current] *= One_third;
      opsq_mean_n[temp_current] *= One_third;
      
      temp_current++;
      temp_current %= n_relax;
    }
  }else if(OPERATION == SNAP_MEAN){
    {
      fprintf(fout, "%g %g %g %g %g %g %g %d"
	      ,NU, A*DX, XI*DX
	      ,MASS[p[0].spec]
	      ,MOI[p[0].spec]
	      ,kBT
	      ,kBT*IMASS[p[0].spec]
	      ,cnt);
      //double icnt = 1./((double)(cnt*Particle_Number));
      double icnt = 1./((double)(cnt));
      for(int d=0;d<DIM;d++){
	double dmy_v,dmy_vsq,dmy_omega,dmy_omegasq;
	dmy_v = particle_v_mean[d] * icnt;
	dmy_vsq = particle_vsq_mean[d] * icnt;
	dmy_omega = particle_omega_mean[d] * icnt;
	dmy_omegasq = particle_omegasq_mean[d] * icnt;
	fprintf(fout, " %g %g %g %g %g %g"
		,dmy_v
		,dmy_vsq
		,dmy_vsq - SQ(dmy_v)
		,dmy_omega
		,dmy_omegasq
		,dmy_omegasq-SQ(dmy_omega));
      }
      fprintf(fout, " %g"
	      ,MASS[p[0].spec]*(
				1. + 0.5 * IMASS_RATIOS[p[0].spec])
	      );
      double dmy_total_potential=0.0; 
      for(int d=0;d<DIM;d++){
	double dmy_potential[DIM];
	dmy_potential[d] = particle_potential_energy[d] * icnt;
	dmy_total_potential += dmy_potential[d];
	fprintf(fout, " %g", dmy_potential[d]);
      }
      fprintf(fout, " %g", dmy_total_potential);
      fprintf(fout, "\n");
    }	
  }else if(OPERATION == MEAN){
    {
      fprintf(fout,"%s\n",line_label);
    }
    {
      fprintf(fout, "%g %g %g %g %g %g %g %d"
	      ,NU, A*DX, XI*DX
	      ,MASS[p[0].spec]
	      ,MOI[p[0].spec]
	      ,kBT
	      ,kBT*IMASS[p[0].spec]
	      ,cnt);
      //double icnt = 1./((double)(cnt*Particle_Number));
      double icnt = 1./((double)(cnt));
      for(int d=0;d<DIM;d++){
	particle_v_mean[d] *= icnt;
	particle_vsq_mean[d] *= icnt;
	particle_omega_mean[d] *= icnt;
	particle_omegasq_mean[d] *= icnt;
	fprintf(fout, " %g %g %g %g %g %g"
		,particle_v_mean[d]
		,particle_vsq_mean[d]
		,particle_vsq_mean[d]-SQ(particle_v_mean[d])
		,particle_omega_mean[d]
		,particle_omegasq_mean[d]
		,particle_omegasq_mean[d]-SQ(particle_omega_mean[d]));
      }
      fprintf(fout, " %g"
	      ,MASS[p[0].spec]*(
				1. + 0.5 * IMASS_RATIOS[p[0].spec])
	      );
      double dmy_potential=0.0; 
      for(int d=0;d<DIM;d++){
	particle_potential_energy[d] *= icnt;
	dmy_potential += particle_potential_energy[d];
	fprintf(fout, " %g", particle_potential_energy[d]);
      }
      fprintf(fout, " %g", dmy_potential);
      fprintf(fout, "\n");
    }
  }else {
    fprintf(stderr, "invalid OPERATION in Count_noisy_particle().\n"); 
    exit_job(EXIT_FAILURE);
  }
}

inline void Count_Slippy_NS(const Particle *p, const Count_SW &OPERATION){
  static double *particle_v_mean;
  static double *particle_vsq_mean;
  static double *particle_omega_mean;
  static double *particle_omegasq_mean;
  static int cnt = 0;

  if(OPERATION == INIT){
    {
      particle_v_mean = alloc_1d_double(DIM);
      particle_vsq_mean = alloc_1d_double(DIM);
      particle_omega_mean = alloc_1d_double(DIM);
      particle_omegasq_mean = alloc_1d_double(DIM);
      for(int d=0;d<DIM;d++){
	particle_v_mean[d] = 0.0;
	particle_vsq_mean[d] = 0.0;
	particle_omega_mean[d] = 0.0;
	particle_omegasq_mean[d] = 0.0;
      }
    }
  }else if(OPERATION == ADD){
    for(int n=0; n<Particle_Number;n++){
      for(int d=0;d<DIM;d++){
	particle_v_mean[d] += p[n].v[d];
	particle_vsq_mean[d] += SQ(p[n].v[d]);
	particle_omega_mean[d] += p[n].omega[d];
	particle_omegasq_mean[d] += SQ(p[n].omega[d]);
      }
    }
    cnt++; 
  }else if(OPERATION == MEAN){
    char *labels[]={""
		    ,"delta_rho"
		    ,"g"
		    ,"nu"
		    ,"nu_ratio"
		    ,"radius"
		    ,"xi"
		    ,"M"
		    ,"I"
		    ,"kBT"
		    ,"kBT/M"
		    ,"#sample"
		    ,"Vx","gamma_x"
		    ,"Vy","gamma_y"
		    ,"Vz","gamma_z"
    };
    {
      fprintf(stderr,"#");
      for(int d=1;d<sizeof(labels)/sizeof(char *);d++){
	fprintf(stderr,"%d:%s ", d, labels[d]);
      }
      fprintf(stderr,"\n");
    }
    {
      fprintf(stderr, "%g %g %g %g %g %g %g %g %g %g %d"
	      ,RHO_particle[p[0].spec]-RHO
	      ,G
	      ,NU
	      ,Nu_ratio
	      ,A*DX, XI*DX
	      ,MASS[p[0].spec]
	      ,MOI[p[0].spec]
	      ,kBT
	      ,kBT*IMASS[p[0].spec]
	      ,cnt);
      double icnt = 1./(double)cnt/(double)Particle_Number;
      double volume = 4./3.*M_PI*POW3(RADIUS);
      for(int d=0;d<DIM;d++){
	particle_v_mean[d] *= icnt;
	fprintf(stderr, " %g %g"
		,particle_v_mean[d]
		,(RHO_particle[p[0].spec]-RHO)*volume*G/particle_v_mean[d]);
      }
      fprintf(stderr, "\n");
    }
  }else {
    fprintf(stderr, "invalid OPERATION in Count_Slippy_NS().\n"); 
    exit_job(EXIT_FAILURE);
  }
}

inline double Particle_Reynolds_number(Particle *p){
  static double vmax;
  vmax = 0.0;
  for(int n=0; n<Particle_Number; n++){
    for(int d=0; d<DIM; d++){  
	vmax = MAX(ABS(p[n].v[d]), vmax);
    }
  }
  return vmax * RADIUS/NU;
}


#endif
