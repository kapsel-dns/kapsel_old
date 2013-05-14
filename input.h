//
// $Id: input.h,v 1.3 2006/11/14 20:07:12 nakayama Exp $
//
#ifndef INPUT_H
#define INPUT_H

#include <stdio.h>
#include <iostream>
#include <string>
#include <cstring>
#include <cfloat>
#include <sys/types.h>
#include <dirent.h>
#include "macro.h"
#include "alloc.h"
#include "udfmanager.h"
#include "parameter_define.h"

/////////////////////
/////////////////////
extern int Fixed_particle;
/////////////////////
/////////////////////

enum Particle_IC {
  None
  ,uniform_random
  ,random_walk
  ,FCC
  ,BCC
  ,user_specify
};
enum SW_time {
    AUTO, MANUAL
};
enum EQ {Navier_Stokes
	 ,Shear_Navier_Stokes
	 ,Shear_Navier_Stokes_Lees_Edwards
	 ,Electrolyte
};
enum PT {spherical_particle
	 ,chain
};
//////  
extern SW_time SW_TIME;
//////  
extern EQ SW_EQ;
extern char *EQ_name[];
//////  
extern PT SW_PT;
extern char *PT_name[];
//////  material parameters
extern double RHO;
extern double IRHO;
extern double ETA;
extern double kBT;
extern double ikBT;
extern double Shear_rate;
extern double Shear_rate_eff;
extern double Shear_strain_realized;
extern double Shear_strain;
extern double Shear_frequency;
extern double Inertia_stress;
extern int Shear_strain_int;
extern double dev_shear_stress[];
extern double &dev_shear_stress_lj;
extern double Delta_ETA;
extern double Nu_ratio;
extern double NU;
extern double *MASS_RATIOS;
extern double *S_surfaces;// pretilt scalar order
extern double *W_surfaces;// spring cst. of anchoring
extern double IMASS_RATIO;
extern double *IMASS_RATIOS;
extern double *RHO_particle;
extern double *MASS;
extern double *MOI;
extern double *IMASS;
extern double *IMOI;
extern double EPSILON ;
extern double T_LJ;
extern int LJ_powers ;
extern int RESUMED ;
extern int last_ts ;
extern double Srate_depend_LJ_cap;
extern double LJ_dia;
//////
/////// avs output;
extern int SW_AVS;
extern char Out_dir[];
extern char Out_name[];
extern int BINARY;
//////
extern int SW_UDF;
////// FFT select
extern int SW_FFT;
//////

/////// 計算条件の設定
extern int Nmax;
extern int Nmin;
extern int Ns[];
extern int Ns_shear[];
extern int HNs[];
extern int N2s[];
extern int HN2s[];
extern int TRNs[];
extern int TRNs_QS[];
extern int &NX;
extern int &NY;
extern int &NZ;
extern int HNZ_;
extern int NZ_;
extern int HN2Z_;
extern int N2Z_;
extern int &HNX;
extern int &HNY;
extern int &HNZ;
extern int &N2X;
extern int &N2Y;
extern int &N2Z;
extern int &HN2X;
extern int &HN2Y;
extern int &HN2Z;
extern int &TRN_X;
extern int &TRN_Y;
extern int &TRN_Z;
extern int &TRN_QS_X;
extern int &TRN_QS_Y;
extern int &TRN_QS_Z;

//////
extern int ROTATION;
extern int STOKES;
extern int LJ_truncate;
extern Particle_IC DISTRIBUTION;
extern int N_iteration_init_distribution;
extern int FIX_CELL;
extern int FIX_CELLxyz[DIM];
//////
extern double DX;
extern double DX3;
extern double A_XI;
extern double XI;
extern double HXI ;
extern double A;
extern double RADIUS;
extern double SIGMA ;
extern double A_R_cutoff ;
extern double R_cutoff ;
//////
extern double G;
extern int G_direction;
//////
extern int Component_Number;
extern int Particle_Number;
extern int *Particle_Numbers;
extern int *Beads_Numbers;
extern int *Chain_Numbers;
extern double VF;
extern double VF_LJ;
extern double Ivolume;
//////
extern int GTS;
extern int Num_snap;
extern int MSTEP;
//////
extern double &LX;
extern double &LY;
extern double &LZ;
extern double &HLX;
extern double &HLY;
extern double &HLZ;
extern double L[];
extern double HL[];
extern double L_particle[];
extern double iL_particle[];
extern double HL_particle[];
extern double WAVE_X;
extern double WAVE_Y;
extern double WAVE_Z;
extern double KMAX2;
//////
extern double Axel;
extern double Tdump;
extern double DT_noise;
extern double DT;
/////// Two_fluid
extern double Mean_Bulk_concentration;
extern int N_spec;
extern double Onsager_solute_coeff;
/////// Electrolyte
extern int Poisson_Boltzmann;
extern int External_field;
extern int AC;
extern int Shear_AC;
extern double *Surface_charge;
extern double *Surface_charge_e;
extern double Elementary_charge;
extern double Valency_counterion;
extern double Valency_positive_ion;
extern double Valency_negative_ion;
extern double Ion_charge;
extern double Onsager_coeff_counterion;
extern double Onsager_coeff_positive_ion;
extern double Onsager_coeff_negative_ion;
extern double Dielectric_cst;
extern double Debye_length;
extern double E_ext[DIM];
extern double Frequency;
extern double Angular_Frequency;
////// Temp Monitor

extern double kT_snap_v; //fix R4A2
extern double kT_snap_o; //fix R4A2

extern double alpha_v;
extern double alpha_o;

//////shear_degree
extern double degree_oblique;

//////

extern char *In_udf,*Sum_udf,*Out_udf,*Def_udf,*Ctrl_udf,*Res_udf;
extern UDFManager *ufin, *ufout, *ufres;
//extern UDFManager *ufsum;
void file_get(const int argc, char *argv[]);
void Gourmet_file_io(const char *infile
		     ,const char *outfile
		     ,const char *sumfile
		     ,const char *deffile
		     ,const char *ctrlfile
		     ,const char *resfile
		     );

#endif
