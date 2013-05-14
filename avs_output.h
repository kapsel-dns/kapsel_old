//
// $Id: avs_output.h,v 1.12 2006/05/15 09:44:04 kin Exp $
//
#ifndef AVS_OUTPUT_H
#define AVS_OUTPUT_H

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "variable.h"
#include "macro.h"
#include "input.h"
#include "operate_omega.h"
#include "make_phi.h"
#include "operate_Qian_Sheng.h"
#include "fluid_solver.h"

enum AVS_Field {
    irregular, uniform, rectilinear
};
//////////////////////////////
typedef struct AVS_parameters{
  int nx;
  int ny;
  int nz;
  char out_fld[128];
  char out_cod[128];
  char out_pfx[128];
  char fld_file[128];
  char cod_file[128];

  char data_file[128];

  char out_pfld[128];
  char out_ppfx[128];
  char pfld_file[128];
  int istart;
  int iend;
  int jstart;
  int jend;
  int kstart;
  int kend;
  int nstep;
} AVS_parameters;

extern AVS_parameters Avs_parameters;

void Init_avs(const AVS_parameters &Avs_parameters);
void Set_avs_parameters(AVS_parameters &Avs_parameters);
void Output_avs(AVS_parameters &Avs_parameters, Value *zeta, double *uk_dc, Particle *p, const CTime &time);
void Output_avs_two_fluid(AVS_parameters &Avs_parameters
			  ,Value zeta[DIM-1]
			  ,double uk_dc[DIM]
			  ,Value u_r[DIM] // working memory
			  ,Value *conc_k
			  ,Value *conc_r // working memory
			  ,Particle *p
			  ,Value &phi // working memory
			  ,const CTime &time);
void Output_udf(UDFManager *ufout
		,AVS_parameters &Avs_parameters
		,Value *zeta
		,double *uk_dc
		,const Particle *p
		,const CTime &time);
void Output_avs_QS(AVS_parameters &Avs_parameters, Value *zeta, double *uk_dc, Value *tensor_order, Particle *p, const CTime &time);
void Output_avs_charge(AVS_parameters &Avs_parameters, Value *zeta, double *uk_dc, Value *Concentration, Particle *p, const CTime &time);

#endif
