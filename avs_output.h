/*!
  \file avs_output.h
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief Output routines for field data in AVS/Express format (header file)
 */
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

/*!
  \brief Write AVS .fld definition files
 */
void Init_avs(const AVS_parameters &Avs_parameters);

/*!
  \brief Set global AVS parameters
 */
void Set_avs_parameters(AVS_parameters &Avs_parameters);

/*!
  \brief Output field data for current configuration in AVS format
 */
void Output_avs(AVS_parameters &Avs_parameters, double **zeta, double *uk_dc, Particle *p, const CTime &time);

/*!
  \brief Output particle data for current configuration in UDF format
 */
void Output_udf(UDFManager *ufout
		,AVS_parameters &Avs_parameters
		,double **zeta
		,double *uk_dc
		,const Particle *p
		,const CTime &time);

/*!
  \brief Output field data of charged system for current configuration
  in AVS format
 */
void Output_avs_charge(AVS_parameters &Avs_parameters, double **zeta, double *uk_dc, double **Concentration, Particle *p, const CTime &time);

#endif
