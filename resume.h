/*!
  \file resume.h
  \brief Routines to read/write restart file (header file)
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
 */

#ifndef RESUME_H
#define RESUME_H

#include "avs_output.h"

/*!
  \brief Save system parameters needed to restart a simulation
 */
void Save_Restart_udf(
		      double **zeta
		      ,double *uk_dc
		      ,const Particle *p
		      ,const CTime &time
		      ,double **conc_k
		      );

/*!
  \brief Set system parameters from restart file
 */
void Force_restore_parameters(
			      double **zeta
			      ,double *uk_dc
			      ,Particle *p
			      ,CTime &time
			      ,double **conc_k
			      );


#endif
