/*!
  \file avs_output_p.h
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief Output routines for particle data in AVS/Express format (header file)
 */
#ifndef AVS_OUTPUT_P_H
#define AVS_OUTPUT_P_H

#include "avs_output.h"

void Init_avs_p(const AVS_parameters &Avs_parameters);
void Output_avs_p(AVS_parameters &Avs_parameters
		  ,Particle *p
		  ,const CTime &time
		  );


#endif
