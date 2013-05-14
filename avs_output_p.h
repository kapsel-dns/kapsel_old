//
// $Id: avs_output_p.h,v 1.1 2005/07/20 16:55:01 nakayama Exp $
//
#ifndef AVS_OUTPUT_P_H
#define AVS_OUTPUT_P_H

#include "avs_output.h"

void Init_avs_p(const AVS_parameters &Avs_parameters);
void Output_avs_p(AVS_parameters &Avs_parameters
		  ,Particle *p
		  ,const CTime &time
		  );


#endif
