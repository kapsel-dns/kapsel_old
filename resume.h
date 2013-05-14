//
// $Id: resume.h,v 1.1 2006/05/15 09:19:02 kin Exp $
//
#ifndef RESUME_H
#define RESUME_H

#include "avs_output.h"

void Save_Restart_udf(UDFManager *ufout
		      ,Value *zeta
		      ,double *uk_dc
		      ,const Particle *p
		      ,const CTime &time
		      ,Value *conc_k
		      );
void Force_restore_parameters(UDFManager *ufout
			      ,Value *zeta
			      ,double *uk_dc
			      ,Particle *p
			      ,CTime &time
			      ,Value *conc_k
			      );

#endif
