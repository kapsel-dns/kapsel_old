//
// $Id: init_particle.h,v 1.1 2006/06/27 18:41:28 nakayama Exp $
//
#ifndef INIT_PARTICLE_H
#define INIT_PARTICLE_H

#include <assert.h> 
#include "macro.h"
#include "variable.h"
#include "md_force.h"
#include "avs_output.h"
#include "input.h"
#include "fluct.h"

void Init_Particle(Particle *p);
void Init_Chain(Particle *p);
void Show_parameter(AVS_parameters Avs_parameters, Particle *p);

inline void Show_particle(Particle *p){
    for(int n=0;n<Particle_Number;n++){
	fprintf(stderr, "%g %g %g %g %g %g\n"
		,p[n].v[0]
		,p[n].v[1]
		,p[n].v[2]
		,p[n].omega[0]
		,p[n].omega[1]
		,p[n].omega[2]
		);
    }
    fprintf(stderr, "\n\n");
}
#endif

