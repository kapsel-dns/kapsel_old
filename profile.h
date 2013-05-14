//
// $Id: profile.h,v 1.1 2006/06/27 18:41:29 nakayama Exp $
//
#ifndef PROFILE_H
#define PROFILE_H

#include <math.h>
#include "input.h"
#include "macro.h"
#include "alloc.h"
#include "interaction.h"
#include "make_phi.h"

inline double H(const double x){
    return x>0?exp(-SQ(DX/x)):0;
    //    return x>0?expm1(-SQ(DX/x))+1:0;
}
inline double Phi(const double &x, const double radius = RADIUS){
    double dmy = H(radius+HXI - x);
    return dmy/( dmy + H(x-radius+HXI) );
}
inline double Phi_exponential(const double &x, const double radius = RADIUS){
  double dmy = ABS(x/radius);
  //dmy = SQ(dmy)*dmy;
  dmy = SQ(dmy);
  return exp(-(SQ(dmy)));
}
inline double Phi_tanh(const double &x, const double radius = RADIUS){
    return 0.5*(tanh( ((radius-x)/XI) )+1.0);
}
inline double DPhi_tanh(const double &x
			,const double radius = RADIUS
			,const double xi = XI){
    double dmy = Phi_tanh(x, radius);
    return 2.0/xi * dmy * (dmy -1.0);
}
inline double DPhi_compact(const double &x
			   ,const double radius = RADIUS
			   ,const double xi = XI){
    static const double DX2 = SQ(DX);
    static const double DX2_2 = DX2 * 2.;
    const double hxi = xi*.5;
    
    double dmy_x = radius - x;
    if(fabs(dmy_x)<hxi){
	double dmy1=(hxi+dmy_x)*(hxi+dmy_x);
	double dmy2=(hxi-dmy_x)*(hxi-dmy_x);
	double dmy11 = exp(-DX2/dmy1);
	double dmy21 = exp(-DX2/dmy2);
	double dmy3 = DX2_2 * dmy11 / (POW3(dmy1)/SQ(dmy1));
	double dmy4 = DX2_2 * dmy21 / (POW3(dmy2)/SQ(dmy2));
	double dmy5 = dmy11 + dmy21;
	return (dmy3*dmy21 + dmy4 * dmy11)/(SQ(dmy5));
    }else{
	return 0.;
    }
}
inline double Phi_compact_sin(const double &x
			      ,const double radius = RADIUS
			      ){
    double dmy_x = radius - x;
    if(fabs(dmy_x)<HXI){
	return 0.5*sin(M_PI*dmy_x/XI)+0.5;
    }else{
	if(dmy_x>HXI){
	    return 1.;
	}else{
	    return 0.;
	}
    }
}
inline double DPhi_compact_sin(const double &x
			       ,const double radius = RADIUS
			       ,const double xi = XI
			       ){
    const double hxi = xi*.5;
    const double ixi = 1./xi;
    
    double dmy_x = radius - x;
    if(fabs(dmy_x)<hxi){
	return M_PI *.5 * ixi * cos(M_PI*dmy_x*ixi);
    }else{
	return 0.;
    }
}

void Particle_domain(
		     double (*profile_func)(const double &x, const double radius)
		     ,int &np_domain
		     ,int** &sekibun_cell
		     ,int &np_domain_interface
		     ,int** &sekibun_cell_interface
		     ,int &np_domain_exponential
		     ,int** &sekibun_cell_exponential
		     );

#endif
