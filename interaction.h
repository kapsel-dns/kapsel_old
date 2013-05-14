//
// $Id: interaction.h,v 1.11 2005/09/08 07:00:51 nakayama Exp $
//
#ifndef INTERACTION_H
#define INTERACTION_H

#include "input.h"
#include "macro.h"

inline void Distance0_shear(const double *x1
			    ,const double *x2
			    ,double &r12
			    ,double *x12){
    double dmy = 0.0;
    int cell[DIM];

    int d;
    {
	d=1;
	x12[d] = x2[d]-x1[d];
	cell[d] = Nint(x12[d]/L_particle[d]);
	x12[d] -= cell[d] * L_particle[d];
    }
    {
	d=0;
	x12[d] = x2[d]-x1[d];
	x12[d] += cell[1] * (-Shear_strain);
	cell[d] = Nint(x12[d]/L_particle[d]);
	x12[d] -= cell[d] * L_particle[d];
    }
    {
	d=2;
	x12[d] = x2[d]-x1[d];
	cell[d] = Nint(x12[d]/L_particle[d]);
	x12[d] -= cell[d] * L_particle[d];
    }
    for(int d=0;d<DIM;d++){
	dmy += SQ( x12[d] );
    }
    r12 = sqrt(dmy);
}
inline void Distance0_shear_hydro(const double *x1
				  ,const double *x2
				  ,double &r12
				  ,double *x12){
    double dmy = 0.0;

    int d;
    {
	d=0;
	x12[d] = x2[d]-x1[d];
	x12[d] -= (double)Nint(x12[d]/L_particle[d]) * L_particle[d];
    }
    {
	d=1;
	x12[d] = x2[d]-x1[d];
    }
    {
	d=2;
	x12[d] = x2[d]-x1[d];
	x12[d] -= (double)Nint(x12[d]/L_particle[d]) * L_particle[d];
    }
    for(int d=0;d<DIM;d++){
	dmy += SQ( x12[d] );
    }
    r12 = sqrt(dmy);
}
inline void Distance0(const double *x1
		      ,const double *x2
		      ,double &r12
		      ,double *x12
		      ){
    double dmy = 0.0;

    for(int d=0;d<DIM;d++){
	x12[d] = x2[d]-x1[d];
	x12[d] -= (double)Nint(x12[d]/L_particle[d]) * L_particle[d];
	dmy += SQ( x12[d] );
    }
    r12 = sqrt(dmy);
}
inline double Distance(const double *x1
		       ,const double *x2
		       ){
    double dmy = 0.0;
    double dmy_x12[DIM];
    Distance0(x1,x2,dmy,dmy_x12);
    return dmy;
}

inline double Lennard_Jones_f(const double &x, const double sigma){
  //    printf("%d\n",LJ_powers);
  // ! x== 0.0 の処理を省略
  double answer=0.0;
  {
    if(LJ_powers==0){
      static const double LJ_coeff1= 24. * EPSILON;
      double dmy = sigma/x;
      dmy = SQ(dmy) * SQ(dmy) * SQ(dmy);
      answer = LJ_coeff1 / SQ(x) * ( 2.0 * SQ(dmy) - dmy );
    }
    if(LJ_powers==1){
      static const double LJ_coeff1= 48. * EPSILON;
      double dmy = sigma/x;
      dmy = SQ(dmy) * SQ(dmy) * SQ(dmy) * SQ(dmy) * SQ(dmy) * SQ(dmy);
      answer = LJ_coeff1 / SQ(x) * ( 2.0 * SQ(dmy) - dmy );
    }
    if(LJ_powers==2){
      static const double LJ_coeff1= 72. * EPSILON;
      double dmy = sigma/x;
      dmy = SQ(dmy) * SQ(dmy) * SQ(dmy) * SQ(dmy) * SQ(dmy) * SQ(dmy) * SQ(dmy) * SQ(dmy) * SQ(dmy);
      answer = LJ_coeff1 / SQ(x) * ( 2.0 * SQ(dmy) - dmy );
    }
  }
  return answer;
}
#endif
