//
// $Id: interaction.h,v 1.1 2006/06/27 18:41:28 nakayama Exp $
//
#ifndef INTERACTION_H
#define INTERACTION_H

#include "input.h"
#include "macro.h"

inline void Distance0_OBL(const double *x1
			  ,const double *x2
			  ,double &r12
			  ,double *x12
    ){
    double dmy = 0.0;
    
    double signY = x2[1] - x1[1];
    x12[1] = x2[1] - x1[1];
    x12[1] -= (double)Nint(x12[1]/L_particle[1]) * L_particle[1];
    signY  -= x12[1];
    int sign = (int) signY;
    if (!(sign == 0)) {
	sign = sign/abs(sign);
    }
    dmy += SQ(x12[1]);
    
    x12[0] = x2[0] - (x1[0] + (double)sign*degree_oblique*L_particle[1]);
    x12[0] -= (double)Nint(x12[0]/L_particle[0]) * L_particle[0];
    dmy += SQ( x12[0] );
    
    x12[2] = x2[2] - x1[2];
    x12[2] -= (double)Nint(x12[2]/L_particle[2]) * L_particle[2];
    dmy += SQ( x12[2] );
    
    r12 = sqrt(dmy);
}
inline int Distance0_OBL_stepover(const double *x1
				  ,const double *x2
				  ,double &r12
				  ,double *x12
    ){
    double dmy = 0.0;
    
    double signY = x2[1] - x1[1];
    x12[1] = x2[1] - x1[1];
    x12[1] -= (double)Nint(x12[1]/L_particle[1]) * L_particle[1];
    signY  -= x12[1];
    int sign = (int) signY;
    if (!(sign == 0)) {
	sign = sign/abs(sign);
    }
    dmy += SQ(x12[1]);
    
    x12[0] = x2[0] - (x1[0] + (double)sign*degree_oblique*L_particle[1]);
    x12[0] -= (double)Nint(x12[0]/L_particle[0]) * L_particle[0];
    dmy += SQ( x12[0] );
    
    x12[2] = x2[2] - x1[2];
    x12[2] -= (double)Nint(x12[2]/L_particle[2]) * L_particle[2];
    dmy += SQ( x12[2] );

    r12 = sqrt(dmy);

    return sign;
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
