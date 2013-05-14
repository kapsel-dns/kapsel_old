/*!
  \file interaction.h
  \brief Compute inter-particle interactions (header file)
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \todo detail lees edwards pbc
 */
#ifndef INTERACTION_H
#define INTERACTION_H

#include "input.h"
#include "macro.h"

/*!
  \brief Distance vector between two points under periodic Lees-Edwards
  boundary conditions
  \details Distance vector from point \f$x1\f$ to point \f$x2\f$
  \f[
  \vec{x}_{12} = \text{LEPBC}\left(\vec{x}_2 - \vec{x}_1\right)
  \f]
  with \f$r_{12} = \sqrt{\vec{x}_{12}\cdot\vec{x}_{12}}\f$
  \param[in] x1 start point
  \param[out] x2 end point
  \param[out] r12 distance
  \param[out] x12 distance vector
 */
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

/*!
  \brief Distance vector between two points under periodic boundary
  conditions
  \details Distance vector from point \f$x1\f$ to point \f$x_2\f$
  \f[
  \vec{x}_{12} = \text{PBC}\left(\vec{x}_2 - \vec{x}_1\right)
  \f]
  with \f$r_{12} = \sqrt{\vec{x}_{12}\cdot\vec{x}_{12}}\f$
  \param[in] x1 start point
  \param[in] x2 end point
  \param[out] r12 distance
  \param[out] x12 distance vector
 */
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

/*!
  \brief Distance between two points under periodic boundary conditions
 */
inline double Distance(const double *x1
		       ,const double *x2
		       ){
    double dmy = 0.0;
    double dmy_x12[DIM];
    Distance0(x1,x2,dmy,dmy_x12);
    return dmy;
}

/*!
  \brief Magnitude of the force between two particles normalized by the
  distance between them
  \details The force on particle i, due to particle j,
  \f$\vec{F}_{i,j}\f$is 
  \f[
  \vec{F}_{i,j} = \vec{r}_{ij}\left(\frac{1}{r_{ij}} \pd{V_{LJ}}{r_{ij}}\right)
  \f]
  with \f$\vec{r}_{ij} = \vec{r}_j - \vec{r}_i\f$. For a Lennard-Jones
  potential with powers \f$2n:n\f$ we have
  \f[
  \vec{F}_{i,j} = -\vec{r}_{ij}\left[\frac{4n\epsilon}{r_{ij}^2}
  \left(2\left(\frac{\sigma}{r_{ij}}\right)^{2n} - \left(\frac{\sigma}{r_{ij}}\right)^n\right)
  \right]
  \f]
  This function returns the quantity in square brackets.
  \param[in] x distance between particles
  \param[in] sigma LJ diameter
 */
inline double Lennard_Jones_f(const double &x, const double sigma){
  //    printf("%d\n",LJ_powers);
  // ! x== 0.0 の処理を省略
  double answer=0.0;
  {
    if(LJ_powers==0){//12:6
      static const double LJ_coeff1= 24. * EPSILON;
      double dmy = sigma/x;
      dmy = SQ(dmy) * SQ(dmy) * SQ(dmy);
      answer = LJ_coeff1 / SQ(x) * ( 2.0 * SQ(dmy) - dmy );
    }
    if(LJ_powers==1){//24:12
      static const double LJ_coeff1= 48. * EPSILON;
      double dmy = sigma/x;
      dmy = SQ(dmy) * SQ(dmy) * SQ(dmy) * SQ(dmy) * SQ(dmy) * SQ(dmy);
      answer = LJ_coeff1 / SQ(x) * ( 2.0 * SQ(dmy) - dmy );
    }
    if(LJ_powers==2){//36:18
      static const double LJ_coeff1= 72. * EPSILON;
      double dmy = sigma/x;
      dmy = SQ(dmy) * SQ(dmy) * SQ(dmy) * SQ(dmy) * SQ(dmy) * SQ(dmy) * SQ(dmy) * SQ(dmy) * SQ(dmy);
      answer = LJ_coeff1 / SQ(x) * ( 2.0 * SQ(dmy) - dmy );
    }
  }
  return answer;
}
#endif
