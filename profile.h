/*!
  \file profile.h
  \brief Smooth particle profile routines (header file)
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
 */
#ifndef PROFILE_H
#define PROFILE_H

#include <math.h>
#include "input.h"
#include "macro.h"
#include "alloc.h"
#include "interaction.h"
#include "make_phi.h"

/*!
  \details
  \f[
  h(x) = 
  \begin{cases}
  \exp{\left(-\Delta^2 / x^2\right)} & x \ge 0 \\
  0 & x < 0
  \end{cases}
  \f]
  \f$x\f$ is the radial distance from the particle center and \f$\Delta\f$ is the grid spacing
 */
inline double H(const double x){
    return x>0?exp(-SQ(DX/x)):0;
    //    return x>0?expm1(-SQ(DX/x))+1:0;
}

/*!
  \brief Main smooth profile function used in the code
  \details
  \f[
  \phi(x) = \frac{h\left(\left[a + \zeta /2\right] -x\right)}
  {h\left[\left(a + \zeta/2\right) - x\right] 
  +h\left[x - \left(a - \zeta/2\right)\right]}
  \f]
  \f$a\f$ is the particle radius and \f$\zeta\f$ the interfatial thickness.
  \param[in] x radial distance from particle center
  \param[in] radius (optional) particle radius
 */
inline double Phi(const double &x, const double radius = RADIUS){
    double dmy = H(radius+HXI - x);
    return dmy/( dmy + H(x-radius+HXI) );
}

/*!
  \brief Gaussian smooth profile function
  \details
  \f[
  \phi_e(x) = \exp{\left(\frac{x^2}{a^2}\right)}
  \f]
  \param[in] x radial distance from particle center
  \param[in] radius (optional) particle radius
 */
inline double Phi_exponential(const double &x, const double radius = RADIUS){
  double dmy = ABS(x/radius);
  //dmy = SQ(dmy)*dmy;
  dmy = SQ(dmy);
  return exp(-(SQ(dmy)));
}

/*!
  \brief Smooth profile function based on hyperbolic tangent function
  \details
  \f[
  \phi_h(x) = \frac{1}{2}\left[ 1 + \tanh\frac{a - x}{\zeta}\right]
  \f]
  \param[in] x radial distance from particle center
  \param[in] radius (optional) particle radius
 */
inline double Phi_tanh(const double &x, const double radius = RADIUS){
    return 0.5*(tanh( ((radius-x)/XI) )+1.0);
}
/*!
  \brief Derivative of the hyperbolic tangent smooth profile function \f$\phi_h\f$
 */
inline double DPhi_tanh(const double &x
			,const double radius = RADIUS
			,const double xi = XI){
    double dmy = Phi_tanh(x, radius);
    return 2.0/xi * dmy * (dmy -1.0);
}

/*!
  \brief Derivative of the smooth profile function \f$\phi(x)\f$
  \details
  \f{eqnarray*}{
  |\nabla\phi(x)| &=& 
  \frac{h'(a + \zeta/2 -x)h(x - a + \zeta/2) + h'(x - a + \zeta/2) h(a + \zeta/2 -r)}{\left[h(a + \zeta/2 - x) + h(x - a + \zeta/2)\right]^2} \\
  &=& \left[\frac{2\Delta^2}{(\zeta/2 + (a-x))^3} h(\zeta/2 + (a-x))h(\zeta/2 - (a-x))
  + \frac{2\Delta^2}{(\zeta/2 - (a - x))^3} h(\zeta/2 - (a - x)) h(\zeta/2 + (a-x))\right]\\
  &&\times
  \left[h(\zeta/2 + (a-x)) + h(\zeta/2 - (a-x))\right]^{-2}
  \f}
 */
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
	double dmy3 = DX2_2 * dmy11 / (dmy1 * (hxi + dmy_x));
	double dmy4 = DX2_2 * dmy21 / (dmy2 * (hxi - dmy_x));
	double dmy5 = dmy11 + dmy21;
	return (dmy3*dmy21 + dmy4 * dmy11)/(SQ(dmy5));
    }else{
	return 0.;
    }
}

/*!
  \brief Compact (sine) smooth profile function which exactly separates the three domains
  \details
  \f[
  \phi_s(x) = 
  \begin{cases}
  1 & x < a - \zeta /2 \\
  \frac{1}{2}\left[1 + \sin\left(\pi\frac{a - x}{\zeta}\right)\right] & |a - x| < \zeta/2 \\
  0 & x > a + \zeta/2  
  \end{cases}
  \f]
  \param[in] x radial distance from particle center
  \param[in] radius (optional) particle radius
  \note The second derivative is discontinuous at the fluid/interface boundary 
  (\f$x = a + \zeta/2)\f$
 */
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

/*!
  \brief Absolute value of the derivative of the compact (sine) smooth profile function \f$\phi_s\f$
  \details
  \f[
  |\phi_s^\prime(x)| =
  \begin{cases}
  \frac{\pi}{2\zeta} \cos(\pi \frac{a - x}{\zeta}) & x < \zeta/2 \\
  0 & x > \zeta/2
  \end{cases}
  \f]
  \param[in] x radial distance from particle center
  \param[in] radius (optional) particle radius
  \param[in] xi (optional) interface thickness
 */
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
inline double DPhi_compact_sin_norm(const double &x,
				    const double radius = RADIUS,
				    const double xi = XI){
  const double hxi = xi * 0.5;
  const double ixi = 1.0/xi;
  double dmy_x = radius - x;
  if(fabs(dmy_x)<hxi){
    return cos(M_PI * dmy_x * ixi);
  }else{
    return 0.0;
  }
}

/*!
  \brief Returns a list of local grid points that can lie within the particle domain
  \details sekibun_cell is a list of np_domain number of grid points (i,j,k) in local coordinates (i.e with respect to the center of the grid particle) that should be considered when performing integrals over the particle domain. 
  \note The center of the sekibun_cell need not coincide with the actual particle center, which is why a distance larger than the particle radius is used when deciding if a point is in/out the particle domain.
  \param[in] profile_func profile function used to determine interface/fluid boundary
  \param[out] np_domain number of points in the sekibun_cell list
  \param[out] sekibun_cell list of local grid points to consider for the particle domain
 */
void Particle_domain(
		     double (*profile_func)(const double &x, const double radius)
		     ,int &np_domain
		     ,int** &sekibun_cell
		     );

#endif
