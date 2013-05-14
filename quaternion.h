/*!
  \file quaternion.h
  \brief Implements simple quaternion algebra
  \author J. Molina
  \date 2012/03/28
  \version 1.0
  \see \ref page_rigid_body for more details
*/
#ifndef QUATERNION_H
#define QUATERNION_H

#include <assert.h>
#include <math.h>
#include <float.h>
#include "macro.h"
#include "parameter_define.h"

const static double QTOL=TOL_MP;
const double QTOL2 = sqrt(QTOL);
const double QTOL4 = sqrt(QTOL2);
const double QTOL6 = pow(QTOL2, 1.0/3.0);

typedef struct quaternion {
  double s; //scalar part
  double v[DIM]; //vector part
} quaternion;

/*!
  \brief Get q0 component of quaternion
 */
inline double qtn_q0(const quaternion &q){
  return q.s;
}
/*!
  \brief Get q1 component of quaternion
 */
inline double qtn_q1(const quaternion &q){
  return q.v[0];
}
/*!
  \brief Get q2 component of quaternion
 */
inline double qtn_q2(const quaternion &q){
  return q.v[1];
}
/*!
  \brief Get q3 component of quaternion
 */
inline double qtn_q3(const quaternion &q){
  return q.v[2];
}
/*!
  \brief Get scalar component of quaternion
 */
inline void qtn_scalar(double &a, const quaternion &q){
  a = q.s;
}
/*!
  \brief Get vector component of quaternion
 */
inline void qtn_vector(double v[DIM], const quaternion &q){
  v[0] = q.v[0];
  v[1] = q.v[1];
  v[2] = q.v[2];
}

/*!
  \brief Initialize by specifying all four components of quaternion
 */
inline void qtn_init(quaternion &q, const double &a0, const double &a1,
		     const double &a2, const double &a3){
  q.s = a0;
  q.v[0] = a1;
  q.v[1] = a2;
  q.v[2] = a3;
}

/*!
  \brief Initialize using 4D vector
 */
inline void qtn_init(quaternion &q, const double v4[4]){
  q.s = v4[0];
  q.v[0] = v4[1];
  q.v[1] = v4[2];
  q.v[2] = v4[3];
}

/*!
 \brief Initialize by specifying scalar, vector parts
 */
inline void qtn_init(quaternion &q, const double &si, const double vi[DIM]){
  q.s = si;
  q.v[0] = vi[0];
  q.v[1] = vi[1];
  q.v[2] = vi[2];
}

/*!
  \brief Initialize by copying qa quaternion
 */
inline void qtn_init(quaternion &q, const quaternion &qa){
  q.s = qa.s;
  q.v[0] = qa.v[0];
  q.v[1] = qa.v[1];
  q.v[2] = qa.v[2];
}

/*!
  \brief Scale quaternion: q = alpha*qa
 */
inline void qtn_scale(quaternion &q, const quaternion &qa, 
		      const double alpha){
  q.s = qa.s*alpha;
  q.v[0] = qa.v[0]*alpha;
  q.v[1] = qa.v[1]*alpha;
  q.v[2] = qa.v[2]*alpha;
}

/*!
  \brief Scale quaternion in place: q = alpha*q
 */
inline void qtn_scale(quaternion &q, const double alpha){
  q.s *= alpha;
  q.v[0] *= alpha;
  q.v[1] *= alpha;
  q.v[2] *= alpha;
}

/*!
  \brief Add two quaternions: q = qa + alpha*qb
 */
inline void qtn_add(quaternion &q, const quaternion &qa, 
		    const quaternion &qb, const double alpha=1.0){
  q.s = qa.s + qb.s*alpha;
  q.v[0] = qa.v[0] + qb.v[0]*alpha;
  q.v[1] = qa.v[1] + qb.v[1]*alpha;
  q.v[2] = qa.v[2] + qb.v[2]*alpha;
}

/*!
  \brief Add two quaternions in place: q = q + alpha*qb
 */
inline void qtn_add(quaternion &q, const quaternion &qb, 
		    const double alpha=1.0){
  q.s += qb.s*alpha;
  q.v[0] += qb.v[0]*alpha;
  q.v[1] += qb.v[1]*alpha;
  q.v[2] += qb.v[2]*alpha;
}

/*!
  \brief Quaternion conjugate: q = qa*
 */
inline void qtn_conj(quaternion &q, const quaternion &qa){
  q.s = qa.s;
  q.v[0] = -qa.v[0];
  q.v[1] = -qa.v[1];
  q.v[2] = -qa.v[2];
}

/*!
  \brief Quaternion conjugate in place: q = q*
 */
inline void qtn_conj(quaternion &q){
  q.v[0] = -q.v[0];
  q.v[1] = -q.v[1];
  q.v[2] = -q.v[2];
}

/*!
  \brief Multiply two quaternions: q = qa.(alpha*qb)
 */
inline void qtn_prod(quaternion &q, const quaternion &qa, 
		     const quaternion &qb,const double alpha=1.0){

  //scalar part
  q.s = alpha*(qa.s*qb.s - (qa.v[0]*qb.v[0] + qa.v[1]*qb.v[1] + qa.v[2]*qb.v[2]));

  //vector part
  q.v[0] = alpha*(qa.s*qb.v[0] + qb.s*qa.v[0] + qa.v[1] * qb.v[2] - qa.v[2] * qb.v[1]);
  q.v[1] = alpha*(qa.s*qb.v[1] + qb.s*qa.v[1] + qa.v[2] * qb.v[0] - qa.v[0] * qb.v[2]);
  q.v[2] = alpha*(qa.s*qb.v[2] + qb.s*qa.v[2] + qa.v[0] * qb.v[1] - qa.v[1] * qb.v[0]);
}

/*!
  \brief Multiply two quaternion in place: q = q.(alpha*qb)
 */
inline void qtn_prod(quaternion &q, const quaternion &qb, 
		     const double alpha=1.0){
  quaternion qa;
  qtn_init(qa, q);
  qtn_prod(q, qa, qb, alpha);
}

/*!
  \brief Compute quaternion norm
 */
inline double qtn_norm(const quaternion &q){
  return sqrt(q.s*q.s + q.v[0]*q.v[0] + q.v[1]*q.v[1] + q.v[2]*q.v[2]);
}

/*!
  \brief Compute quaternion square norm
 */
inline double qtn_sqnorm(const quaternion &q){
  return q.s*q.s + q.v[0]*q.v[0] + q.v[1]*q.v[1] + q.v[2]*q.v[2];
}

/*!
  \brief Compute norm of scalar part
 */
inline double qtn_norm_s(const quaternion &q){
  return ABS(q.s);
}

/*!
  \brief Compute norm of imaginary part
 */
inline double qtn_norm_v(const quaternion &q){
  return sqrt(q.v[0]*q.v[0] + q.v[1]*q.v[1] + q.v[2]*q.v[2]);
}

/*!
  \brief Check normal quaternion
 */
inline void qtn_isnormal(const quaternion &q, const double &rtol=LARGE_TOL_MP){
  if(!equal_tol(qtn_norm(q), 1.0, rtol)){
    fprintf(stderr, "#Error: unnormalized rotation quaternion %12.10E\n", 
	    qtn_norm(q));
    assert(equal_tol(qtn_norm(q), 1.0, rtol));
  }
}

/*!
  \brief Normalize quaternion
 */
inline void qtn_normalize(quaternion &q){
  double qnorm = qtn_norm(q);
  if(zero_mp(qnorm)){
    fprintf(stderr, "Erorr: qnorm = %.10f\n", qnorm);
    exit_job(EXIT_FAILURE);
  }
  qnorm = 1.0 / qnorm;
  qtn_scale(q, qnorm);
}

/*!
  \brief Compute quaternion inverse: q = alpha*qa^(-1)
 */
inline void qtn_inv(quaternion &q, const quaternion qa, const double alpha=1.0){
  double q2i = qtn_sqnorm(qa);
  assert(non_zero_mp(q2i));

  q2i = 1.0 / q2i;
  qtn_conj(q, qa);
  qtn_scale(q, q2i*alpha);
}

/*!
  \brief Compute quaternion inverse in place: q = alpha*q^(-1)
 */
inline void qtn_inv(quaternion &q, const double alpha=1.0){
  double q2i = qtn_sqnorm(q);
  assert(non_zero_mp(q2i));

  q2i = 1.0 / q2i;
  qtn_conj(q);
  qtn_scale(q, q2i*alpha);
}

/*!
  \brief Compare two quaternions
 */
inline int qtn_cmp(const quaternion &qa, const quaternion &qb, 
		    const double tol=LARGE_TOL_MP){
  bool equalq = equal_tol(qa.s, qb.s, tol) && 
    equal_tol(qa.v[0], qb.v[0], tol) &&
    equal_tol(qa.v[1], qb.v[1], tol) &&
    equal_tol(qa.v[2], qb.v[2], tol);
  return (equalq ? 1 : 0);
}

/*!
  \brief Return sinc(x)=sin(x)/x function\n
  sinc(x) = 1 - x^2 / 6 + x^4 / 120 - x^6 / 5040 + O(x^8)
 */
inline double sinc(const double &x){
  double sincx;
  if(ABS(x) > QTOL6){
    sincx = sin(x)/x;
  }else{
    double x2 = x*x;
    sincx = 1.0;
    if(ABS(x) > QTOL){
      sincx -= x2/6.0;
      if(ABS(x) > QTOL2){
	sincx += x2*x2/120.0;
	if(ABS(x) > QTOL4){
	  sincx -= x2*x2*x2/5040.0;
	}
      }
    }
  }
  return sincx;
}

/*!
  \brief Compute rotation quaternion given angle and (normal) axis vector
 */
inline void rv_rqtn(quaternion &q, const double &phi, const double n[DIM]){
  
  double dn = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if(!equal_tol(dn, 1.0, HUGE_TOL_MP)){
    fprintf(stderr, "Error: n not normal vector, |n|= %.15f\n", dn);
    exit_job(EXIT_FAILURE);
  }

  //use sinc function for increased accuracy at small angles
  double phi_half = phi/2.0;
  double sinc_phi_half = sinc(phi_half)*phi_half;
  q.s = cos(phi_half);
  q.v[0] = sinc_phi_half * n[0];
  q.v[1] = sinc_phi_half * n[1];
  q.v[2] = sinc_phi_half * n[2];
  qtn_normalize(q);
}


/*!
  \brief Compute rotation quaternion given rotation vector
  \details Rotation angle is given directly by the magnitude of the vector
 */
inline void rv_rqtn(quaternion &q, const double v[DIM]){
  //use sinc function for increased accuracy at small angles
  double phi_half = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])/2.0;
  double sinc_phi_half = sinc(phi_half)/2.0; //v=n*phi
  q.s = cos(phi_half);
  q.v[0] = sinc_phi_half * v[0];
  q.v[1] = sinc_phi_half * v[1];
  q.v[2] = sinc_phi_half * v[2];
  qtn_normalize(q);
}

/*!
  \brief Compute rotation angle and (normal) axis vector from rotation quaternion
 */
inline void rqtn_rv(double &phi, double n[DIM], const quaternion &q){

  qtn_isnormal(q);
  qtn_vector(n, q);
  double ds =  sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  phi = 2.0*asin(ds);

  if(positive_mp(ds)){
    ds = 1.0/ds;
    n[0] *= ds;
    n[1] *= ds;
    n[2] *= ds;
  }else{
    phi = 0.0;
    n[0] = 0.0;
    n[1] = 0.0;
    n[2] = 1.0;
  }
}

/*!
  \brief Computes rotation vector from rotation quaternion
 */
inline void rqtn_rv(double v[DIM], const quaternion &q){
  qtn_isnormal(q);
  qtn_vector(v, q);
  double ds = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  double phi = 2.0*asin(ds);

  if(positive_mp(ds)){
    ds = phi/ds;
    v[0] *= ds;
    v[1] *= ds;
    v[2] *= ds;
  }else{
    phi = 0.0;
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 1.0;
  }
}

/*!
  \brief Computes euler angles from rotation quaternion (z-x-z / 3-1-3 convention)
 */
inline void rqtn_euler(double &psi, double &theta, double &phi,
                       const quaternion &q){
  const double &q0 = q.s;
  const double &q1 = q.v[0];
  const double &q2 = q.v[1];
  const double &q3 = q.v[2];
  double dmy_cos;

  dmy_cos = (1.0 - 2.0*(q2*q2 + q1*q1));

  if(!equal_tol(ABS(dmy_cos), 1.0, TOL_MP)){ // check for gimbal lock
    psi = atan2(2.0*(q1*q3 - q0*q2), 2.0*(q0*q1 + q2*q3));
    theta = acos(dmy_cos);
    phi = atan2(2.0*(q1*q3 + q0*q2), 2.0*(q0*q1 - q2*q3));
  }else {//rotation with respect to z-axis only, set phi=0!
    if(dmy_cos > 0){ //theta = 0
      psi = atan2(2.0*(q0*q3 - q1*q2), 1.0 - 2.0 * (q2*q2 + q3*q3));
      theta = 0.0;
      phi = 0.0;
    }else{ //theta = pi
      psi = atan2(2.0*(q1*q2 - q0*q3), 1.0 - 2.0 * (q2*q2 + q3*q3));
      theta = M_PI;
      phi = 0.0;
    }
  }
}

/*!
  \brief Compute random rotation quaternion
 */
inline void random_rqtn(quaternion &q){
  double q0,q1,q2,q3;
  RA_on_sphere4D(q0,q1,q2,q3);
  qtn_init(q, q0, q1, q2, q3);
  qtn_isnormal(q);
}

/*!
  \brief Compute rotation matrix given rotation quaternion
 */
inline void rqtn_rm(double R[DIM][DIM], const quaternion &q){
  const double &q0 = q.s;
  const double &q1 = q.v[0];
  const double &q2 = q.v[1];
  const double &q3 = q.v[2];

  R[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
  R[0][1] = 2.0*(q1*q2 - q0*q3);
  R[0][2] = 2.0*(q1*q3 + q0*q2);
  
  R[1][0] = 2.0*(q1*q2 + q0*q3);
  R[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
  R[1][2] = 2.0*(q2*q3 - q0*q1);

  R[2][0] = 2.0*(q1*q3 - q0*q2);
  R[2][1] = 2.0*(q2*q3 + q0*q1);
  R[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;
}

/*!
  \brief Compute rotation quaternion given rotation matrix
 */
inline void rm_rqtn(quaternion &q, const double R[DIM][DIM]){
  double qq[4];

  //q_i^2 components
  qq[0] = 1.0 + R[0][0] + R[1][1] + R[2][2];
  qq[1] = 1.0 + R[0][0] - R[1][1] - R[2][2];
  qq[2] = 1.0 - R[0][0] + R[1][1] - R[2][2];
  qq[3] = 1.0 - R[0][0] - R[1][1] + R[2][2];
  assert(qq[0] >= 0 && qq[1] >= 0 && qq[2] >= 0 && qq[3] >= 0);

  double dmy;
  double mm = qq[0];
  int i = 0;
  for(int d = 1; d < 4; d++){ // find component with largest absolute magnitude
    dmy = qq[d];
    if(dmy > mm){
      mm = dmy;
      i = d;
    }
  }
  switch (i){
  case 0:
    qq[0] = sqrt(qq[0]) / 2.0;
    dmy = 1.0/(4.0 * qq[0]);

    qq[1] = (R[2][1] - R[1][2]) * dmy;
    qq[2] = (R[0][2] - R[2][0]) * dmy;
    qq[3] = (R[1][0] - R[0][1]) * dmy;
    break;
  case 1:
    qq[1] = sqrt(qq[1]) / 2.0;
    dmy = 1.0/(4.0 * qq[1]);

    qq[0] = (R[2][1] - R[1][2]) * dmy;
    qq[2] = (R[1][0] + R[0][1]) * dmy;
    qq[3] = (R[0][2] + R[2][0]) * dmy;
    break;
  case 2:
    qq[2] = sqrt(qq[2]) / 2.0;
    dmy = 1.0/(4.0 * qq[2]);

    qq[0] = (R[0][2] - R[2][0]) * dmy;
    qq[1] = (R[1][0] + R[0][1]) * dmy;
    qq[3] = (R[2][1] + R[1][2]) * dmy;
    break;
  case 3:
    qq[3] = sqrt(qq[3]) / 2.0;
    dmy = 1.0/(4.0 * qq[3]);

    qq[0] = (R[1][0] - R[0][1]) * dmy;
    qq[1] = (R[2][0] + R[0][2]) * dmy;
    qq[2] = (R[2][1] + R[1][2]) * dmy;
    break;
  }
  qtn_init(q, qq);
}

#endif
