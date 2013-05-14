/*!
  \file rigid_body.h
  \brief Auxiliary routines to solve equations of motion for rigid
  bodies (header file)
  \author J. Molina
  \date 2012/03/29
  \version 1.0
 */
#ifndef RIGID_BODY_H
#define RIGID_BODY_H

#include "quaternion.h"
#include "lad3.h"

enum COORD_SYSTEM {BODY_FRAME, SPACE_FRAME};
enum COORD_TRANS {BODY2SPACE, SPACE2BODY};

/*!
  \brief Compute random rotation matrix
 */
inline void random_rotation(double QR[DIM][DIM]){
  quaternion dmy_q;
  random_rqtn(dmy_q);
  rqtn_rm(QR, dmy_q);
  M_isValidRotation(QR);
}


//
// Get skew anti-symmetric matrix
inline void skew(double ws[DIM][DIM], 
		 const double w[DIM]){
  assert(DIM == 3);
  ws[0][0] = 0.0;
  ws[0][1] = -w[2];
  ws[0][2] = w[1];

  ws[1][0] = w[2];
  ws[1][1] = 0.0;
  ws[1][2] = -w[0];

  ws[2][0] = -w[1];
  ws[2][1] = w[0];
  ws[2][2] = 0.0;
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////
// Vector Transformation: Space <--> Body Coordinate Frames
//////////////////////////////////////////////////
//////////////////////////////////////////////////
/*!
  \brief Transform between body/space and space/body frames given the
  current orientation QUATERNION
 */
void rigid_body_rotation(double rotated[DIM], 
				const double original[DIM], 
				const quaternion &q, 
				const COORD_TRANS &transform);

/*!
  \brief Transform between body/space and space/body frames given the
  current orientation QUATERNION (in place)
 */
inline void rigid_body_rotation(double rotated[DIM], 
				const quaternion &q, 
				const COORD_TRANS &transform){
  double original[DIM];
  v_copy(original, rotated);
  rigid_body_rotation(rotated, original, q, transform);
}


/*!
  \brief Transform between body/space and space/body frames given the 
  current orientation MATRIX
 */
void rigid_body_rotation(double rotated[DIM], 
			 const double original[DIM],
			 const double QR[DIM][DIM], 
			 const COORD_TRANS &transform);
/*!
  \brief Transform between body/space and space/body frames given the 
  current orientation MATRIX (in place)
 */
inline void rigid_body_rotation(double rotated[DIM], 
				const double QR[DIM][DIM], 
				const COORD_TRANS &transform){
  double original[DIM];
  v_copy(original, rotated);
  rigid_body_rotation(rotated, original, QR, transform);
}

//////////////////////////////////////////////////
//
// Kinematics Calculations
//

// Quaternion Representation
/*!
  \brief Compute time derivative of orientation quaternion
  \param[out] dqdt time derivative of q
  \param[in] q current orientation quaternion
  \param[in] omega current angular velocity
  \param[in] coord coordinate system of given angular velocity
 */
void qdot(quaternion &dqdt, 
	  const quaternion &q, 
	  const double omega[DIM], 
	  const COORD_SYSTEM &coord);
/*!
  \brief Compute time derivative of orientation quaternion (in place)
 */
inline void qdot(quaternion &dqdt,
		 const double omega[DIM],
		 const COORD_SYSTEM &coord){
  quaternion q;
  qtn_init(q, dqdt);
  qdot(dqdt, q, omega, coord);
}

/*!
  \brief Compute time derivative of orientation matrix
  \param[out] dQRdt time derivative of QR
  \param[in] QR current orientation matrix
  \param[in] omega current angular velocity
  \param[in] coord coordinate system for given angular velocity
 */
void Qdot(double dQRdt[DIM][DIM],
	  const double QR[DIM][DIM],
	  const double omega[DIM],
	  const COORD_SYSTEM &coord);
/*!
  \brief Compute time derivative of orientation matrix (in place)
 */
inline void Qdot(double dQRdt[DIM][DIM],
		 const double omega[DIM],
		 const COORD_SYSTEM &coord){
  double QR[DIM][DIM];
  M_copy(QR, dQRdt);
  Qdot(dQRdt, QR, omega, coord);
}

/*!
  \brief Compute time derivative of angular velocity given
  the current velocity and torque on the body, by solving the Euler
  Equations for rigid body motion
  \param[out] dwdt time derivative of angular velocity (body frame)
  \param[in] w angular velocity (body frame)
  \param[in] tau torque (body frame)
  \param[in] QR orientation matrix
  \param[in] I inertia tensor (body frame)
 */
/*void wdot(double dwdt[DIM],
	  const double w[DIM],
	  const double tau[DIM],
	  const double I[DIM][DIM]);
*/
/*!
  \brief Compute time derivative of angular velocity (in place)
 */
/*inline void wdot(double dwdt[DIM],
		 const double tau[DIM],
		 const double I[DIM][DIM]);
*/
/*!
  \brief Update orientation/ angular velocity using fourth order RK scheme
  of a free rigid body (torque=0)
  \param[in,out] q orientation quaternion
  \param[in,out] omega angular velocity (BODY FRAME)
  \param[in] I intertia tensor
  \param[in] dt time step
 */
/*void propagate_wq_RK4(quaternion &q, 
		      double omega[DIM],
		      const double I[DIM][DIM],
		      const double &dt);
*/
//Rotation Matrix - Lie Group representation
/*!
  \brief Update orientation/ angular velocity using fourth order RK scheme 
  for omega and first order Magnus series expansion for orientation matrix
  for free rigid body (torque = 0)
  \param[in,out] QR orientation matrix
  \param[in,out] omega angular velocity
  \param[in] I inertia tensor
  \param[in] dt time step
 */
/*void propagate_w_RK4_Q_Euler(double QR[DIM][DIM], 
			     double omega[DIM],
			     const double I[DIM][DIM],
			     const double &dt);
*/
#endif
