#ifndef LAD3_H
#define LAD3_H
/*!
  \file lad3.h
  \brief Basic Matrix / Vector Euclidean 3D routines (not very optimized)
  \author J. Molina
  \date 2012/01/08
  \version 1.0
 */

/*!
  \brief distance between vectors a and b
 */
inline double v_rms(const double a[DIM], const double b[DIM]){
  double c[DIM];
  for(int d = 0; d < DIM; d++){
    c[d] = a[d] - b[d];
  }
  return sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
}

// Copy routines
/*!
  \brief copy vector
 */
inline void v_copy(double copy[DIM], const double original[DIM]){
  for(int i = 0; i < DIM; i++)
    copy[i] = original[i];
}
/*!
  \brief copy matrix
 */
inline void M_copy(double copy[DIM][DIM], const double original[DIM][DIM]){
  for(int i = 0; i < DIM; i++)
    for(int j = 0; j < DIM; j++)
      copy[i][j] = original[i][j];
}

/*!
  \brief Compare two matrices (user specified relative error)
 */
inline int M_cmp(const double A[DIM][DIM], const double B[DIM][DIM],
		 const double tol=LARGE_TOL_MP){
  bool mequal = true;
  for(int i = 0; i < DIM; i++){
    for(int j = 0; j < DIM; j++){
      mequal = mequal && equal_tol(A[i][j], B[i][j], tol);
    }
  }
  return (mequal ? 1 : 0);
}

/*!
  \brief Scale vector
 */
inline void v_scale(double v[DIM], const double scale){
  for(int i = 0; i < DIM; i++)
    v[i] *= scale;
}
/*!
  \brief Scale matrix
 */
inline void M_scale(double A[DIM][DIM], const double scale){
  for(int i = 0; i < DIM; i++)
    for(int j = 0; j < DIM; j++)
      A[i][j] *= scale;
}


/*! 
  \brief Vector add : c = a + alpha * b
 */
inline void v_add(double c[DIM], const double a[DIM], const double b[DIM],
		  const double alpha=1.0){
  for(int i = 0; i < DIM; i++)
    c[i] = a[i] + alpha*b[i]; 
}
/*!
  \brief Vector add: c += alpha * b
 */
inline void v_add(double c[DIM], const double b[DIM], 
		  const double alpha=1.0){
  double a[DIM];
  v_copy(a, c);
  v_add(c, a, b, alpha);
}

/*! 
  \brief Matrix add : C = A + alpha*B
 */
inline void M_add(double C[DIM][DIM], const double A[DIM][DIM], 
		  const double B[DIM][DIM], const double alpha=1.0){
  for(int i = 0; i < DIM; i++)
    for(int j = 0; j < DIM; j++)
      C[i][j] = A[i][j] + alpha * B[i][j];
}
/*!
  \brief Matrix add : C += alpha*B
 */
inline void M_add(double C[DIM][DIM], const double B[DIM][DIM],
		  const double alpha=1.0){
  double A[DIM][DIM];
  M_copy(A, C);
  M_add(C, A, B, alpha);
}

/*!
  \brief Matrix Multiply : C = alpha*(A.B)
 */
inline void M_prod(double C[DIM][DIM], 
		   const double A[DIM][DIM],
		   const double B[DIM][DIM],
		   const double alpha=1.0){
  double dmy;
  for(int i = 0; i < DIM; i++){
    for(int j = 0; j < DIM; j++){
      dmy = 0.0;
      for(int k = 0; k < DIM; k++){
	dmy += A[i][k]*B[k][j];
      }
      C[i][j] = alpha*dmy;
    }
  }
}
/*!
  \brief Matrix Multiply : C = alpha*(C*B)
 */
inline void M_prod(double C[DIM][DIM], const double B[DIM][DIM],
		   const double alpha=1.0){
  double A[DIM][DIM];
  M_copy(A, C);
  M_prod(C, A, B, alpha);
}

/*!
  \brief Matrix determinant
 */
inline double M_det(const double A[DIM][DIM]){
  assert(DIM == 3);
  return -A[0][2]*A[1][1]*A[2][0] + A[0][1]*A[1][2]*A[2][0] + 
    A[0][2]*A[1][0]*A[2][1] - A[0][0]*A[1][2]*A[2][1] - 
    A[0][1]*A[1][0]*A[2][2] + A[0][0]*A[1][1]*A[2][2];
}

/*!
  \brief Matrix Inverse : B = alpha*(A^(-1))
 */
inline void M_inv(double B[DIM][DIM], const double A[DIM][DIM],
		  const double alpha=1.0){
  assert(DIM == 3);
  double idetA = M_det(A);
  assert(non_zero_mp(idetA));
  idetA = alpha/idetA;

  B[0][0] = (A[1][1]*A[2][2] - A[1][2]*A[2][1])*idetA;
  B[0][1] = (A[0][2]*A[2][1] - A[0][1]*A[2][2])*idetA;
  B[0][2] = (A[0][1]*A[1][2] - A[0][2]*A[1][1])*idetA;

  B[1][0] = (A[1][2]*A[2][0] - A[1][0]*A[2][2])*idetA;
  B[1][1] = (A[0][0]*A[2][2] - A[0][2]*A[2][0])*idetA;
  B[1][2] = (A[0][2]*A[1][0] - A[0][0]*A[1][2])*idetA;

  B[2][0] = (A[1][0]*A[2][1] - A[1][1]*A[2][0])*idetA;
  B[2][1] = (A[0][1]*A[2][0] - A[0][0]*A[2][1])*idetA;
  B[2][2] = (A[0][0]*A[1][1] - A[0][1]*A[1][0])*idetA;
}

/*!
  \brief Matrix Inverse : B = alpha*(B^(-1))
 */
inline void M_inv(double B[DIM][DIM], const double alpha=1.0){
  double A[DIM][DIM];
  M_copy(A, B);
  M_inv(B, A, alpha);
}

/*!
  \brief Matrix Transpose : B = alpha*A^T
 */
inline void M_trans(double B[DIM][DIM], const double A[DIM][DIM],
		    const double alpha=1.0){
  for(int i = 0; i < DIM; i++)
    for(int j = 0; j < DIM; j++)
      B[i][j] = alpha*A[j][i];
}

/*!
  \brief Matrix Transpose : B = alpha*B^T
 */
inline void M_trans(double B[DIM][DIM], const double alpha=1.0){
  double A[DIM][DIM];
  M_copy(A, B);
  M_trans(B, A, alpha);
}

/*!
  \brief Verify that matrix is a valid rotation matrix
 */
inline void M_isValidRotation(double QR[DIM][DIM]){
  double ID[DIM][DIM] = {{1.0, 0.0, 0.0},
		      {0.0, 1.0, 0.0},
		      {0.0, 0.0, 1.0}};
  double iQR[DIM][DIM];
  double tQR[DIM][DIM];
  int det, test_det, test_inv;

  test_det = (equal_tol(M_det(QR), 1.0, LARGE_TOL_MP) ? 1 : 0);
  M_trans(tQR, QR);
  M_prod(iQR, tQR, QR);
  test_inv = M_cmp(iQR, ID, LARGE_TOL_MP);
  assert(test_det && test_inv);
}

/*! 
  \brief Vector cross product: c = a x (alpha * b)
 */
inline void v_cross(double c[DIM], 
		    const double a[DIM],
		    const double b[DIM], 
		    const double alpha=1.0){
  assert(DIM == 3);
  c[0] = alpha*(a[1]*b[2] - a[2]*b[1]);
  c[1] = alpha*(a[2]*b[0] - a[0]*b[2]);
  c[2] = alpha*(a[0]*b[1] - a[1]*b[0]);
}
/*!
 \brief Vector cross product: c = c x (alpha*b)
 */
inline void v_cross(double c[DIM],
		    const double b[DIM],
		    const double alpha=1.0){
  double a[DIM];
  v_copy(a, c);
  v_cross(c, a, b, alpha);
}

/*! 
 \brief Vector inner product
 */
inline double v_inner_prod(const double a[DIM],
		    const double b[DIM]){
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
/*!
  \brief Vector norm
 */
inline double v_norm(const double a[DIM]){
  return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

/*!
  \brief Vector outer product: (AB)_ij = alpha*(a_i*b_j)
 */
inline void v_outer_prod(double AB[DIM][DIM],
		  const double a[DIM],
		  const double b[DIM],
		  const double alpha=1.0){
  for(int i = 0; i < DIM; i++)
    for(int j = 0; j < DIM; j++)
      AB[i][j] = alpha*a[i]*b[j];
}

/*!
  \brief Matrix - Vector product : y = alpha*(A.x)
 */
inline void M_v_prod(double y[DIM],
		     const double A[DIM][DIM],
		     const double x[DIM],
		     const double alpha=1.0){
  double dmy;
  for(int i = 0; i < DIM; i++){
    dmy = 0.0;
    for(int j = 0; j < DIM; j++){
      dmy += A[i][j] * x[j];
    }
    y[i] = alpha*dmy;
  }
}

/*!
  \brief Matrix - Vector product : y = alpha*(A.y)
 */
inline void M_v_prod(double y[DIM],
		     const double A[DIM][DIM],
		     const double alpha=1.0){
  double x[DIM];
  v_copy(x, y);
  M_v_prod(y, A, x, alpha);
}

/*!
  \brief Matrix - Vector product and addition : y += alpha*(A.x)
 */
inline void M_v_prod_add(double y[DIM],
			 const double A[DIM][DIM],
			 const double x[DIM],
			 const double alpha=1.0){
  double dmy;
  for(int i = 0; i < DIM; i++){
    dmy = 0.0;
    for(int j = 0; j < DIM; j++){
      dmy += A[i][j] * x[j];
    }
    y[i] += alpha*dmy;
  }
}


/*!
  \brief Vector - Matrix product : y = alpha*(x.A) =  alpha*(Transpose(A).x)
 */
inline void v_M_prod(double y[DIM],
		     const double x[DIM],
		     const double A[DIM][DIM],
		     const double alpha=1.0){
  double dmy;  
  for(int i = 0; i < DIM; i++){
    dmy = 0.0;
    for(int j = 0; j < DIM; j++){
      dmy += x[j]*A[j][i];
    }
    y[i] = alpha*dmy;
  }
}

/*!
  \brief Vector - Matrix product : y = alpha*(y.A)
 */
inline void v_M_prod(double y[DIM],
		     const double A[DIM][DIM],
		     const double alpha=1.0){
  double x[DIM];
  v_copy(x, y);
  v_M_prod(y, x, A, alpha);
}

#endif
