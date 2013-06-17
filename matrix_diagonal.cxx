#include "matrix_diagonal.h"

inline void jacobi_rotation(double **A, const double &s, const double &tau,
                            const int &i, const int &j, 
                            const int &k, const int &l){
  double g = A[i][j];
  double h = A[k][l];
  A[i][j] = g - s * (h + g * tau);
  A[k][l] = h + s * (g - h * tau);
}

inline void transpose(double **A, const int &n){
  for(int i = 0; i < n-1; i++){
    for(int j = i+1; j < n; j++){
      double swap = A[i][j];
      A[i][j] = A[j][i];
      A[j][i] = swap;
    }
  }
}

inline void jacobi_ordered(double **V, double *d, const int &n){
  int k;
  double dmy;
  double dmy_v;
  for(int i = 0; i < n - 1; i++){
    k = i;
    dmy = d[k];
    for(int j = i; j < n; j++){
      if(d[j] >= dmy){
        k = j;
        dmy = d[k];
      }
    }

    if(k != i){
      d[k] = d[i];
      d[i] = dmy;

      for(int j = 0; j < n; j++){
        dmy_v = V[k][j];
        V[k][j] = V[i][j];
        V[i][j] = dmy_v;
      }
    }
  }
}
void jacobi(double **A_input, double **V, double *d, int &rot, const int &n){
  const int max_step = 50;
  const double epsilon = 1.0E-17;
  double sm, c, s, t, h, g, tau, theta, tresh;
  double **A;
  double *b;
  double *z;
  A = alloc_2d_double(n,n);
  b = alloc_1d_double(n);
  z = alloc_1d_double(n);
  
  double mval = 0.0;
  double mval_inv = 0.0;
  for(int ip = 0; ip < n; ip++){
    for(int iq = 0; iq < n; iq++){
      mval = MAX(mval, ABS(A_input[ip][iq]));
    }
  }
  mval_inv = (mval > 0.0 ? 1.0/mval : 1.0);

  //initialization
  for(int ip = 0; ip < n; ip++){
    for(int iq = 0; iq < n; iq++) {
      A[ip][iq] = A_input[ip][iq] * mval_inv;
      V[ip][iq] = 0.0;
    }
    V[ip][ip] = 1.0;
  }
  for(int ip = 0; ip < n; ip++){
    b[ip] = d[ip] = A[ip][ip];
    z[ip] = 0.0;
  }
  rot = 0;

  //iteration
  int i,j;
  for(i = 1; i <= max_step; i++){
    sm = 0.0;
    for(int ip = 0; ip < n-1; ip++){
      for(int iq = ip + 1; iq < n; iq++) sm += ABS(A[ip][iq]);
    }

    //check convergence
    //machine underflow should be set to zero
    if(zero_mp(sm)){
      transpose(V, n);
      jacobi_ordered(V, d, n);
      for(int ip = 0; ip < n; ip++){
        d[ip] *= mval;
      }
      free_2d_double(A);
      free_1d_double(b);
      free_1d_double(z);
      return;
    }

    //set threshold
    if(i < 4){ 
      tresh = 0.2 * sm / SQ(n);
    }else{
      tresh = 0.0;
    }

    //jacobi sweep
    for(int ip = 0; ip < n-1; ip++){
      for(int iq = ip+1; iq < n; iq++){
        g = 100.0 * ABS(A[ip][iq]);

        //skip rotation if off-diagonal elements are negligible
        if( i > 4 && g <= epsilon * ABS(d[ip]) &&
            g <= epsilon * ABS(d[iq])){
          A[ip][iq] = 0.0;
        }else if( ABS(A[ip][iq]) > tresh){
          h = d[iq] - d[ip];

          //tangent of rotation angle
          if(g <= epsilon * ABS(h)){
            t = A[ip][iq] / h;
          }else{
            theta = 0.5 * h / A[ip][iq];
            t = 1.0 / (ABS(theta) + sqrt(1.0 + SQ(theta)));
            if(theta < 0.0) t = -t;
          }

          c = 1.0/sqrt(1 + SQ(t));
          s = t*c;
          tau = s / (1.0 + c);
          h = t * A[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          A[ip][iq] = 0.0;
          for(j = 0; j < ip; j++){
            jacobi_rotation(A, s, tau, j, ip, j, iq);
          }
          for(j = ip+1; j < iq; j++){
            jacobi_rotation(A, s, tau, ip, j, j, iq);
          }
          for(j = iq+1; j < n; j++){
            jacobi_rotation(A, s, tau, ip, j, iq, j);
          }
          for(j = 0; j < n; j++){
            jacobi_rotation(V, s, tau, j, ip, j, iq);
          }
          rot++;
        }

      }//iq
    }//ip

    for(int ip = 0; ip < n; ip++){
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }// i iterations

  free_2d_double(A);
  free_1d_double(b);
  free_1d_double(z);
  fprintf(stderr, "# Jacobi diagonalization did not converge\n");
  exit_job(EXIT_FAILURE);

}
