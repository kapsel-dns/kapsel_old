//
// $Id: operate_qij.h,v 1.6 2005/05/27 08:55:17 nakayama Exp $
//
#ifndef OPERATE_QIJ_H_
#define OPERATE_QIJ_H_

#include "variable.h"
#include "operate_omega.h"

extern Value *Tensor_order;

void Mem_alloc_qij(void);
void Q2q_k(Value q[QDIM]);
void Q_k2q(Value q[QDIM]);
void Q_k2q_out(const Value q[QDIM], Value q_out[QDIM]);
void Q_k2h_linear_k(const Value q[QDIM], Value h[QDIM]);

inline void Q2q_full(const Value *q, const int i,const int j, const int k, double q_full[][DIM]){
    q_full[0][0] = q[0][i][j][k]; 
    q_full[0][1] = q[1][i][j][k]; 
    q_full[0][2] = q[2][i][j][k]; 
    q_full[1][1] = q[3][i][j][k]; 
    q_full[1][2] = q[4][i][j][k]; 
    
    q_full[1][0] = q_full[0][1];
    q_full[2][0] = q_full[0][2];
    q_full[2][1] = q_full[1][2];
    q_full[2][2] = -q_full[0][0]-q_full[1][1];
}
inline double Q2qq(const double q_full[DIM][DIM]){
    double dmy=0.0;
    for(int m=0;m<DIM;m++){
	dmy += SQ(q_full[m][m]);
	for(int n=m+1;n<DIM;n++){
	    dmy += 2.*SQ(q_full[m][n]);
	}
    }
    return dmy; 
}
inline void AdotB(const double a_full[DIM][DIM], const double b_full[DIM][DIM], double c_full[DIM][DIM]){
    for(int m=0;m<DIM;m++){
	for(int n=0;n<DIM;n++){
	    c_full[m][n] = 0.0;
	    for(int l=0;l<DIM;l++){
		c_full[m][n] += a_full[m][l] * b_full[l][n];
	    }
	}
    }
}
inline void Director2Q(const double &s, const double &polar, const double &azimuthal, Value *q, const int &i, const int &j, const int &k){
    
    assert(polar >= 0.0);
    assert(polar <= M_PI);
    assert(azimuthal >= 0.0);
    assert(azimuthal <= PI2);
    double sin_polar=sin(polar);
    double cos_polar=cos(polar);
    double sin_azimuthal=sin(azimuthal);
    double cos_azimuthal=cos(azimuthal);
    
    double director[DIM];
    director[0] = sin_polar * cos_azimuthal;
    director[1] = sin_polar * sin_azimuthal;
    director[2] = cos_polar;

    double dmy_s = s;// * Q2S_prefactor;
    q[0][i][j][k] = dmy_s * ( SQ(director[0]) - One_third );
    q[1][i][j][k] = dmy_s * ( director[0] * director[1] );
    q[2][i][j][k] = dmy_s * ( director[0] * director[2] );
    q[3][i][j][k] = dmy_s * ( SQ(director[1])-One_third );
    q[4][i][j][k] = dmy_s * ( director[1] * director[2] );
}

inline void Truncate_four_seconds_rule(Value &a){
    for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	    for(int k=NZ; k<NZ_; k++){
		a[i][j][k] = 0.0;
	    }
	}
    }
    {
	int i=HNX;
	for(int j=0; j<NY; j++){
	    for(int k=0; k<NZ; k++){
		a[i][j][k] = 0.0;
	    }
	}
    }
    for(int i=0; i<NX; i++){
	int j=HNY;
	for(int k=0; k<NZ; k++){
	    //if(i == HNX || j == HNY || k >= NZ){
	    a[i][j][k] = 0.0;
	}
    }
}

#endif
