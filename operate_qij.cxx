//
// $Id: operate_qij.cxx,v 1.7 2005/05/27 08:55:54 nakayama Exp $
//

#include "operate_qij.h"

Value *Tensor_order;

void Mem_alloc_qij(void){
  Tensor_order = (Value *) malloc(sizeof(Value *) * QDIM);
  for(int d=0;d<QDIM;d++){
    Tensor_order[d] = alloc_3d_double(NX, NY, NZ_);
  }
}
void Q2q_k(Value q[QDIM]){
  for(int d=0;d<QDIM;d++){
    A2a_k(q[d]);
  }
}
void Q_k2q(Value q[QDIM]){
  for(int d=0;d<QDIM;d++){
    A_k2a(q[d]);
  }
}
void Q_k2q_out(const Value q[QDIM], Value q_out[QDIM]){
  for(int d=0;d<QDIM;d++){
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ_; k++){
	  q_out[d][i][j][k] = q[d][i][j][k]; 
	}
      }
    }
  }
  Q_k2q(q_out);
}

void Q_k2h_linear_k(const Value q[QDIM], Value h[QDIM]){
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ_; k++){
	double ks[DIM];
	ks[0] =KX_int[i][j][k] * WAVE_X;
	ks[1] =KY_int[i][j][k] * WAVE_Y;
	ks[2] =KZ_int[i][j][k] * WAVE_Z;

	double q_full[DIM][DIM];
	Q2q_full(q, i, j, k, q_full);
	double divergence[DIM]={0.,0.,0.};
	for(int n=0;n<DIM;n++){
	  for(int m=0;m<DIM;m++){
	    divergence[n] +=  ks[m]*q_full[n][m];
	  }
	}
	double l2_trace=0.0;
	{
	  for(int m=0;m<DIM;m++){
	    l2_trace += ks[m] * divergence[m];
	  }
	  l2_trace *= One_third; 
	}
	  

	double dmy1 = - (L1tilde * K2[i][j][k]);

	h[0][i][j][k] = dmy1 * q_full[0][0]
	  -L2tilde * (ks[0]*divergence[0]- l2_trace);
	h[1][i][j][k] = dmy1 * q_full[0][1]
	  -L2tilde * (ks[0]*divergence[1]+ks[1]*divergence[0])*.5;
	h[2][i][j][k] = dmy1 * q_full[0][2]
	  -L2tilde * (ks[0]*divergence[2]+ks[2]*divergence[0])*.5;
	h[3][i][j][k] = dmy1 * q_full[1][1]
	  -L2tilde * (ks[1]*divergence[1]- l2_trace);
	h[4][i][j][k] = dmy1 * q_full[1][2]
	  -L2tilde * (ks[1]*divergence[2]+ks[2]*divergence[1])*.5;
      }
    }
  }
}
