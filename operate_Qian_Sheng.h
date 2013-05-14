//
// $Id: operate_Qian_Sheng.h,v 1.15 2005/06/27 10:18:46 nakayama Exp $
//
#ifndef OPERATE_QIAN_SHENG_H_
#define OPERATE_QIAN_SHENG_H_

#include <float.h>
#include "variable.h"
#include "operate_qij.h"
#include "particle_solver.h"
#include "f_particle.h"

extern Value *Stress;
extern Value *Q_rhs;
extern Value *tmp1_QDIM;
extern Value *tmp2_QDIM;

void Init_tensor_order(Value *q, Particle *p);
void Mem_alloc_QS(void);
void QS_rhs(Value zeta[DIM-1], const double uk_dc[DIM], Value q_k_orig[QDIM], Value stress[DIM*DIM], Value q_rhs[QDIM], const Index_range *ijk_range, const int &n_ijk_range);
void AC_rhs(Value Q_k[QDIM],Value q_rhs[QDIM],Particle *p);
void OG_rhs(Value zeta[DIM-1], const double uk_dc[DIM], Value q_k_orig[QDIM], Value stress[DIM*DIM], Value q_rhs[QDIM], const Index_range *ijk_range, const int &n_ijk_range);
void Advection_slippy_NS(Value zeta[DIM-1], const double uk_dc[DIM], Value rhs[DIM-1], const Index_range *ijk_range, const int &n_ijk_range, Particle *p);
void Rhs_slippy_NS_stress(Value zeta[DIM-1], const double uk_dc[DIM], Value rhs[DIM-1], const Index_range *ijk_range, const int &n_ijk_range, Particle *p);
void Rhs_slippy_NS_surface_normal(Value zeta[DIM-1], double force_dc[DIM], Value rhs[DIM-1], const Index_range *ijk_range, const int &n_ijk_range, Particle *p);

inline void AC_rhs_test(Value Q_k[QDIM], Value q_rhs[QDIM]){
    
    ///// Q(r)
    {
	for(int d=0;d<QDIM;d++){
	    Truncate_four_seconds_rule(Q_k[d]);
	    for(int i=0; i<NX; i++){
		for(int j=0; j<NY; j++){
		    for(int k=0; k<NZ_; k++){
			tmp1_QDIM[d][i][j][k]=Q_k[d][i][j][k];
		    }
		}
	    }
	}
	///////////////////
	for(int d=0;d<QDIM;d++){
	    K_space_extension(tmp1_QDIM[d]);
	    // test for result of extension
	    if(0){
		for(int i=0; i<N2X; i++){
		    for(int j=0; j<N2Y; j++){
			for(int k=0; k<N2Z_; k++){
			    if(ABS(tmp1_QDIM[d][i][j][k])>0){
				fprintf(stderr, "(%d %d %d) %g\n"
					,i,j,k
					,tmp1_QDIM[d][i][j][k]);
				
			    }
			    int kxe = Calc_KX_extension_Ooura(i,j,k);
			    int kye = Calc_KY_extension_Ooura(i,j,k);
			    int kze = Calc_KZ_extension_Ooura(i,j,k);
			    if((kxe>=HNX) || kxe <= -HNX 
			       || (kye>=HNY) || (kye <= -HNY)
			       ||(kze>=HNZ) || (kze <= -HNZ)){
				if(tmp1_QDIM[d][i][j][k] != 0.0){
				    fprintf(stderr
					    ,"a: %d: %d %d %d: %g\n"
					    ,d,i,j,k
					    ,tmp1_QDIM[d][i][j][k]);
				}
			    }else{
				assert(k<NZ);
				assert( (i< HNX) || (i> N2X-HNX));
				assert( (j< HNY) || (j> N2Y-HNY));
				int i2=i,j2=j,k2=k;
				if(i > HNX){
				    i2 = i -NX;
				}
				if(j > HNY){
				    j2 = j- NY;
				}
				assert(Calc_KX_extension_Ooura(i,j,k) == Calc_KX(i2,j2,k2));
				assert(Calc_KY_extension_Ooura(i,j,k) == Calc_KY(i2,j2,k2));
				assert(Calc_KZ_extension_Ooura(i,j,k) == Calc_KZ(i2,j2,k2));
				if(tmp1_QDIM[d][i][j][k]!=Q_k[d][i2][j2][k2]){
				    //if(fabs(tmp1_QDIM[d][i][j][k]-Q_k[d][i2][j2][k2])> 1.e-5){
				    fprintf(stderr
					    ,"b: %d: %d %d %d: %g %g\n"
					    ,d,i,j,k
					    ,tmp1_QDIM[d][i][j][k]
					    ,Q_k[d][i][j][k]);
				};
			    }
			}
		    }
		    fprintf(stderr, "\n");
		}
	    }
	    A_k2a_extended(tmp1_QDIM[d]);
	}
	
      cerr << "AAAAAA\n" << endl;
      for(int d=0;d<QDIM;d++){
	      A2a_k_extended(tmp1_QDIM[d]);
	  if(0){
	      for(int i=0; i<N2X; i++){
		  for(int j=0; j<N2Y; j++){
		      for(int k=0; k<N2Z_; k++){
			  if(ABS(tmp1_QDIM[d][i][j][k])>1.e-10){
			      fprintf(stderr, "(%d %d %d) %g\n"
				      ,i,j,k
				      ,tmp1_QDIM[d][i][j][k]);
			  }
		      }
		  }
		  fprintf(stderr, "\n");
	      }
	  }
	      K_space_extension_reverse(tmp1_QDIM[d]);
      }

      cerr << "BBBBBB\n" << endl;
      if(1){
	  double error = 0.0;
	  for(int d=0;d<QDIM;d++){
	      for(int i=0; i<NX; i++){
		  for(int j=0; j<NY; j++){
		      for(int k=0; k<NZ_; k++){
			  if(Q_k[d][i][j][k]-tmp1_QDIM[d][i][j][k] > 1.e-9){
			      fprintf(stderr, "%d (%d %d %d) %g %g\n"
				      ,d
				      ,i,j,k
				      ,Q_k[d][i][j][k]
				      ,tmp1_QDIM[d][i][j][k]);
			  }
			  
			  error = MAX(error,fabs(Q_k[d][i][j][k]-tmp1_QDIM[d][i][j][k]));
		      }
		  }
	      }
	      fprintf(stderr, "%g\n", error);
	      abort();
	  }
      }
      ///////////////////
  }
  ///// molcular field terms
  Q_k2h_linear_k(Q_k, q_rhs);
  for(int d=0;d<QDIM;d++){
      K_space_extension(q_rhs[d]);
      A_k2a_extended(q_rhs[d]);
  }

  double linear_max =-DBL_MAX;
  double linear_min =DBL_MAX;
  double qq_trace_max = 0.0;
  double qq_trace_min = DBL_MAX;
  for(int i=0; i<N2X; i++){
    for(int j=0; j<N2Y; j++){
      for(int k=0; k<N2Z; k++){
	
	double q_full[DIM][DIM];
	Q2q_full(tmp1_QDIM, i, j, k, q_full);
	double qqij[DIM][DIM];
	AdotB(q_full, q_full, qqij);
	double qq_trace= 0.;
	for(int m=0;m<DIM;m++){
	  qq_trace += qqij[m][m];
	}
	double dmy_c = -C_LdG * qq_trace;
	if(0){
	  double qq_trace2= Q2qq(q_full);
	  fprintf(stderr, "%g %g %g\n",qq_trace, qq_trace2,qq_trace-qq_trace2);
	}
	{
	  double dmy = A_LdG - dmy_c;
	  linear_min = MIN(dmy, linear_min);
	  linear_max = MAX(dmy, linear_max);
	  qq_trace_max = MAX(qq_trace_max, qq_trace); 
	  qq_trace_min = MIN(qq_trace_min, qq_trace); 
	}			 
	//fprintf(stderr, "%g %g %g\n", qq_trace, qq_trace2, dmy_c);
	qq_trace *= (One_third * B_LdG);
	
	double h_full[DIM][DIM];
	Q2q_full(q_rhs, i, j, k, h_full);
	{
	  for(int m=0;m<DIM;m++){
	    for(int n=m;n<DIM;n++){
	      h_full[m][n] += (dmy_c * q_full[m][n]+B_LdG *qqij[m][n]);
	    }
	  }
	  for(int m=0;m<DIM;m++){
	    h_full[m][m] -= qq_trace;
	  }
	  for(int m=0;m<DIM;m++){
	    for(int n=m+1;n<DIM;n++){
	      h_full[n][m] = h_full[m][n];
	    }
	  }
	}
	{
	  q_rhs[0][i][j][k] += One_overMu1 * h_full[0][0];
	  q_rhs[1][i][j][k] += One_overMu1 * h_full[0][1];
	  q_rhs[2][i][j][k] += One_overMu1 * h_full[0][2];
	  q_rhs[3][i][j][k] += One_overMu1 * h_full[1][1];
	  q_rhs[4][i][j][k] += One_overMu1 * h_full[1][2];
	}
      }
    }
  }
  fprintf(stdout, "%g %g %g %g %g %g\n", -linear_max, -linear_min
	  ,qq_trace_max
	  ,qq_trace_min
	  , sqrt(qq_trace_max*Q2S_prefactor)/A_eq
	  , sqrt(qq_trace_min*Q2S_prefactor)/A_eq
	  );
  {
    for(int d=0;d<QDIM;d++){
      A2a_k_extended(q_rhs[d]);
      K_space_extension_reverse(q_rhs[d]);
    }    
  }
}

#endif
