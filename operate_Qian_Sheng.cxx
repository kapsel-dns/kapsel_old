//
// $Id: operate_Qian_Sheng.cxx,v 1.24 2005/07/20 17:39:44 nakayama Exp $
//
#include "operate_Qian_Sheng.h"

Value *Stress;
Value *Q_rhs;
Value *tmp1_QDIM;
Value *tmp2_QDIM;


void Init_pretilt(Value Qr_extended[QDIM]
		  ,const Value DPHI_extended[DIM]
		  ,const Particle *p
		  ){
  const double dx = 0.5 * DX;
  int **sekibun_cell = Sekibun_cell_extended;
  int *Nlattice = N2s;
  int np_domain = NP_domain_extended;
  
  
  for(int n=0; n < Particle_Number; n++){
    double xp[DIM];
    for(int d=0;d<DIM;d++){
      xp[d] = p[n].x[d];
    }
    const double anchor_scalar_order = S_surfaces[p[n].spec];
    
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell 
      = Particle_cell(xp, dx, x_int, residue);// {1,0} が返ってくる
    sw_in_cell = 1;
    int r_mesh[DIM];
    double r[DIM];
    for(int mesh=0; mesh < np_domain; mesh++){
      Relative_coord(sekibun_cell[mesh]
		     , x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
      double dphi[DIM];
      double amp2_dphi = 0.;
      for(int m=0;m<DIM;m++){
	dphi[m] = DPHI_extended[m][r_mesh[0]][r_mesh[1]][r_mesh[2]];
	amp2_dphi += SQ(dphi[m]);
      }
      if(amp2_dphi > 0.0){
	double amp1_dphi = sqrt(amp2_dphi);
	double iamp2_dphi = 1./amp2_dphi;

	double q_full_anchor[DIM][DIM];

	{// pretilt tensor order 
	  for(int m=0;m<DIM;m++){
	    for(int n=m;n<DIM;n++){
	      q_full_anchor[m][n] = dphi[m] * dphi[n] * iamp2_dphi;
	    }
	  }
	  for(int m=0;m<DIM;m++){
	    q_full_anchor[m][m] -= One_third;
	  }
	  for(int m=0;m<DIM;m++){
	    for(int n=m;n<DIM;n++){
	      q_full_anchor[m][n] *= anchor_scalar_order;
	    }
	  }
	}

	Qr_extended[0][r_mesh[0]][r_mesh[1]][r_mesh[2]] = q_full_anchor[0][0];
	Qr_extended[1][r_mesh[0]][r_mesh[1]][r_mesh[2]] = q_full_anchor[0][1];
	Qr_extended[2][r_mesh[0]][r_mesh[1]][r_mesh[2]] = q_full_anchor[0][2];
	Qr_extended[3][r_mesh[0]][r_mesh[1]][r_mesh[2]] = q_full_anchor[1][1];
	Qr_extended[4][r_mesh[0]][r_mesh[1]][r_mesh[2]] = q_full_anchor[1][2];
      }
    }
  }
}
void Init_tensor_order(Value *q, Particle *p){
  double sinit=0.1 * A_eq;//.8;//.1;
  SRA(GIVEN_SEED,10);
  if(1){
    for(int i=0; i<N2X; i++){
      for(int j=0; j<N2Y; j++){
	for(int k=0; k<N2Z; k++){
	  double polar = RAx(M_PI);
	  double azimuthal = RAx(PI2);
	  double sdmy= sinit;
	  
	  Director2Q(sdmy, polar, azimuthal, tmp1_QDIM, i, j, k); 
	}
      }
    }  
    if(Particle_Number> 0){
      {
	// nabla phi(r)
	Reset_phi_extended(tmp2_QDIM[0]);
	Make_phi_particle_extended(tmp2_QDIM[0],p);
	A2a_k_extended(tmp2_QDIM[0]);
	A_k2dxa_k_extended(tmp2_QDIM[0], tmp2_QDIM[1]);
	A_k2dya_k_extended(tmp2_QDIM[0], tmp2_QDIM[2]);
	A_k2dza_k_extended(tmp2_QDIM[0], tmp2_QDIM[3]);
	for(int d=0;d<DIM;d++){
	  A_k2a_extended(tmp2_QDIM[d+1]);
	}
      }
      Init_pretilt(tmp1_QDIM, tmp2_QDIM+1,p);
    }
    for(int d=0;d<QDIM;d++){
      A2a_k_extended(tmp1_QDIM[d]);
      K_space_extension_reverse(tmp1_QDIM[d]);
      if(SW_EQ == nematic_Allen_Cahn){
	Truncate_four_seconds_rule(tmp1_QDIM[d]);
      }else if(SW_EQ == Olmsted_Goldbart || SW_EQ == Qian_Sheng){
	for(int d=0;d<QDIM;d++){
	  Truncate_four_fifths_rule(q[d]);
	}
      }
      for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	  for(int k=0; k<NZ_; k++){
	    q[d][i][j][k] = tmp1_QDIM[d][i][j][k];
	  }
	}
      }
    }
  }else{
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ; k++){
	  double polar = RAx(M_PI);
	  double azimuthal = RAx(PI2);
	  double sdmy= sinit * POW3(2.0);
	  Director2Q(sdmy, polar, azimuthal, q, i, j, k); 
	}
      }
    }
    Q2q_k(q);
    if(SW_EQ == nematic_Allen_Cahn){
      for(int d=0;d<QDIM;d++){
	Truncate_four_seconds_rule(q[d]);
      }
    }else if(SW_EQ == Olmsted_Goldbart || SW_EQ == Qian_Sheng){
      for(int d=0;d<QDIM;d++){
	Truncate_four_fifths_rule(q[d]);
      }
    }
  }
}
void Mem_alloc_QS(void){
  tmp1_QDIM = (Value *) malloc(sizeof(Value *) * QDIM);
  tmp2_QDIM = (Value *) malloc(sizeof(Value *) * QDIM);
  Q_rhs = (Value *) malloc(sizeof(Value *) * QDIM);
  Stress = (Value *) malloc(sizeof(Value *) * (DIM*DIM));
  for(int d=0;d<QDIM;d++){
    tmp1_QDIM[d] = alloc_3d_double(N2X, N2Y, N2Z_);
    tmp2_QDIM[d] = alloc_3d_double(N2X, N2Y, N2Z_);
    Q_rhs[d] = alloc_3d_double(N2X, N2Y, N2Z_);
  }
  for(int d=0;d<DIM*DIM;d++){
    Stress[d] = alloc_3d_double(N2X, N2Y, N2Z_);
  }

}

inline void Stress_k2zeta_rhs_k(Value stress[DIM*DIM], const Index_range *ijk_range, const int &n_ijk_range){
  int offset=DIM-1;
  assert(offset < DIM*DIM);
  for(int n=0;n<n_ijk_range;n++){
    for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
      for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
    	for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	  double ks[DIM];
	  ks[0] =KX_int[i][j][k] * WAVE_X;
	  ks[1] =KY_int[i][j][k] * WAVE_Y;
	  ks[2] =KZ_int[i][j][k] * WAVE_Z;
	  
	  double sigma[DIM][DIM];
	  for(int m=0;m<DIM;m++){
	    for(int n=0;n<DIM;n++){
	      int l= m*DIM + n;
	      sigma[m][n] = stress[l][i][j][k];
	    }
	  }
	  double div_stress[DIM];
	  for(int n=0;n<DIM;n++){
	    for(int m=0;m<DIM;m++){
	      div_stress[n] = ks[m] * sigma[m][n];
	    }
	  }
	  stress[offset+0][i][j][k]= -IRHO*(ks[1] * div_stress[2] - ks[2] * div_stress[1]);
	  stress[offset+1][i][j][k]= -IRHO*(ks[2] * div_stress[0] - ks[0] * div_stress[2]);
	  stress[offset+2][i][j][k]= -IRHO*(ks[0] * div_stress[1] - ks[1] * div_stress[0]);
	}
      }
    }
  }
  Omega_k2zeta_k(stress+offset, stress);
}

inline double vdotnablaQ(const double v[DIM], const double dqij[][DIM][DIM], const int &i, const int&j){
  double dmy=0.0;
  for(int l=0;l<DIM;l++){
    dmy += (v[l] * dqij[l][i][j]);
  }
  return dmy;
}

inline void QS_rhs1(const Value zeta[DIM-1], const double uk_dc[DIM], const Value q_k_orig[QDIM], Value stress[DIM*DIM], Value q_rhs[QDIM]){
  {
    {
      Zeta_k2u_k(zeta, uk_dc, tmp2_QDIM); ///// making v(x)
      for(int d=0;d<DIM;d++){
	K_space_extension(tmp2_QDIM[d]);
	A_k2a_extended(tmp2_QDIM[d]);
      }
    }
    
    Value nabla_q[][QDIM] =
      {{stress[0],stress[1],stress[2],stress[3],stress[4]},
       {stress[5],stress[6],stress[7],stress[8],q_rhs[0]},
       {q_rhs[1],q_rhs[2],tmp1_QDIM[0],tmp1_QDIM[1],tmp1_QDIM[2]}};
    for(int d=0;d<QDIM;d++){///// making \nabla Q(x)
      A_k2dxa_k(q_k_orig[d],nabla_q[0][d]);
      A_k2dya_k(q_k_orig[d],nabla_q[1][d]);
      A_k2dza_k(q_k_orig[d],nabla_q[2][d]);

      for(int l=0;l<DIM;l++){
	K_space_extension(nabla_q[l][d]);
	A_k2a_extended(nabla_q[l][d]);
      }
    }
    ///// site-local operation
    for(int i=0; i<N2X; i++){
      for(int j=0; j<N2Y; j++){
	for(int k=0; k<N2Z; k++){

	  double v[DIM];
	  double dqij[DIM][DIM][DIM];
	  for(int d=0;d<DIM;d++){
	    v[d] = tmp2_QDIM[d][i][j][k];
	    Q2q_full(nabla_q[d], i,j,k, dqij[d]);
	  }

	  double dmy_q_rhs[QDIM];
	  { ///// advection of Q
	    dmy_q_rhs[0]=-vdotnablaQ(v,dqij,0,0);
	    dmy_q_rhs[1]=-vdotnablaQ(v,dqij,0,1);
	    dmy_q_rhs[2]=-vdotnablaQ(v,dqij,0,2);
	    dmy_q_rhs[3]=-vdotnablaQ(v,dqij,1,1);
	    dmy_q_rhs[4]=-vdotnablaQ(v,dqij,1,2);
	  }
	  double sigma[DIM][DIM];
	  {///// advection of v
	    sigma[0][0] = -RHO*SQ(v[0]);
	    sigma[1][1] = -RHO*SQ(v[1]);
	    sigma[2][2] = -RHO*SQ(v[2]);

	    sigma[0][1] = -RHO*v[0]*v[1];
	    sigma[0][2] = -RHO*v[0]*v[2];
	    sigma[1][2] = -RHO*v[1]*v[2];

	    sigma[1][0] = sigma[0][1];
	    sigma[2][0] = sigma[0][2];
	    sigma[2][1] = sigma[1][2];
	  }
	  {/////// distortion stress
	    double dmy_sigma[DIM][DIM]={{0.,0.,0.},
					{0.,0.,0.},
					{0.,0.,0.}};
	    assert(L1tilde >= 0.0); 
	    if(L1tilde > 0.0){/////// L1 term
	      for(int m=0;m<DIM;m++){
		for(int n=0;n<DIM;n++){
		  dmy_sigma[0][0] += SQ(dqij[0][m][n]);
		  dmy_sigma[1][1] += SQ(dqij[1][m][n]);
		  dmy_sigma[2][2] += SQ(dqij[2][m][n]);
		  
		  dmy_sigma[0][1] += dqij[0][m][n]*dqij[1][m][n];
		  dmy_sigma[0][2] += dqij[0][m][n]*dqij[2][m][n];
		  dmy_sigma[1][2] += dqij[1][m][n]*dqij[2][m][n];
		}
	      }
	      {
		for(int m=0;m<DIM;m++){
		  dmy_sigma[m][m] *= -L1tilde;
		  for(int n=m+1;n<DIM;n++){
		    dmy_sigma[m][n] *= -L1tilde;
		    dmy_sigma[n][m] = dmy_sigma[m][n];
		  }
		}
	      }
	      for(int m=0;m<DIM;m++){
		for(int n=0;n<DIM;n++){
		  sigma[m][n] += dmy_sigma[m][n];
		}
	      }
	    }
	    assert(L2tilde >= 0.0); 
	    if(L2tilde > 0.0){ /////// L2 term
	      for(int d=0;d<DIM;d++){
		v[d] = 0.0;
		for(int l=0;l<DIM;l++){
		  v[d] += dqij[l][l][d];
		}
	      }
	      for(int m=0;m<DIM;m++){
		for(int n=0;n<DIM;n++){
		  dmy_sigma[m][n] = 0.0;
		  for(int d=0;d<DIM;d++){
		    dmy_sigma[m][n] += v[d] * dqij[n][m][d];
		  }
		}
	      }
	      for(int m=0;m<DIM;m++){
		for(int n=0;n<DIM;n++){
		  sigma[m][n] += (-L2tilde*dmy_sigma[m][n]);
		}
	      }
	    }
	  }
	  {  // end of distortion and advections 
	    for(int m=0;m<DIM;m++){
	      for(int n=0;n<DIM;n++){
		stress[m*DIM+n][i][j][k] = sigma[m][n];
	      }
	    }
	    for(int d=0;d<QDIM;d++){
	      q_rhs[d][i][j][k] = dmy_q_rhs[d];
	    }
	  }
	}
      }
    }
  }
}
inline void QS_rhs2_QA(const Value Q_r[QDIM], const Value A_r[QDIM], Value stress[DIM*DIM], Value q_rhs[QDIM]){
  ///// QQ:A term and QA+AQ term and A term
  for(int i=0; i<N2X; i++){
    for(int j=0; j<N2Y; j++){
      for(int k=0; k<N2Z; k++){
	double q_full[DIM][DIM];
	Q2q_full(Q_r, i, j, k, q_full);
	double rate_of_strain[DIM][DIM];
	Q2q_full(A_r,i,j,k, rate_of_strain);
	
	double qa[DIM][DIM];
	AdotB(q_full, rate_of_strain, qa);
	double qddota=0.;
	for(int m=0;m<DIM;m++){
	  qddota += qa[m][m];
	}
	{
	  double dmy1= Beta_1 * qddota;
	  for(int m=0;m<DIM;m++){
	    stress[m*DIM+m][i][j][k] += (
					 dmy1 * q_full[m][m]
					 + Beta_56 * qa[m][m]
					 +Isotropic_viscosity 
					 * rate_of_strain[m][m]
					 );
	    for(int n=m+1;n<DIM;n++){
	      double dmy =(
			   dmy1 * q_full[m][n]
			   + Beta_56 * .5* (qa[m][n] + qa[n][m])
			   +Isotropic_viscosity * rate_of_strain[m][n]
			   );
	      stress[m*DIM+n][i][j][k] += dmy;
	      stress[n*DIM+m][i][j][k] += dmy;
	    }
	  }
	}
	{
	  q_rhs[0][i][j][k] -= Mu_2overMu_1_half * rate_of_strain[0][0];
	  q_rhs[1][i][j][k] -= Mu_2overMu_1_half * rate_of_strain[0][1];
	  q_rhs[2][i][j][k] -= Mu_2overMu_1_half * rate_of_strain[0][2];
	  q_rhs[3][i][j][k] -= Mu_2overMu_1_half * rate_of_strain[1][1];
	  q_rhs[4][i][j][k] -= Mu_2overMu_1_half * rate_of_strain[1][2];
	}
      }
    }
  }
}
inline void OG_rhs2_A(const Value A_r[QDIM], Value stress[DIM*DIM], Value q_rhs[QDIM]){
  ///// A term
  for(int i=0; i<N2X; i++){
    for(int j=0; j<N2Y; j++){
      for(int k=0; k<N2Z; k++){
	double rate_of_strain[DIM][DIM];
	Q2q_full(A_r,i,j,k, rate_of_strain);
	{
	  for(int m=0;m<DIM;m++){
	    stress[m*DIM+m][i][j][k] += 
	      Isotropic_viscosity * rate_of_strain[m][m];
	    for(int n=m+1;n<DIM;n++){
	      double dmy = Isotropic_viscosity * rate_of_strain[m][n];
	      stress[m*DIM+n][i][j][k] += dmy;
	      stress[n*DIM+m][i][j][k] += dmy;
	    }
	  }
	}
	{
	  q_rhs[0][i][j][k] -= Mu_2overMu_1_half * rate_of_strain[0][0];
	  q_rhs[1][i][j][k] -= Mu_2overMu_1_half * rate_of_strain[0][1];
	  q_rhs[2][i][j][k] -= Mu_2overMu_1_half * rate_of_strain[0][2];
	  q_rhs[3][i][j][k] -= Mu_2overMu_1_half * rate_of_strain[1][1];
	  q_rhs[4][i][j][k] -= Mu_2overMu_1_half * rate_of_strain[1][2];
	}
      }
    }
  }
}

inline void QS_rhs2_hQ(const Value Q_k[QDIM], const Value Q_r[QDIM], Value stress[DIM*DIM], Value q_rhs[QDIM]){
  ///// molcular field terms
  Q_k2h_linear_k(Q_k, tmp2_QDIM);
  for(int d=0;d<QDIM;d++){
    K_space_extension(tmp2_QDIM[d]);
    A_k2a_extended(tmp2_QDIM[d]);
  }
  
  double linear_max =-DBL_MAX;
  double linear_min =DBL_MAX;
  double qq_trace_max = 0.0;
  double qq_trace_min = DBL_MAX;
  for(int i=0; i<N2X; i++){
    for(int j=0; j<N2Y; j++){
      for(int k=0; k<N2Z; k++){
	double q_full[DIM][DIM];
	Q2q_full(Q_r, i, j, k, q_full);
	double qqij[DIM][DIM];
	AdotB(q_full, q_full, qqij);
	double qq_trace= 0.;
	for(int m=0;m<DIM;m++){
	  qq_trace += qqij[m][m];
	}
	double dmy_c = -A_LdG-C_LdG * qq_trace;
	{
	  double dmy = - dmy_c;
	  linear_min = MIN(dmy, linear_min);
	  linear_max = MAX(dmy, linear_max);
	  qq_trace_max = MAX(qq_trace_max, qq_trace); 
	  qq_trace_min = MIN(qq_trace_min, qq_trace); 
	}			 
	qq_trace *= (One_third * B_LdG);
	
	double h_full[DIM][DIM];
	Q2q_full(tmp2_QDIM, i, j, k, h_full);
	{
	  for(int m=0;m<DIM;m++){
	    for(int n=m;n<DIM;n++){
	      h_full[m][n] += (dmy_c * q_full[m][n]+B_LdG *qqij[m][n]);
	    }
	  }
	  for(int m=0;m<DIM;m++){
	    h_full[m][m] -= qq_trace;
	  }
	  if(1){
	    for(int m=0;m<DIM;m++){
	      for(int n=m+1;n<DIM;n++){
		h_full[n][m] = h_full[m][n];
	      }
	    }
	  }
	}
	if(1){
	  AdotB(h_full, q_full, qqij);
	  
	  double sigma[DIM][DIM]={{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
	  {///// anti-symmetric part of viscous stress
	    sigma[1][0] = sigma[0][1] = (qqij[0][1] - qqij[1][0]); 
	    sigma[2][0] = sigma[0][2] = (qqij[0][2] - qqij[2][0]); 
	    sigma[2][1] = sigma[1][2] = (qqij[1][2] - qqij[2][1]); 
	  }
	  for(int m=0;m<DIM;m++){
	    sigma[m][m] += Mu_2overMu_1_half * h_full[m][m];
	    for(int n=m+1;n<DIM;n++){
	      double dmy = Mu_2overMu_1_half * h_full[m][n];
	      sigma[m][n] += dmy;
	      sigma[n][m] += dmy;
	    }
	  }
	  for(int m=0;m<DIM;m++){
	    for(int n=0;n<DIM;n++){
	      stress[m*DIM+n][i][j][k] += sigma[m][n];
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
  if(1){
    fprintf(stdout, "%g %g %g %g %g %g\n", -linear_max, -linear_min
	    ,qq_trace_max
	    ,qq_trace_min
	    , sqrt(qq_trace_max*Q2S_prefactor)/A_eq
	    , sqrt(qq_trace_min*Q2S_prefactor)/A_eq
	    );
    //fflush(stdout);
  }
}

void Anchoring_energy_harmonic_wall_gradient(Value q_rhs_extended[QDIM]
					     ,const Value Qr_extended[QDIM]
					     ){
  double dmy = 0.;
  for(int n=0;n<Component_Number;n++){
    dmy = MAX(W_surfaces[n],dmy);
  }
  const double anchor_coeff = dmy;
  
  double q_full_anchor[DIM][DIM];
  const double anchor_director[DIM] = {0.,0.,1.};
  const double anchor_scalar_order = 1.;
  
  {// pretilt tensor order 
    for(int m=0;m<DIM;m++){
      for(int n=m;n<DIM;n++){
	q_full_anchor[m][n] = anchor_director[m] * anchor_director[n];
      }
    }
    for(int m=0;m<DIM;m++){
      q_full_anchor[m][m] -= One_third;
    }
    for(int m=0;m<DIM;m++){
      for(int n=m;n<DIM;n++){
	q_full_anchor[m][n] *= anchor_scalar_order;
      }
    }
  }
  
  for(int i=0; i<N2X; i++){
    for(int j=0; j<N2Y; j++){
      for(int k=0; k<N2Z; k++){
	if( i == 0 || j == 0 ||k== 0){
	  double q_full[DIM][DIM];
	  Q2q_full(Qr_extended, i,j,k, q_full);
	  
	  double h_anchor[DIM][DIM];
	  for(int m=0;m<DIM;m++){
	    for(int n=m;n<DIM;n++){
	      h_anchor[m][n] = 
		-anchor_coeff * (
				 q_full[m][n] - q_full_anchor[m][n]
				 );
	    }
	  }
	  
	  q_rhs_extended[0][i][j][k] += h_anchor[0][0];
	  q_rhs_extended[1][i][j][k] += h_anchor[0][1];
	  q_rhs_extended[2][i][j][k] += h_anchor[0][2];
	  q_rhs_extended[3][i][j][k] += h_anchor[1][1];
	  q_rhs_extended[4][i][j][k] += h_anchor[1][2];
	}
      }
    }
  }
}

void Anchoring_energy_harmonic_gradient(Value q_rhs_extended[QDIM]
			       ,const Value Qr_extended[QDIM]
			       ,const Value DPHI_extended[DIM]
			       ,const Particle *p
			       ){
  const double dx = 0.5 * DX;
  int **sekibun_cell = Sekibun_cell_extended;
  int *Nlattice = N2s;
  int np_domain = NP_domain_extended;


  for(int n=0; n < Particle_Number; n++){
    double xp[DIM];
    for(int d=0;d<DIM;d++){
      xp[d] = p[n].x[d];
    }
    const double anchor_coeff = W_surfaces[p[n].spec];
    const double anchor_scalar_order = S_surfaces[p[n].spec];

    
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell 
      = Particle_cell(xp, dx, x_int, residue);// {1,0} が返ってくる
    sw_in_cell = 1;
    int r_mesh[DIM];
    double r[DIM];
    for(int mesh=0; mesh < np_domain; mesh++){
      Relative_coord(sekibun_cell[mesh]
		     , x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
      double dphi[DIM];
      double amp2_dphi = 0.;
      for(int m=0;m<DIM;m++){
	dphi[m] = DPHI_extended[m][r_mesh[0]][r_mesh[1]][r_mesh[2]];
	amp2_dphi += SQ(dphi[m]);
      }
      if(amp2_dphi > 0.0){
	double amp1_dphi = sqrt(amp2_dphi);
	double iamp2_dphi = 1./amp2_dphi;

	double q_full[DIM][DIM];
	Q2q_full(Qr_extended, r_mesh[0], r_mesh[1], r_mesh[2], q_full);

	double q_full_anchor[DIM][DIM];

	{// pretilt tensor order 
	  for(int m=0;m<DIM;m++){
	    for(int n=m;n<DIM;n++){
	      q_full_anchor[m][n] = dphi[m] * dphi[n] * iamp2_dphi;
	    }
	  }
	  for(int m=0;m<DIM;m++){
	    q_full_anchor[m][m] -= One_third;
	  }
	  for(int m=0;m<DIM;m++){
	    for(int n=m;n<DIM;n++){
	      q_full_anchor[m][n] *= anchor_scalar_order;
	    }
	  }
	}

	double h_anchor[DIM][DIM];
	const double dmy = anchor_coeff * amp1_dphi;
	for(int m=0;m<DIM;m++){
	  for(int n=m;n<DIM;n++){
	    h_anchor[m][n] = 
	      -dmy * (
		      q_full[m][n] - q_full_anchor[m][n]
		      );
	  }
	}
	q_rhs_extended[0][r_mesh[0]][r_mesh[1]][r_mesh[2]] += h_anchor[0][0];
	q_rhs_extended[1][r_mesh[0]][r_mesh[1]][r_mesh[2]] += h_anchor[0][1];
	q_rhs_extended[2][r_mesh[0]][r_mesh[1]][r_mesh[2]] += h_anchor[0][2];
	q_rhs_extended[3][r_mesh[0]][r_mesh[1]][r_mesh[2]] += h_anchor[1][1];
	q_rhs_extended[4][r_mesh[0]][r_mesh[1]][r_mesh[2]] += h_anchor[1][2];
      }
    }
  }
}
void Anchoring_energy_director_gradient(Value q_rhs_extended[QDIM]
					,const Value Qr_extended[QDIM]
					,const Value DPHI_extended[DIM]
					,const Particle *p
					){
  const double dx = 0.5 * DX;
  int **sekibun_cell = Sekibun_cell_extended;
  int *Nlattice = N2s;
  int np_domain = NP_domain_extended;


  for(int n=0; n < Particle_Number; n++){
    double xp[DIM];
    for(int d=0;d<DIM;d++){
      xp[d] = p[n].x[d];
    }
    const double anchor_coeff = W_surfaces[p[n].spec];
    
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell 
      = Particle_cell(xp, dx, x_int, residue);// {1,0} が返ってくる
    sw_in_cell = 1;
    int r_mesh[DIM];
    double r[DIM];
    for(int mesh=0; mesh < np_domain; mesh++){
      Relative_coord(sekibun_cell[mesh]
		     , x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
      double dphi[DIM];
      double amp2_dphi = 0.;
      for(int m=0;m<DIM;m++){
	dphi[m] = DPHI_extended[m][r_mesh[0]][r_mesh[1]][r_mesh[2]];
	amp2_dphi += SQ(dphi[m]);
      }
      if(amp2_dphi > 0.0){
	double amp1_dphi = sqrt(amp2_dphi);
	double iamp2_dphi = 1./amp2_dphi;

	double q_full[DIM][DIM];
	Q2q_full(Qr_extended, r_mesh[0], r_mesh[1], r_mesh[2], q_full);
	double qqtrace = 0.0;
	for(int m=0;m<DIM;m++){
	  qqtrace += SQ(q_full[m][m]);
	  for(int n=m+1;n<DIM;n++){
	    qqtrace += 2.*SQ(q_full[m][n]);
	  }
	}
	if(qqtrace > 0.0){
	  double iscalar2 = 1./qqtrace;
	  double iscalar1 = sqrt(iscalar2);
	  
	  double dmy1 = 0.0;
	  for(int m=0;m<DIM;m++){
	    dmy1 += SQ(dphi[m])* q_full[m][m];
	    for(int n=m+1;n<DIM;n++){
	      dmy1 += (2.* dphi[m] * dphi[n] * q_full[m][n]);
	    }
	  }
	  dmy1 *= (0.5 * iscalar2 * iamp2_dphi);
	  
	  double h_anchor[DIM][DIM];
	  const double dmy = anchor_coeff * amp1_dphi * iscalar1;
	  for(int m=0;m<DIM;m++){
	    for(int n=m;n<DIM;n++){
	      h_anchor[m][n] = -dmy * (dphi[m] * dphi[n] * iamp2_dphi
				      - dmy1 * q_full[m][n]);
	    }
	  }
	  q_rhs_extended[0][r_mesh[0]][r_mesh[1]][r_mesh[2]] += h_anchor[0][0];
	  q_rhs_extended[1][r_mesh[0]][r_mesh[1]][r_mesh[2]] += h_anchor[0][1];
	  q_rhs_extended[2][r_mesh[0]][r_mesh[1]][r_mesh[2]] += h_anchor[0][2];
	  q_rhs_extended[3][r_mesh[0]][r_mesh[1]][r_mesh[2]] += h_anchor[1][1];
	  q_rhs_extended[4][r_mesh[0]][r_mesh[1]][r_mesh[2]] += h_anchor[1][2];
	}
      }
    }
  }
}
void Anchoring_energy_sdirector_gradient(Value q_rhs_extended[QDIM]
					 ,const Value Qr_extended[QDIM]
					 ,const Value DPHI_extended[DIM]
					 ,const Particle *p
					){
  const double dx = 0.5 * DX;
  int **sekibun_cell = Sekibun_cell_extended;
  int *Nlattice = N2s;
  int np_domain = NP_domain_extended;
  
  
  for(int n=0; n < Particle_Number; n++){
    double xp[DIM];
    for(int d=0;d<DIM;d++){
      xp[d] = p[n].x[d];
    }
    const double anchor_coeff = W_surfaces[p[n].spec];
    
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell 
      = Particle_cell(xp, dx, x_int, residue);// {1,0} が返ってくる
    sw_in_cell = 1;
    int r_mesh[DIM];
    double r[DIM];
    for(int mesh=0; mesh < np_domain; mesh++){
      Relative_coord(sekibun_cell[mesh]
		     , x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);
      double dphi[DIM];
      double amp2_dphi = 0.;
      for(int m=0;m<DIM;m++){
	dphi[m] = DPHI_extended[m][r_mesh[0]][r_mesh[1]][r_mesh[2]];
	amp2_dphi += SQ(dphi[m]);
      }
      if(amp2_dphi > 0.0){
	double amp1_dphi = sqrt(amp2_dphi);
	double iamp2_dphi = 1./amp2_dphi;
	
	{
	  double h_anchor[DIM][DIM];
	  const double dmy = anchor_coeff * amp1_dphi;
	  for(int m=0;m<DIM;m++){
	    for(int n=m;n<DIM;n++){
	      h_anchor[m][n] = -dmy * dphi[m] * dphi[n] * iamp2_dphi;
	    }
	  }
	  q_rhs_extended[0][r_mesh[0]][r_mesh[1]][r_mesh[2]] += h_anchor[0][0];
	  q_rhs_extended[1][r_mesh[0]][r_mesh[1]][r_mesh[2]] += h_anchor[0][1];
	  q_rhs_extended[2][r_mesh[0]][r_mesh[1]][r_mesh[2]] += h_anchor[0][2];
	  q_rhs_extended[3][r_mesh[0]][r_mesh[1]][r_mesh[2]] += h_anchor[1][1];
	  q_rhs_extended[4][r_mesh[0]][r_mesh[1]][r_mesh[2]] += h_anchor[1][2];
	}
      }
    }
  }
}

void AC_rhs(Value Q_k[QDIM]
	    ,Value q_rhs[QDIM]
	    ,Particle *p){

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
    for(int d=0;d<QDIM;d++){
      K_space_extension(tmp1_QDIM[d]);
      A_k2a_extended(tmp1_QDIM[d]);
    }
  }
  ///// molcular field terms
  Q_k2h_linear_k(Q_k, q_rhs);
  for(int d=0;d<QDIM;d++){
    K_space_extension(q_rhs[d]);
    A_k2a_extended(q_rhs[d]);
  }

  if(Particle_Number> 0){
    {
      // nabla phi(r)
      Reset_phi_extended(tmp2_QDIM[0]);
      Make_phi_particle_extended(tmp2_QDIM[0],p);
      A2a_k_extended(tmp2_QDIM[0]);
      A_k2dxa_k_extended(tmp2_QDIM[0], tmp2_QDIM[1]);
      A_k2dya_k_extended(tmp2_QDIM[0], tmp2_QDIM[2]);
      A_k2dza_k_extended(tmp2_QDIM[0], tmp2_QDIM[3]);
      for(int d=0;d<DIM;d++){
	A_k2a_extended(tmp2_QDIM[d+1]);
      }
    }
    Anchoring_energy_harmonic_gradient(q_rhs, tmp1_QDIM, tmp2_QDIM+1, p);
    Anchoring_energy_harmonic_wall_gradient(q_rhs, tmp1_QDIM);

    //Anchoring_energy_director_gradient(q_rhs, tmp1_QDIM, tmp2_QDIM+1, p);
    //Anchoring_energy_sdirector_gradient(q_rhs, tmp1_QDIM, tmp2_QDIM+1, p);
    { // phi(r)
      Reset_phi_extended(tmp2_QDIM[0]);
      Make_phi_particle_extended(tmp2_QDIM[0],p);
    }
  }else {
    Reset_phi_extended(tmp2_QDIM[0]);
  }

  double linear_max =-DBL_MAX;
  double linear_min =DBL_MAX;
  double qq_trace_max = 0.0;
  double qq_trace_min = DBL_MAX;
  double qq_trace_mean = 0.0;
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
	double dmy_c = -A_LdG-C_LdG * qq_trace;
	{
	  double dmy = - dmy_c;
	  {
	    double dmy_phif = (1.-tmp2_QDIM[0][i][j][k]);
	    dmy_c *= dmy_phif;
	    qq_trace *= dmy_phif;
	    if(dmy_phif > 0.){
	      linear_min = MIN(dmy, linear_min);
	      linear_max = MAX(dmy, linear_max);
	      qq_trace_max = MAX(qq_trace_max, qq_trace); 
	      qq_trace_min = MIN(qq_trace_min, qq_trace); 
	      qq_trace_mean += qq_trace;
	    }
	  }	
	}		 
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
	  if(0){
	    for(int m=0;m<DIM;m++){
	      for(int n=m+1;n<DIM;n++){
		h_full[n][m] = h_full[m][n];
	      }
	    }
	  }
	}
	if(1){
	  double dmy = One_overMu1 * (1.-tmp2_QDIM[0][i][j][k]);
	  q_rhs[0][i][j][k] += dmy * h_full[0][0];
	  q_rhs[1][i][j][k] += dmy * h_full[0][1];
	  q_rhs[2][i][j][k] += dmy * h_full[0][2];
	  q_rhs[3][i][j][k] += dmy * h_full[1][1];
	  q_rhs[4][i][j][k] += dmy * h_full[1][2];
	}
      }
    }
  }
  if(1){
    fprintf(stdout, "%g %g %g %g %g %g %g\n"
	    , -linear_max, -linear_min
	    ,qq_trace_max
	    ,qq_trace_min
	    , sqrt(qq_trace_max*Q2S_prefactor)/A_eq
	    , sqrt(qq_trace_min*Q2S_prefactor)/A_eq
	    ,qq_trace_mean /(double)(N2X*N2Y*N2Z)
	    );
    //fflush(stdout);
  }
  {
    for(int d=0;d<QDIM;d++){
      A2a_k_extended(q_rhs[d]);
      K_space_extension_reverse(q_rhs[d]);
    }    
  }
}
inline void omega2Omega_full(const Value *omega, const int i,const int j, const int k, double Omega_full[][DIM]){
    Omega_full[0][0] = 0.0;
    Omega_full[1][1] = 0.0;
    Omega_full[2][2] = 0.0;
    Omega_full[1][2] = .5*omega[0][i][j][k]; 
    Omega_full[2][0] = .5*omega[1][i][j][k]; 
    Omega_full[0][1] = .5*omega[2][i][j][k]; 

    Omega_full[2][1] = -Omega_full[1][2];
    Omega_full[0][2] = -Omega_full[2][0];
    Omega_full[1][0] = -Omega_full[0][1];
}
inline void QS_rhs2_QOmega(const Value omega[DIM], const Value Q_r[QDIM], Value q_rhs[QDIM]){
  
  for(int i=0; i<N2X; i++){
    for(int j=0; j<N2Y; j++){
      for(int k=0; k<N2Z; k++){
	double q_full[DIM][DIM];
	Q2q_full(Q_r, i, j, k, q_full);
	double Omega_full[DIM][DIM];
	omega2Omega_full(omega,i,j,k,Omega_full);

	double q_omega[DIM][DIM];
	AdotB(q_full, Omega_full, q_omega);
	for(int m=0;m<DIM;m++){
	  q_omega[m][m] *= 2.0;
	  for(int n=m+1;n<DIM;n++){
	    double dmy = q_omega[m][n]+q_omega[n][m];
	    q_omega[m][n] = dmy;
	    q_omega[n][m] = dmy;
	  }
	}
	double trace_q_omega = 0.0;
	for(int m=0;m<DIM;m++){
	  trace_q_omega += q_omega[m][m];
	}
	trace_q_omega *= One_third;
	for(int m=0;m<DIM;m++){
	  q_omega[m][m] -= trace_q_omega;
	}
	{
	  q_rhs[0][i][j][k] += q_omega[0][0];
	  q_rhs[1][i][j][k] += q_omega[0][1];
	  q_rhs[2][i][j][k] += q_omega[0][2];
	  q_rhs[3][i][j][k] += q_omega[1][1];
	  q_rhs[4][i][j][k] += q_omega[1][2];
	}
      }
    }
  }
}
inline void OG_rhs2(const Value zeta[DIM-1], const Value q_k_orig[QDIM], const Value tensor_order_r[QDIM], Value stress[DIM*DIM], Value q_rhs[QDIM]){

  Zeta_k2Strain_k(zeta, tmp2_QDIM);
  for(int d=0;d<QDIM;d++){
    K_space_extension(tmp2_QDIM[d]);
    A_k2a_extended(tmp2_QDIM[d]);
  }
  
  OG_rhs2_A(tmp2_QDIM, stress, q_rhs);
  QS_rhs2_hQ(q_k_orig, tensor_order_r, stress, q_rhs);

  Zeta_k2omega_k(zeta, tmp2_QDIM);
  for(int d=0;d<DIM;d++){
    K_space_extension(tmp2_QDIM[d]);
    A_k2a_extended(tmp2_QDIM[d]);
  }
  QS_rhs2_QOmega(tmp2_QDIM, tensor_order_r, q_rhs);

}
inline void QS_rhs2(const Value zeta[DIM-1], const Value q_k_orig[QDIM], const Value tensor_order_r[QDIM], Value stress[DIM*DIM], Value q_rhs[QDIM]){

  Zeta_k2Strain_k(zeta, tmp2_QDIM);
  for(int d=0;d<QDIM;d++){
    K_space_extension(tmp2_QDIM[d]);
    A_k2a_extended(tmp2_QDIM[d]);
  }

  QS_rhs2_QA(tensor_order_r, tmp2_QDIM, stress, q_rhs);
  QS_rhs2_hQ(q_k_orig, tensor_order_r, stress, q_rhs);

  Zeta_k2omega_k(zeta, tmp2_QDIM);
  for(int d=0;d<DIM;d++){
    K_space_extension(tmp2_QDIM[d]);
    A_k2a_extended(tmp2_QDIM[d]);
  }
  QS_rhs2_QOmega(tmp2_QDIM, tensor_order_r, q_rhs);

}
inline void Truncate_Qian_Sheng(Value zeta[DIM-1], Value q[QDIM]){
  for(int d=0;d<DIM-1;d++){
    Truncate_four_fifths_rule(zeta[d]);
  }
  for(int d=0;d<QDIM;d++){
    Truncate_four_fifths_rule(q[d]);
  }
}
void QS_rhs(Value zeta[DIM-1], const double uk_dc[DIM], Value q_k_orig[QDIM], Value stress[DIM*DIM], Value q_rhs[QDIM], const Index_range *ijk_range, const int &n_ijk_range){
  Truncate_Qian_Sheng(zeta, q_k_orig);
  QS_rhs1(zeta, uk_dc, q_k_orig, stress, q_rhs);
  {
    for(int d=0;d<QDIM;d++){
      for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	  for(int k=0; k<NZ_; k++){
	    tmp1_QDIM[d][i][j][k]=q_k_orig[d][i][j][k];
	  }
	}
      }
    }
    for(int d=0;d<QDIM;d++){
      K_space_extension(tmp1_QDIM[d]);
      A_k2a_extended(tmp1_QDIM[d]);
    }
    QS_rhs2(zeta, q_k_orig, tmp1_QDIM, stress, q_rhs);
  }
  {
    for(int d=0;d<QDIM;d++){
      A2a_k_extended(q_rhs[d]);
      K_space_extension_reverse(q_rhs[d]);
    }    
    for(int d=0;d<DIM*DIM;d++){
      A2a_k_extended(stress[d]);
      K_space_extension_reverse(stress[d]);
    }
  }
  Stress_k2zeta_rhs_k(stress, ijk_range, n_ijk_range);
}

void OG_rhs(Value zeta[DIM-1], const double uk_dc[DIM], Value q_k_orig[QDIM], Value stress[DIM*DIM], Value q_rhs[QDIM], const Index_range *ijk_range, const int &n_ijk_range){
  Truncate_Qian_Sheng(zeta, q_k_orig);
  QS_rhs1(zeta, uk_dc, q_k_orig, stress, q_rhs);
  {
    for(int d=0;d<QDIM;d++){
      for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	  for(int k=0; k<NZ_; k++){
	    tmp1_QDIM[d][i][j][k]=q_k_orig[d][i][j][k];
	  }
	}
      }
    }
    for(int d=0;d<QDIM;d++){
      K_space_extension(tmp1_QDIM[d]);
      A_k2a_extended(tmp1_QDIM[d]);
    }
    OG_rhs2(zeta, q_k_orig, tmp1_QDIM, stress, q_rhs);
  }
  {
    for(int d=0;d<QDIM;d++){
      A2a_k_extended(q_rhs[d]);
      K_space_extension_reverse(q_rhs[d]);
    }    
    for(int d=0;d<DIM*DIM;d++){
      A2a_k_extended(stress[d]);
      K_space_extension_reverse(stress[d]);
    }
  }
  Stress_k2zeta_rhs_k(stress, ijk_range, n_ijk_range);
}



////////////////// FP stress
void Advection_slippy_NS(Value zeta[DIM-1], const double uk_dc[DIM], Value rhs[DIM-1], const Index_range *ijk_range, const int &n_ijk_range, Particle *p){
  {
    Truncate_vector_two_third_rule(zeta, DIM-1);
    Zeta_k2u(zeta, uk_dc, u);

    if(STOKES){
      Reset_phi_u(phi, f_particle);
      Make_phi_u_advection(phi, f_particle, p);
      Solenoidal_u(f_particle);
    }else {
      Reset_phi_u(phi, up);
      Make_phi_u_advection(phi, up, p);
      Make_f_particle_dt_sole(f_particle, u, up, phi);
      Add_f_particle(f_particle, u);
    }
  }

  {
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ; k++){
	  double v[DIM];
	  double v_advection[DIM];
	  for(int d=0;d<DIM;d++){
	    v[d] = u[d][i][j][k];
	    v_advection[d] = f_particle[d][i][j][k];
	  }
	  static double sigma[DIM][DIM];
	  {///// advection of v
	    for(int m=0;m<DIM;m++){
	      for(int n=0;n<DIM;n++){
		sigma[m][n] = -RHO* v_advection[m] * v[n];
	      }
	    }
	  }
	  {
	    for(int m=0;m<DIM-1;m++){
	      sigma[m][m] -= sigma[DIM-1][DIM-1];
	    }
	    sigma[DIM-1][DIM-1] = 0.0;
	    for(int m=0;m<DIM;m++){
	      for(int n=0;n<DIM;n++){
		Stress[m*DIM+n][i][j][k] = sigma[m][n];
	      }
	    }
	  }
	}
      }
    }
    {
      Reset_phi(Stress[DIM*DIM-1]);
      for(int d=0;d<DIM*DIM-1;d++){
	A2a_k(Stress[d]);
      }
      Stress_k2zeta_rhs_k(Stress, ijk_range, n_ijk_range);
    }
  }
  { 
    for(int n=0;n<n_ijk_range;n++){
      for(int d=0;d<DIM-1;d++){
	for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	  for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	    for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	      rhs[d][i][j][k] += Stress[d][i][j][k];
	    }
	  }
	}
      }
    }
  }
}
void Rhs_slippy_NS_stress(Value zeta[DIM-1], const double uk_dc[DIM], Value rhs[DIM-1], const Index_range *ijk_range, const int &n_ijk_range, Particle *p){
  {// making phi
    Reset_phi(phi);
    Make_phi_u_particle(phi, up, p);
  }
  {
    Zeta_k2Strain_k(zeta, tmp2_QDIM);
    for(int d=0;d<QDIM;d++){
      A_k2a(tmp2_QDIM[d]);
    }
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ; k++){
	  static double sigma[DIM][DIM];
	  Q2q_full(tmp2_QDIM,i,j,k, sigma);
	  double fpd_viscosity = 2.* (ETA + phi[i][j][k]* Delta_ETA) ;
	  { // FP viscos stress
	    for(int m=0;m<DIM;m++){
	      for(int n=0;n<DIM;n++){
		sigma[m][n] *= fpd_viscosity;
	      }
	    }
	  }
	  {
	    for(int m=0;m<DIM;m++){
	      sigma[m][m] -= sigma[DIM-1][DIM-1];
	    }
	    sigma[DIM-1][DIM-1] = 0.0;
	    for(int m=0;m<DIM;m++){
	      for(int n=0;n<DIM;n++){
		Stress[m*DIM+n][i][j][k] = sigma[m][n];
	      }
	    }
	  }
	}
      }
    }
    {
      Reset_phi(Stress[DIM*DIM-1]);
      for(int d=0;d<DIM*DIM-1;d++){
	A2a_k(Stress[d]);
      }
      Stress_k2zeta_rhs_k(Stress, ijk_range, n_ijk_range);
    }
  }
  { 
    for(int n=0;n<n_ijk_range;n++){
      for(int d=0;d<DIM-1;d++){
	for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	  for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	    for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	      rhs[d][i][j][k] += Stress[d][i][j][k];
	    }
	  }
	}
      }
    }
  }
}
void Rhs_slippy_NS_surface_normal(Value zeta[DIM-1], double force_dc[DIM], Value rhs[DIM-1], const Index_range *ijk_range, const int &n_ijk_range, Particle *p){
  {// making phi
    Reset_phi(phi);
    Make_phi_u_particle(phi, up, p);

    A2a_k(phi);
    A_k2da_k(phi, up);
    U_k2u(up);
  }
  
  {
    Zeta_k2Strain_k(zeta, tmp2_QDIM);
    for(int d=0;d<QDIM;d++){
      A_k2a(tmp2_QDIM[d]);
    }
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
	for(int k=0; k<NZ; k++){
	  static double sigma[DIM][DIM];
	  Q2q_full(tmp2_QDIM,i,j,k, sigma);
	  {// FP surface normal stress
	    double normal[DIM]; 
	    double norm2 =0.0;
	    for(int m=0;m<DIM;m++){
	      normal[m] = up[m][i][j][k];
	      norm2 += SQ(normal[m]);
	    }
	    double dmy = 0.0;
	    for(int m=0;m<DIM;m++){
	      for(int n=0;n<DIM;n++){
		dmy += sigma[m][n] * normal[m] * normal[n];
	      }
	    }
	    if(norm2 > 0.0){
	      dmy /= norm2;
	    }
	    dmy *= (Delta_ETA/RHO * 2.0)*1.e2;
	    for(int m=0;m<DIM;m++){
	      //up[m][i][j][k] *= dmy;
	      up[m][i][j][k] *= Delta_ETA;//*1.e2;
	    }
	  }
	}
      }
    }
  }
  { 
    double dmy[DIM];
    U2zeta_k(u, dmy, up);
    for(int d=0;d<DIM;d++){
      force_dc[d] += dmy[d];
    }
    for(int n=0;n<n_ijk_range;n++){
      for(int d=0;d<DIM-1;d++){
	for(int i=ijk_range[n].istart; i<=ijk_range[n].iend; i++){
	  for(int j=ijk_range[n].jstart; j<=ijk_range[n].jend; j++){
	    for(int k=ijk_range[n].kstart; k<=ijk_range[n].kend; k++){
	      rhs[d][i][j][k] += u[d][i][j][k];
	    }
	  }
	}
      }
    }
  }
}
