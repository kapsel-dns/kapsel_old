//
// $Id: fft_wrapper.cxx,v 1.21 2006/05/19 04:59:49 kin Exp $
//

#include "fft_wrapper.h"

int *ip;
double *w;
double *t;

int *ip_extended;
double *w_extended;
double *t_extended;

Value u[DIM], up[DIM], phi;
Value_int KX_int, KY_int, KZ_int;
Value_int KX_int_extended, KY_int_extended, KZ_int_extended;
Value K2,IK2;

int (*Calc_KX)( const int &i, const int &j, const int &k);
int (*Calc_KY)( const int &i, const int &j, const int &k);
int (*Calc_KZ)( const int &i, const int &j, const int &k);
void (*Truncate_two_third_rule)(Value &a);
void (*Truncate_four_fifths_rule)(Value &a);



//// function definition
inline int Calc_KX_Ooura_primitive(const int &i
			 ,const int &j
			 ,const int &k
			 ,const int &nx,const int &hnx
			 ,const int &ny,const int &hny
			 ,const int &nz_,const int &hnz_
			 ){
  // i,j,k: 配列の添字, 
  return (i>hnx) ? i-nx:i;
}
inline int Calc_KY_Ooura_primitive(const int &i
			 ,const int &j
			 ,const int &k
			 ,const int &nx,const int &hnx
			 ,const int &ny,const int &hny
			 ,const int &nz_,const int &hnz_
			 ){
  // i,j,k: 配列の添字, 
  assert(i < nx);
  assert(j < ny);
  
  return (j>hny) ? j-ny:j;
}
inline int Calc_KZ_Ooura_primitive(const int &i
			 ,const int &j
			 ,const int &k
			 ,const int &nx,const int &hnx
			 ,const int &ny,const int &hny
			 ,const int &nz_,const int &hnz_
			 ){
  // i,j,k: 配列の添字, 
  assert(i < nx);
  assert(j < ny);
  
  return k/2;
}
inline int Calc_KX_Ooura(const int &i
			 ,const int &j
			 ,const int &k
			 ){
  return Calc_KX_Ooura_primitive(i,j,k,NX,HNX,NY,HNY,NZ_,HNZ_);
}
inline int Calc_KY_Ooura(const int &i
			 ,const int &j
			 ,const int &k
			 ){
  return Calc_KY_Ooura_primitive(i,j,k,NX,HNX,NY,HNY,NZ_,HNZ_);
}
inline int Calc_KZ_Ooura(const int &i
			 ,const int &j
			 ,const int &k
			 ){
  return Calc_KZ_Ooura_primitive(i,j,k,NX,HNX,NY,HNY,NZ_,HNZ_);
}

inline int Calc_KX_Ooura_extended(const int &i
				  ,const int &j
				  ,const int &k
				  ){
  return Calc_KX_Ooura_primitive(i,j,k,N2X,HN2X,N2Y,HN2Y,N2Z_,HN2Z_);
}
inline int Calc_KY_Ooura_extended(const int &i
				  ,const int &j
				  ,const int &k
				  ){
  return Calc_KY_Ooura_primitive(i,j,k,N2X,HN2X,N2Y,HN2Y,N2Z_,HN2Z_);
}
inline int Calc_KZ_Ooura_extended(const int &i
				  ,const int &j
				  ,const int &k
				  ){
  return Calc_KZ_Ooura_primitive(i,j,k,N2X,HN2X,N2Y,HN2Y,N2Z_,HN2Z_);
}
inline void Init_K(void){
  KX_int = alloc_3d_int(NX, NY, NZ_);
  KY_int = alloc_3d_int(NX, NY, NZ_);
  KZ_int = alloc_3d_int(NX, NY, NZ_);
  K2 = alloc_3d_double(NX, NY, NZ_);
  IK2 = alloc_3d_double(NX, NY, NZ_);
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ_; k++){
	int kx = Calc_KX(i,j,k);
	int ky = Calc_KY(i,j,k);
	int kz = Calc_KZ(i,j,k);
	
	KX_int[i][j][k] = kx;
	KY_int[i][j][k] = ky;
	KZ_int[i][j][k] = kz;
	
	K2[i][j][k] = SQ(WAVE_X*(double)kx) 
	  + SQ(WAVE_Y*(double)ky)
	  + SQ(WAVE_Z*(double)kz);
	if(K2[i][j][k] > 0.0){
	  IK2[i][j][k] = 1./K2[i][j][k];
	}else{
	  //IK2[i][j][k] = 1.0;
	  IK2[i][j][k] = 0.0;
	}
      }
    }
  }
}
inline void Init_K_extended(void){
  KX_int_extended = alloc_3d_int(N2X, N2Y, N2Z_);
  KY_int_extended = alloc_3d_int(N2X, N2Y, N2Z_);
  KZ_int_extended = alloc_3d_int(N2X, N2Y, N2Z_);
  for(int i=0; i<N2X; i++){
    for(int j=0; j<N2Y; j++){
      for(int k=0; k<N2Z_; k++){
	int kx = Calc_KX_Ooura_extended(i,j,k);
	int ky = Calc_KY_Ooura_extended(i,j,k);
	int kz = Calc_KZ_Ooura_extended(i,j,k);
	
	KX_int_extended[i][j][k] = kx;
	KY_int_extended[i][j][k] = ky;
	KZ_int_extended[i][j][k] = kz;
      }
    }
  }
}
inline void Truncate_general(Value &a, const Index_range &ijk_range){
  for(int i=ijk_range.istart; i<=ijk_range.iend; i++){
    for(int j=ijk_range.jstart; j<=ijk_range.jend; j++){
      for(int k=ijk_range.kstart; k<=ijk_range.kend; k++){
	//assert( abs(Calc_KY_Ooura(i,j,k))>= TRN_Y); 
	//assert( abs(Calc_KX_Ooura(i,j,k))>= TRN_X); 
	//assert( Calc_KZ_Ooura(i,j,k)>= TRN_Z); 
	assert( (abs(Calc_KY_Ooura(i,j,k))>= TRN_Y || abs(Calc_KX_Ooura(i,j,k))>= TRN_X || Calc_KZ_Ooura(i,j,k)>= TRN_Z)); 
	a[i][j][k] = 0.0;
      }
    }
  }
}
inline void Truncate_two_third_rule_ooura(Value &a){
  static Index_range dmy_range;
  const int trn_z2=2*TRN_Z;
  {
    dmy_range.istart=0;
    dmy_range.iend=NX-1;
    dmy_range.jstart=0;
    dmy_range.jend=NY-1;
    dmy_range.kstart=trn_z2;
    dmy_range.kend=NZ_-1;
    Truncate_general(a, dmy_range);
  }
  {
    dmy_range.istart=0;
    dmy_range.iend=NX-1;
    dmy_range.jstart=TRN_Y;
    dmy_range.jend=NY-TRN_Y;
    dmy_range.kstart=0;
    dmy_range.kend=trn_z2-1;
    Truncate_general(a, dmy_range);
  }
  {
    dmy_range.istart=TRN_X;
    dmy_range.iend=NX-TRN_X;
    dmy_range.jstart=0;
    dmy_range.jend=TRN_Y-1;
    dmy_range.kstart=0;
    dmy_range.kend=trn_z2-1;
    Truncate_general(a, dmy_range);
  }
  {
    dmy_range.istart=TRN_X;
    dmy_range.iend=NX-TRN_X;
    dmy_range.jstart=NY-TRN_Y+1;
    dmy_range.jend=NY-1;
    dmy_range.kstart=0;
    dmy_range.kend=trn_z2-1;
    Truncate_general(a, dmy_range);
  }
}
inline void Truncate_four_fifths_rule_ooura(Value &a){
  /////////////////////////////////////////
  static Index_range dmy_range;
  const int trn_z2=2*TRN_QS_Z;
  {
    dmy_range.istart=0;
    dmy_range.iend=NX-1;
    dmy_range.jstart=0;
    dmy_range.jend=NY-1;
    dmy_range.kstart=trn_z2;
    dmy_range.kend=NZ_-1;
    Truncate_general(a, dmy_range);
  }
  {
    dmy_range.istart=0;
    dmy_range.iend=NX-1;
    dmy_range.jstart=TRN_QS_Y;
    dmy_range.jend=NY-TRN_QS_Y;
    dmy_range.kstart=0;
    dmy_range.kend=trn_z2-1;
    Truncate_general(a, dmy_range);
  }
  {
    dmy_range.istart=TRN_QS_X;
    dmy_range.iend=NX-TRN_QS_X;
    dmy_range.jstart=0;
    dmy_range.jend=TRN_QS_Y-1;
    dmy_range.kstart=0;
    dmy_range.kend=trn_z2-1;
    Truncate_general(a, dmy_range);
  }
  {
    dmy_range.istart=TRN_QS_X;
    dmy_range.iend=NX-TRN_QS_X;
    dmy_range.jstart=NY-TRN_QS_Y+1;
    dmy_range.jend=NY-1;
    dmy_range.kstart=0;
    dmy_range.kend=trn_z2-1;
    Truncate_general(a, dmy_range);
  }
}
inline void Init_fft_ooura(void){
  fprintf(stderr,"# Ooura rdft3d is selected.\n");
  
  {
    int n, nt;
    nt=MAX(NX,NY);
    n=MAX(nt,NZ/2);
    ip = alloc_1d_int(2 + (int) sqrt((double)n + 0.5));
    t = alloc_1d_double( 8 * nt );
    w = alloc_1d_double( n/2 + NZ/4 ) ; 
    ip[0] = 0;
  }
  if(SW_EQ==Qian_Sheng || SW_EQ == nematic_Allen_Cahn || SW_EQ == Olmsted_Goldbart){
    int n, nt;
    nt=MAX(N2X,N2Y);
    n=MAX(nt,N2Z/2);
    ip_extended = alloc_1d_int(2 + (int) sqrt((double)n + 0.5));
    t_extended = alloc_1d_double( 8 * nt );
    w_extended = alloc_1d_double( n/2 + N2Z/4 ) ; 
    ip_extended[0] = 0;
  }
  
  Calc_KX = Calc_KX_Ooura;
  Calc_KY = Calc_KY_Ooura;
  Calc_KZ = Calc_KZ_Ooura;
  Truncate_two_third_rule = Truncate_two_third_rule_ooura;
  Truncate_four_fifths_rule = Truncate_four_fifths_rule_ooura;
  Init_K();
  //if(SW_EQ==Qian_Sheng || SW_EQ == nematic_Allen_Cahn || SW_EQ == Olmsted_Goldbart){
  if( SW_EQ == nematic_Allen_Cahn ){
      Init_K_extended();
  }
}

void Init_fft(void){
  if(SW_FFT==Ooura){ // Ooura fft
    Init_fft_ooura();
  }else {
    fprintf(stderr,"specify SW_FFT correctly.\n");
    exit_job(EXIT_FAILURE);
  }
}

inline void A_k2dja_k_primitive(const Value &a
			 ,Value &da
			 ,const int &nx
			 ,const int &ny
			 ,const int &hnz_
			 ,const Value_int k_int
			 ,const double &wave_base
			 ){
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      for(int k=0; k<hnz_; k++){
	int k2=2*k;
	double wavenumber = k_int[i][j][k2] * wave_base;
	da[i][j][k2] = wavenumber * a[i][j][k2+1];
	da[i][j][k2+1] = -wavenumber * a[i][j][k2];
      }
    }
  }
}
void A_k2dxa_k(const Value &a, Value &da){
  A_k2dja_k_primitive(a,da,NX,NY,HNZ_,KX_int, WAVE_X);
}
void A_k2dya_k(const Value &a, Value &da){
  A_k2dja_k_primitive(a,da,NX,NY,HNZ_,KY_int, WAVE_Y);
}
void A_k2dza_k(const Value &a, Value &da){
  A_k2dja_k_primitive(a,da,NX,NY,HNZ_,KZ_int, WAVE_Z);
}
void A_k2dxa_k_extended(const Value &a, Value &da){
  A_k2dja_k_primitive(a,da,N2X,N2Y,HN2Z_,KX_int_extended, WAVE_X);
}
void A_k2dya_k_extended(const Value &a, Value &da){
  A_k2dja_k_primitive(a,da,N2X,N2Y,HN2Z_,KY_int_extended, WAVE_Y);
}
void A_k2dza_k_extended(const Value &a, Value &da){
  A_k2dja_k_primitive(a,da,N2X,N2Y,HN2Z_,KZ_int_extended, WAVE_Z);
}
void Omega_k2zeta_k(const Value omega[DIM], Value zeta[DIM-1]){
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ_; k++){
	if(KX_int[i][j][k] != 0){
	  zeta[0][i][j][k] = omega[1][i][j][k];
	  zeta[1][i][j][k] = omega[2][i][j][k];
	}else if(KY_int[i][j][k] != 0){
	  zeta[0][i][j][k] = omega[2][i][j][k];
	  zeta[1][i][j][k] = omega[0][i][j][k];
	}else{
	  zeta[0][i][j][k] = omega[0][i][j][k];
	  zeta[1][i][j][k] = omega[1][i][j][k];
	}
      }
    }
  }
  assert(zeta[0][0][0][0] == 0.0); 
  assert(zeta[1][0][0][0] == 0.0); 
}
void U_k2zeta_k(const Value u[DIM], Value zeta[DIM-1], double uk_dc[DIM]){
  for(int d=0; d<DIM; d++){
    uk_dc[d] = u[d][0][0][0]; 
  }
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<HNZ_; k++){
	int k2=2*k;
	double ks[DIM];
	ks[0] = KX_int[i][j][k2] * WAVE_X;
	ks[1] = KY_int[i][j][k2] * WAVE_Y;
	ks[2] = KZ_int[i][j][k2] * WAVE_Z;
	
	double omega_re[DIM];
	double omega_im[DIM];
	omega_re[0] = 
	  (ks[1] * u[2][i][j][k2+1] 
	    - ks[2] * u[1][i][j][k2+1]);
	omega_im[0] = 
	  -(ks[1] * u[2][i][j][k2] 
	   - ks[2] * u[1][i][j][k2]
	   );
	
	omega_re[1] = 
	  (ks[2] * u[0][i][j][k2+1] 
	    - ks[0] * u[2][i][j][k2+1]);
	omega_im[1] = 
	  -(ks[2] * u[0][i][j][k2] 
	   - ks[0] * u[2][i][j][k2]
	   );
	
	omega_re[2] = 
	  (ks[0] * u[1][i][j][k2+1] 
	    - ks[1] * u[0][i][j][k2+1]);
	omega_im[2] = 
	  -(ks[0] * u[1][i][j][k2] 
	   - ks[1] * u[0][i][j][k2]
	   );
	if(KX_int[i][j][k2] != 0){
	  zeta[0][i][j][k2] = omega_re[1];
	  zeta[0][i][j][k2+1] = omega_im[1];
	  zeta[1][i][j][k2] = omega_re[2];
	  zeta[1][i][j][k2+1] = omega_im[2];
	}else if(KY_int[i][j][k2] != 0){
	  zeta[0][i][j][k2] = omega_re[2];
	  zeta[0][i][j][k2+1] = omega_im[2];
	  zeta[1][i][j][k2] = omega_re[0];
	  zeta[1][i][j][k2+1] = omega_im[0];
	}else{
	  zeta[0][i][j][k2] = omega_re[0];
	  zeta[0][i][j][k2+1] = omega_im[0];
	  zeta[1][i][j][k2] = omega_re[1];
	  zeta[1][i][j][k2+1] = omega_im[1];
	}
      }
    }
  }
  assert(zeta[0][0][0][0] == 0.0); 
  assert(zeta[1][0][0][0] == 0.0); 
}
void Zeta_k2u_k(const Value zeta[DIM-1], const double uk_dc[DIM], Value u[DIM]){
  double omega_re[DIM];
  double omega_im[DIM];
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<HNZ_; k++){
	int k2=2*k;
	double dmy1_re, dmy2_re;
	double dmy1_im, dmy2_im;
	double kx = WAVE_X * KX_int[i][j][k2];
	double ky = WAVE_Y * KY_int[i][j][k2];
	double kz = WAVE_Z * KZ_int[i][j][k2];
	double ik2 = IK2[i][j][k2];

	dmy1_re =zeta[0][i][j][k2];
	dmy1_im =zeta[0][i][j][k2+1];
	dmy2_re =zeta[1][i][j][k2];
	dmy2_im =zeta[1][i][j][k2+1];

	if(KX_int[i][j][k2] != 0){
	  omega_re[0] = -(1./kx)*(ky * dmy1_re + kz * dmy2_re);
	  omega_im[0] = -(1./kx)*(ky * dmy1_im + kz * dmy2_im);
	  omega_re[1] = dmy1_re;
	  omega_im[1] = dmy1_im;
	  omega_re[2] = dmy2_re;
	  omega_im[2] = dmy2_im;
	  
	  
	}else if(KY_int[i][j][k2] != 0){
	  double dmy = -(kz/ky);
	  omega_re[0] = dmy2_re;
	  omega_im[0] = dmy2_im;
	  omega_re[1] = dmy * dmy1_re;
	  omega_im[1] = dmy * dmy1_im;
	  omega_re[2] = dmy1_re;
	  omega_im[2] = dmy1_im;
	}else{
	  omega_re[0] = dmy1_re;
	  omega_im[0] = dmy1_im;
	  omega_re[1] = dmy2_re;
	  omega_im[1] = dmy2_im;
	  omega_re[2] = 0.0;
	  omega_im[2] = 0.0;
	}
	kx *= ik2;
	ky *= ik2;
	kz *= ik2;
	u[0][i][j][k2] = 
	  (ky * omega_im[2] - kz * omega_im[1]); 
	u[0][i][j][k2+1] = 
	  -(ky * omega_re[2] - kz * omega_re[1]); 
	u[1][i][j][k2] = 
	  (kz * omega_im[0] - kx * omega_im[2]); 
	u[1][i][j][k2+1] = 
	  -(kz * omega_re[0] - kx * omega_re[2]); 
	u[2][i][j][k2] = 
	  (kx * omega_im[1] - ky * omega_im[0]); 
	u[2][i][j][k2+1] = 
	  -(kx * omega_re[1] - ky * omega_re[0]); 
      }
    }
  }
  for(int d=0; d<DIM; d++){
    u[d][0][0][0] = uk_dc[d];  
  }
}
void Zeta_k2omega_k(const Value zeta[DIM-1], Value omega[DIM]){
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<NZ_; k++){
	double dmy1, dmy2;
	double kx = WAVE_X * KX_int[i][j][k];
	double ky = WAVE_Y * KY_int[i][j][k];
	double kz = WAVE_Z * KZ_int[i][j][k];
	
	dmy1 =zeta[0][i][j][k];
	dmy2 =zeta[1][i][j][k];
		
	if(KX_int[i][j][k] != 0){
	  omega[0][i][j][k] = -(1./kx)*(ky * dmy1 + kz * dmy2);
	  omega[1][i][j][k] = dmy1;
	  omega[2][i][j][k] = dmy2;

	}else if(KY_int[i][j][k] != 0){
	  omega[0][i][j][k] = dmy2;
	  omega[1][i][j][k] = -(kz/ky)*dmy1;
	  omega[2][i][j][k] = dmy1;
	}else{
	  omega[0][i][j][k] = dmy1;
	  omega[1][i][j][k] = dmy2;
	  omega[2][i][j][k] = 0.0;
	}
      }
    }
  }
  for(int d=0; d<DIM; d++){
    omega[d][0][0][0] = 0.0; 
  }
}


void U_k2omega_k_inplace(Value u[DIM]){
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<HNZ_; k++){
	int k2=2*k;
	double ks[DIM];
	ks[0] = KX_int[i][j][k2] * WAVE_X;
	ks[1] = KY_int[i][j][k2] * WAVE_Y;
	ks[2] = KZ_int[i][j][k2] * WAVE_Z;

	double dmy_u_re[DIM];
	double dmy_u_im[DIM];
	double dmy_rot_re[DIM];
	double dmy_rot_im[DIM];
	for(int d=0; d<DIM; d++){
	  dmy_u_re[d] =u[d][i][j][k2];
	  dmy_u_im[d] =u[d][i][j][k2+1];
	}

	dmy_rot_re[0] = (ks[1] * dmy_u_im[2] - ks[2] * dmy_u_im[1]);
	dmy_rot_im[0] = -(ks[1] * dmy_u_re[2] - ks[2] * dmy_u_re[1]);

	dmy_rot_re[1] = (ks[2] * dmy_u_im[0] - ks[0] * dmy_u_im[2]);
	dmy_rot_im[1] = -(ks[2] * dmy_u_re[0] - ks[0] * dmy_u_re[2]);

	dmy_rot_re[2] = (ks[0] * dmy_u_im[1] - ks[1] * dmy_u_im[0]);
	dmy_rot_im[2] = -(ks[0] * dmy_u_re[1] - ks[1] * dmy_u_re[0]);

	for(int d=0; d<DIM; d++){
	  u[d][i][j][k2] = dmy_rot_re[d];
	  u[d][i][j][k2+1] = dmy_rot_im[d];
	}
      }
    }
  }
}

void Zeta_k2Strain_k(const Value zeta[DIM-1], Value strain_k[QDIM]){
  // Strain_k は 5成分
  double dmy[DIM]={0.,0.,0.};
  Zeta_k2u_k(zeta, dmy, strain_k);
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<HNZ_; k++){
	int k2=2*k;
	double ks[DIM];
	ks[0] = KX_int[i][j][k2] * WAVE_X * .5;
	ks[1] = KY_int[i][j][k2] * WAVE_Y * .5;
	ks[2] = KZ_int[i][j][k2] * WAVE_Z * .5;
	double u_dmy[DIM][2];
	for(int d=0;d<DIM;d++){
	  u_dmy[d][0] = u[d][i][j][k2];
	  u_dmy[d][1] = u[d][i][j][k2+1];
	}

	strain_k[0][i][j][k2] = ks[0] * u_dmy[0][1] * 2.;
	strain_k[0][i][j][k2+1] = -ks[0] * u_dmy[0][0] * 2.;

	strain_k[1][i][j][k2] = 
	  ks[0] * u_dmy[1][1]
	  +ks[1] * u_dmy[0][1];
	strain_k[1][i][j][k2+1] = 
	  - ks[0] * u_dmy[1][0] 
	  - ks[1] * u_dmy[0][0];

	strain_k[2][i][j][k2] = 
	  ks[0] * u_dmy[2][1]
	  +ks[2] * u_dmy[0][1];
	strain_k[2][i][j][k2+1] = 
	  - ks[0] * u_dmy[2][0] 
	  - ks[2] * u_dmy[0][0];

	strain_k[3][i][j][k2] = ks[1] * u_dmy[1][1] * 2.;
	strain_k[3][i][j][k2+1] = -ks[1] * u_dmy[1][0] * 2.;

	strain_k[4][i][j][k2] = 
	  ks[1] * u_dmy[2][1]
	  +ks[2] * u_dmy[1][1];
	strain_k[4][i][j][k2+1] = 
	  -ks[1] * u_dmy[2][0] 
	  -ks[2] * u_dmy[1][0];
      }
    }
  }
}
void U_k2divergence_k(const Value u[DIM], Value &div){
  for(int i=0; i<NX; i++){
    for(int j=0; j<NY; j++){
      for(int k=0; k<HNZ_; k++){
	int k2=2*k;
	double ks[DIM];
	ks[0] = KX_int[i][j][k2] * WAVE_X;
	ks[1] = KY_int[i][j][k2] * WAVE_Y;
	ks[2] = KZ_int[i][j][k2] * WAVE_Z;
	div[i][j][k2] = 0.0;
	div[i][j][k2+1] = 0.0;
	for(int d=0;d<DIM;d++){
	  div[i][j][k2] += ks[d] * u[d][i][j][k2+1];
	  div[i][j][k2+1] += -ks[d] * u[d][i][j][k2];
	}
      }
    }
  }
}
