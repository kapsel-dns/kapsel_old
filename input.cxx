//
// $Id: input.cxx,v 1.77 2006/06/08 05:56:03 nakayama Exp $
//
#include "input.h"

/////////////////////
/////////////////////
int Lubrication_measurement = 0;
int Terminal_velocity_measurement = 0;
int DKT_measurement = 0;
int Check_PB = 0;
int Fixed_particle = 0;
/////////////////////
/////////////////////

//////
EQ SW_EQ;
char *EQ_name[]={"Navier_Stokes"
		 ,"Qian_Sheng"
		 ,"nematic_Allen_Cahn"
		 ,"Olmsted_Goldbart"
		 ,"Slippy_Navier_Stokes"
		 ,"Shear_Navier_Stokes"
		 ,"Electrolyte"
		 ,"Two_fluid"
};
//////
SW_time SW_TIME;
//////
int SW_AVS;
char Out_dir[128];
char Out_name[128];;
int BINARY;
//////
int SW_UDF;
/////// FFT
//int SW_FFT=MPI_RFFTW; 
//int SW_FFT=RFFTW;
int SW_FFT=Ooura;
//////
/////// 計算条件の設定
int Nmax;
int Nmin;
int Ns[DIM];
int Ns_shear[DIM];
int HNs[DIM];
int TRNs[DIM];
int N2s[DIM];
int HN2s[DIM];
int TRNs_QS[DIM];
//////
int &NX=Ns[0];
int &NY=Ns[1];
int &NZ=Ns[2];
int &HNX=HNs[0];
int &HNY=HNs[1];
int &HNZ=HNs[2];
int &N2X=N2s[0];
int &N2Y=N2s[1];
int &N2Z=N2s[2];
int &HN2X=HN2s[0];
int &HN2Y=HN2s[1];
int &HN2Z=HN2s[2];
int &TRN_X=TRNs[0];
int &TRN_Y=TRNs[1];
int &TRN_Z=TRNs[2];
int &TRN_QS_X=TRNs_QS[0];
int &TRN_QS_Y=TRNs_QS[1];
int &TRN_QS_Z=TRNs_QS[2];
int HNZ_;
int NZ_;
int HN2Z_;
int N2Z_;
int ROTATION;
int HYDRO_int;
int STOKES;
WALL_type WALL;
double U_wall[DIM];
//const double A_radius_wall = 2.;
const double A_radius_wall = .5;
double Radius_wall;
int LJ_truncate;
Particle_IC DISTRIBUTION;
int N_iteration_init_distribution;
int FIX_CELL;
int FIX_CELLxyz[DIM];
//////
double EPSILON;
double T_LJ;
int LJ_powers;
int RESUMED;
int last_ts;
double Srate_depend_LJ_cap;
//////
double RHO;
double ETA;
double kBT;
double Shear_rate;
double Shear_strain;
int Shear_strain_int;
int NY_shear;
double LY_shear;
double dev_shear_stress[3];
double &dev_shear_stress_total = dev_shear_stress[0];
double &dev_shear_stress_lub = dev_shear_stress[1];
double &dev_shear_stress_lj = dev_shear_stress[2];
//////
double Delta_ETA;
double Nu_ratio;
//////
double Axel;
double DT;
//////
//////
double *MASS_RATIOS;
double *S_surfaces;// pretilt scalar order
double *W_surfaces;// spring cst. of anchoring
double DX;
double A_XI;
double A;
//////
double G;
int G_direction;
//////
int Component_Number;
int Particle_Number;
int *Particle_Numbers;
int GTS;
int Num_snap;

//
double NU;
double IRHO;
double *RHO_particle;
double *MASS;
//////
double *IMASS_RATIOS;
double *MOI;
double *IMASS;
double *IMOI;
//////
double Tdump;
double DT_noise;
//////
double DX3;
//////
double L[DIM];
double HL[DIM];
double L_particle[DIM];
double iL_particle[DIM];
double HL_particle[DIM];
double &LX=L[0];
double &LY=L[1];
double &LZ=L[2];
double &HLX=HL[0];
double &HLY=HL[1];
double &HLZ=HL[2];
double WAVE_X;
double WAVE_Y;
double WAVE_Z;
double KMAX2;
//////
double RADIUS;
double SIGMA;
double R_cutoff;
double XI;
double HXI;
double VF;
double VF_LJ;
double Ivolume;
//////
int MSTEP;
//////
double A_R_cutoff;
double LJ_dia;
/////// nematics
double A_LdG;//=-1.;
double B_LdG;//=3.*ABS(A_LdG);
double C_LdG;//=3.*ABS(A_LdG);
double L1;//=1.e-1;//1.;// in 0 or [1, \infty)
double L2;//=0;//1.e0; // in 0 or [1,\infty)
double L1tilde;
double L2tilde;
double Q2S_prefactor= DIM/(DIM-1.);
double IQ2S_prefactor=1./Q2S_prefactor;
double A_eq;// = (B_LdG + sqrt(SQ(B_LdG)-24.*A_LdG*C_LdG))/C_LdG*.25;
/////// Qian--Sheng viscosities
double Beta_1=1.;
double Beta_4=1.;
double Beta_5=1.;
double Beta_6=1.;
double Beta_56= (Beta_5 + Beta_6);
double Mu_1=1.;
double Mu_2=-1.;
double Mu_2overMu_1_half=Mu_2/Mu_1*.5;
double One_overMu1=1.0/Mu_1;
double Isotropic_viscosity=Beta_4-SQ(Mu_2)/(4.*Mu_1);
/////// Two_fluid
double Mean_Bulk_concentration;
int N_spec;
double Onsager_solute_coeff;
/////// Electrolyte
int Poisson_Boltzmann;
int External_field;
int AC;
double *Surface_charge;
double *Surface_charge_e;
double Elementary_charge=1.;
double Valency_counterion;
double Valency_positive_ion;
double Valency_negative_ion;
double Onsager_coeff_counterion;
double Onsager_coeff_positive_ion;
double Onsager_coeff_negative_ion;
double Dielectric_cst;
double Debye_length;
double E_ext[DIM];
double Frequency;
double Angular_Frequency;
//////
inline void Set_global_parameters(void){
  Particle_Number=0;
  for(int i=0; i<Component_Number; i++){
    Particle_Number += Particle_Numbers[i];
    IMASS_RATIOS[i] = 1.0/MASS_RATIOS[i];
    RHO_particle[i] = MASS_RATIOS[i] * RHO;
    MASS[i] = RHO_particle[i] * 4./3.*M_PI * POW3(RADIUS);
    IMASS[i] = 1./MASS[i];
    MOI[i] = 2./5. * MASS[i] * SQ(RADIUS);
    IMOI[i] = 1./MOI[i];
  }
  //  printf("%d\n",Particle_Number);
  IRHO=1./RHO;
  NU=ETA * IRHO;
  //////
  DX3= DX*DX*DX;
  //////
  HNZ_= NZ / 2 + 1;
  NZ_= 2 * HNZ_;
  for (int d=0;d<DIM;d++){
    HNs[d] = Ns[d]/2;
    TRNs[d] = (Ns[d]+2)/3;
    N2s[d] = Ns[d] * 2;
    HN2s[d] = Ns[d];
    TRNs_QS[d] = (N2s[d]+4)/5;
  }
  HN2Z_= N2s[2] / 2 + 1;
  N2Z_= 2 * HN2Z_;
  //////
  {
    for(int d=0;d<DIM;d++){
      L[d] = Ns[d]*DX; // real-dimension
      L_particle[d] = L[d];
      iL_particle[d] = 1./L_particle[d];
      HL[d] = L[d] * .5;
      HL_particle[d] = L[d] * .5;
    }
    //if(0){
    if( SW_EQ == Shear_Navier_Stokes){
      L_particle[1] = LY_shear;
      HL_particle[1] = LY_shear * .5;
    }
  }
  WAVE_X= PI2/LX;
  WAVE_Y= PI2/LY;
  WAVE_Z= PI2/LZ;
  if( SW_EQ == Navier_Stokes){
    KMAX2 = SQ(WAVE_X * TRN_X) 
      + SQ(WAVE_Y * TRN_Y) 
      + SQ(WAVE_Z * TRN_Z);
  }else if( SW_EQ == Qian_Sheng ){
    KMAX2 = SQ(WAVE_X * TRN_QS_X) 
      + SQ(WAVE_Y * TRN_QS_Y) 
      + SQ(WAVE_Z * TRN_QS_Z);
  }else if( SW_EQ == nematic_Allen_Cahn ){
    KMAX2 = SQ(WAVE_X * (HNX-1)) 
      + SQ(WAVE_Y * (HNY-1)) 
      + SQ(WAVE_Z * (HNZ-1));
  }else if( SW_EQ == Olmsted_Goldbart){
    KMAX2 = SQ(WAVE_X * TRN_QS_X) 
      + SQ(WAVE_Y * TRN_QS_Y) 
      + SQ(WAVE_Z * TRN_QS_Z);
  }else if( SW_EQ == Slippy_Navier_Stokes){
    KMAX2 = SQ(WAVE_X * TRN_X) 
      + SQ(WAVE_Y * TRN_Y) 
      + SQ(WAVE_Z * TRN_Z);
  }else if( SW_EQ == Shear_Navier_Stokes){
    KMAX2 = SQ(WAVE_X * TRN_X) 
      + SQ(WAVE_Y * TRN_Y) 
      + SQ(WAVE_Z * TRN_Z);
  }else if( SW_EQ == Electrolyte){
    KMAX2 = SQ(WAVE_X * TRN_X) 
      + SQ(WAVE_Y * TRN_Y) 
      + SQ(WAVE_Z * TRN_Z);
  }else if( SW_EQ == Two_fluid){
    KMAX2 = SQ(WAVE_X * TRN_X) 
      + SQ(WAVE_Y * TRN_Y) 
      + SQ(WAVE_Z * TRN_Z);
  }
  L1tilde=ABS(A)*L1*SQ(PI2)/KMAX2;
  L2tilde=ABS(A)*L2*SQ(PI2)/KMAX2;
  //////
  {
    if(SW_EQ == Slippy_Navier_Stokes ){
      double nu = NU;
      if(Delta_ETA < 0.){
	nu = IRHO * ETA;
      }else{
	nu = IRHO * ( ETA + Delta_ETA);
      }
      Tdump=1./(nu * KMAX2);
    }else if(SW_EQ == Navier_Stokes ){
      Tdump=1./(NU * KMAX2);
    }else if(SW_EQ == Shear_Navier_Stokes ){
      
      double shear_CFL_time = DX/(Shear_rate*LY_shear);
      //double shear_CFL_time = SQ(M_PI)/(Shear_rate*LY_shear);
      //double shear_stokes_time = RADIUS/(Shear_rate*LY_shear*0.5);
      double shear_stokes_time = XI/(Shear_rate*LY_shear*0.5);

      double mass_min = DBL_MAX;
      {
	for(int i=0; i<Component_Number; i++){
	  mass_min = MIN(mass_min, MASS[i]);
	}
      }
      double LJ_stokes_time = sqrt(mass_min * XI/Srate_depend_LJ_cap);
      Tdump=1./(NU * KMAX2);
      //fprintf(stderr, "aaaaa1:vis_time 2:shearCFLtime 3:shearstokestime 4:LJstokestime\n");
      //fprintf(stderr, "aaaaa%g %g %g %g\n" , Tdump, shear_CFL_time, shear_stokes_time, LJ_stokes_time);
      Tdump = MIN(Tdump, shear_CFL_time);
      Tdump = MIN(Tdump, shear_stokes_time);
    }else if(SW_EQ == Qian_Sheng){
      Tdump=1./(Isotropic_viscosity * ( SQ(WAVE_X * HNX)+SQ(WAVE_Y * HNY)+SQ(WAVE_Z * HNZ)));
      double rotational_diffusion_time 
	= Mu_1/fabs(A_LdG + L1tilde *KMAX2 + L1tilde *KMAX2);
      Tdump = MIN(Tdump, rotational_diffusion_time);
      Tdump *= 1.e-1;
    }else if(SW_EQ == nematic_Allen_Cahn){
      //fprintf(stderr, "%g %g\n" ,-A_LdG ,A_LdG + L1tilde *KMAX2 + L1tilde *KMAX2);
      //Axel = 5.e-1;
      double W=0;
      for(int i=0;i<Component_Number;i++){
	W = MAX(W_surfaces[i],W);
      } 
      double rotational_diffusion_time 
	//=MAX(fabs(A_LdG + L1tilde *KMAX2 + L1tilde *KMAX2),fabs(A_LdG));
	=MAX(fabs(A_LdG + L1tilde *KMAX2 + L1tilde *KMAX2+W/XI),fabs(A_LdG));
      rotational_diffusion_time 
	//= Mu_1/fabs(A_LdG + L1tilde *KMAX2 + L1tilde *KMAX2);
	= Mu_1/rotational_diffusion_time;
      Tdump = rotational_diffusion_time;
      //fprintf(stderr," aaa: %g %g\n", Tdump, rotational_diffusion_time);
    }else if(SW_EQ == Olmsted_Goldbart){
      Tdump=1./(Isotropic_viscosity * ( SQ(WAVE_X * HNX)+SQ(WAVE_Y * HNY)+SQ(WAVE_Z * HNZ)));
      double rotational_diffusion_time 
	= Mu_1/fabs(A_LdG + L1tilde *KMAX2 + L1tilde *KMAX2);
      Tdump = MIN(Tdump, rotational_diffusion_time);
    }else if(SW_EQ == Electrolyte){
      Tdump=1./(NU * KMAX2);
      double KMAX=sqrt(KMAX2);
      double dmy_onsager_coeff;
      if(N_spec==1){
	dmy_onsager_coeff=Onsager_coeff_counterion;
      }else if(N_spec==2){
	dmy_onsager_coeff=MAX(Onsager_coeff_positive_ion,Onsager_coeff_negative_ion);
      }
      double diffusion_time
	= 1./(kBT * dmy_onsager_coeff * KMAX2);
      Tdump = MIN(Tdump, diffusion_time);
      if(External_field){
	if(AC){
	  double dmy=1.e-2;
	  double frequency_time = 1./Frequency;
	  Angular_Frequency = PI2 * Frequency;
	  Tdump = MIN(Tdump, dmy*frequency_time);
	}
      }
    }else if(SW_EQ == Two_fluid ){
      Tdump=1./(NU * KMAX2);
      double Tdiff = 1./(kBT * Onsager_solute_coeff * KMAX2);
      //fprintf(stderr, "#aaa:%g %g\n", Tdump, Tdiff);
      Tdump = MIN(Tdump, Tdiff);
    }else{
      Tdump=1./(NU * KMAX2);
    }
    if(SW_TIME == AUTO){
      DT = Axel * Tdump;
      if(SW_EQ == Navier_Stokes ){
	if(kBT > 0){
	  if(1){
	    double sdv2 = (NX*NY*NZ)/POW3(DX) * ETA * kBT;
	      //* POW3(3./2.);
	    //DT_noise= 1.e3/sdv2;
	    DT_noise= 1.e0/sdv2;
	  }else{
	    double mass_max = 0.0;
	    for(int i=0; i<Component_Number; i++){
	      mass_max = MAX(mass_max, MASS[i]);
	    }
	    double themal_speed = kBT/mass_max;
	    DT_noise = XI/themal_speed;
	  }
	  //DT = Axel * MIN(Tdump, DT_noise);
	}
      }
      if(SW_EQ == Slippy_Navier_Stokes){
	DT *= 1.e-1;
	
	double rho_ratio_max=0.0; 
	for(int n=0;n<Component_Number;n++){
	  rho_ratio_max = MAX( rho_ratio_max, MASS_RATIOS[n]);
	}
	double vterm = 2./9. *SQ(RADIUS) *G/NU* (rho_ratio_max-1.);
	double dt_interface_stokes = XI/vterm;
	//    fprintf(stderr, "%g %g\n",dt_interface_stokes,DT);
      }
    }else if(SW_TIME == MANUAL){
      ;
    }
  }
  //////
  HXI = XI * 0.5;
  //////
  //////
  MSTEP= GTS * Num_snap;
  //////
  if(Lubrication_measurement){
    //LJ_dia = SIGMA;
    LJ_dia = SIGMA/pow(2.0,1./6.);
  }else if(DKT_measurement){
    //LJ_dia = SIGMA;
    LJ_dia = SIGMA/pow(2.0,1./6.);
  }else {
    if(SW_EQ == Shear_Navier_Stokes ){
      //LJ_dia = SIGMA;
      LJ_dia = MIN((SIGMA+XI)/pow(2.0,1./6.),SIGMA);
    }else {
      LJ_dia = SIGMA;
    }
  }
  R_cutoff = A_R_cutoff * LJ_dia;
  {
    double radius_dmy = pow(2.0,1./6.)*LJ_dia*.5;
    double lz;
    if(WALL == z_dirichlet){
      lz = LZ - 2.*Radius_wall;
    }else {
      lz = LZ;
    }
    if(SW_EQ == Shear_Navier_Stokes){
      Ivolume = 1./(LX * LY_shear * lz);
    }else {
      Ivolume = 1./(LX * LY * lz);
    }
    double dmy = (double)Particle_Number * 4./3.*M_PI * Ivolume;
    VF = dmy * POW3(RADIUS);
    VF_LJ = dmy * POW3(radius_dmy);
  }
}

inline void Set_LdG_parameters(UDFManager *ufin, UDFManager *ufout, Location &target){
  {
    ufin->get(target.sub("DX"),DX);
    target.down("Landau_de_Gennes_coefficients");
    ufin->get(target.sub("A"),A_LdG);
    ufin->get(target.sub("l1"),L1);
    ufin->get(target.sub("l2"),L2);
    B_LdG=3.*ABS(A_LdG);
    C_LdG=3.*ABS(A_LdG);
    //B_LdG=4.*ABS(A_LdG);
    //C_LdG=2.*ABS(A_LdG);
    A_eq = (B_LdG + sqrt(SQ(B_LdG)-24.*A_LdG*C_LdG))/C_LdG*.25;
  }
  {
    target.up();
    ufout->put(target.sub("DX"),DX);
    target.down("Landau_de_Gennes_coefficients");
    ufout->put(target.sub("A"),A_LdG);
    ufout->put(target.sub("l1"),L1);
    ufout->put(target.sub("l2"),L2);
  }
} 

UDFManager *ufin;
UDFManager *ufout;
UDFManager *ufsum;
UDFManager *ufres;
void Gourmet_file_io(const char *infile
		     ,char *outfile
		     ,char *sumfile
		     ,char *deffile
		     ,char *ctrlfile
		     ,char *resfile
		     ){

  if(file_check(infile)) ufin=new UDFManager(infile);

  // --------------------------------------------------------
  // レコードを追加するか新規にするか決めるため, 重複するけど
  // resumed or not は outfile開く前にも見とく
  {
    string str;
    ufin->get("resume.Calculation",str);
    if(str == "NEW"){
      RESUMED = 0;
    }else if(str == "CONTINUE"){
      RESUMED = 1;
    }else {
      fprintf(stderr, "invalid Calculation\n"); 
      exit_job(EXIT_FAILURE);
    }
  }
  if(!RESUMED){
    if(file_check(deffile)) ufout= new UDFManager(outfile,deffile,true);
  }
  else{
    if(file_check(deffile)) ufout= new UDFManager(outfile);
  }
  if(file_check(resfile)) ufres= new UDFManager(resfile);
  /*
  if(file_check(deffile)) ufout= new UDFManager(outfile,deffile,true);
  if(file_check(resfile)) ufres= new UDFManager(resfile);
  */
  // --------------------------------------------------------

  /////// resumed or not
  {
    string str;
    ufin->get("resume.Calculation",str);
    if(str == "NEW"){
      RESUMED = 0;
    }else if(str == "CONTINUE"){
      RESUMED = 1;
    }else {
      fprintf(stderr, "invalid Calculation\n"); 
      exit_job(EXIT_FAILURE);
    }
    ufout->put("Calculation",str);  
  }
  {
    ufin->get("resume.CONTINUE.Saved_Data.jikan.ts",last_ts);
    ufout->put("resume.CONTINUE.Saved_Data.jikan.ts",last_ts);
  }

  /////// select constitutive eq
  {
    Location target("constitutive_eq");
    string str;
    ufin->get(target.sub("type"),str);
    ufout->put(target.sub("type"),str);
    if(str == EQ_name[Navier_Stokes]){
      SW_EQ=Navier_Stokes;
      {
	target.down(EQ_name[SW_EQ]);
	{
	  ufin->get(target.sub("DX"),DX);
	  ufin->get(target.sub("RHO"),RHO);
	  ufin->get(target.sub("ETA"),ETA);
	  ufin->get(target.sub("kBT"),kBT);
	}
	{
	  ufout->put(target.sub("DX"),DX);
	  ufout->put(target.sub("RHO"),RHO);
	  ufout->put(target.sub("ETA"),ETA);
	  ufout->put(target.sub("kBT"),kBT);
	}
      }
    }else if(str == EQ_name[Qian_Sheng]){
      SW_EQ=Qian_Sheng;
      {
	target.down(EQ_name[SW_EQ]);
	Set_LdG_parameters(ufin, ufout, target);
      }
    }else if(str == EQ_name[nematic_Allen_Cahn]){
      SW_EQ=nematic_Allen_Cahn;
      {
	target.down(EQ_name[SW_EQ]);
	Set_LdG_parameters(ufin, ufout, target);
      }
    }else if(str == EQ_name[Olmsted_Goldbart]){
      SW_EQ=Olmsted_Goldbart;
      {
	target.down(EQ_name[SW_EQ]);
	Set_LdG_parameters(ufin, ufout, target);
      }
    }else if(str == EQ_name[Slippy_Navier_Stokes]){
      SW_EQ=Slippy_Navier_Stokes;
      {
	target.down(EQ_name[SW_EQ]);
	{
	  ufin->get(target.sub("DX"),DX);
	  ufin->get(target.sub("RHO"),RHO);
	  ufin->get(target.sub("ETA"),ETA);
	  ufin->get(target.sub("Nu_ratio"),Nu_ratio);
	  Delta_ETA = ETA * (Nu_ratio -1.);
	}
	{
	  ufout->put(target.sub("DX"),DX);
	  ufout->put(target.sub("RHO"),RHO);
	  ufout->put(target.sub("ETA"),ETA);
	  ufout->put(target.sub("Nu_ratio"),Nu_ratio);
	}
      }
    }else if(str == EQ_name[Shear_Navier_Stokes]){
      SW_EQ=Shear_Navier_Stokes;
      {
	target.down(EQ_name[SW_EQ]);
	{
	  ufin->get(target.sub("DX"),DX);
	  ufin->get(target.sub("RHO"),RHO);
	  ufin->get(target.sub("ETA"),ETA);
	  ufin->get(target.sub("kBT"),kBT);
	  ufin->get(target.sub("Shear_rate"),Shear_rate);
	}
	Shear_strain = 0.0;
	Shear_strain_int = 0;
	{
	  //Srate_depend_LJ_cap = 1.e3;
	  //Srate_depend_LJ_cap = Shear_rate * 1.e3;
	  Srate_depend_LJ_cap = DBL_MAX;
	}
	{
	  ufout->put(target.sub("DX"),DX);
	  ufout->put(target.sub("RHO"),RHO);
	  ufout->put(target.sub("ETA"),ETA);
	  ufout->put(target.sub("kBT"),kBT);
	  ufout->put(target.sub("Shear_rate"),Shear_rate);
	}
      }
    }else if(str == EQ_name[Electrolyte]){
      SW_EQ=Electrolyte;
      {
	target.down(EQ_name[Electrolyte]);
	{
	  ufin->get(target.sub("DX"),DX);
	  ufin->get(target.sub("RHO"),RHO);
	  ufin->get(target.sub("ETA"),ETA);
	  ufin->get(target.sub("kBT"),kBT);
	  ufin->get(target.sub("Dielectric_cst"),Dielectric_cst);

	  ufout->put(target.sub("DX"),DX);
	  ufout->put(target.sub("RHO"),RHO);
	  ufout->put(target.sub("ETA"),ETA);
	  ufout->put(target.sub("kBT"),kBT);
	  ufout->put(target.sub("Dielectric_cst"),Dielectric_cst);
	  {
	    ufin->get(target.sub("INIT_profile"),str);
	    ufout->put(target.sub("INIT_profile"),str);
	    if(str == "Uniform"){
	      Poisson_Boltzmann=0;
	    }else if(str == "Poisson_Boltzmann"){
	      Poisson_Boltzmann=1;
	    }else{
	      fprintf(stderr, "invalid INIT_profile\n"); 
	      exit_job(EXIT_FAILURE);
	    }
	  }
	  {
	    Location target("constitutive_eq.Electrolyte.Add_salt");
	    ufin->get(target.sub("type"),str);
	    ufout->put(target.sub("type"),str);
	    if(str == "saltfree"){
	      N_spec =1;
	      target.down("saltfree");
	      {
		ufin->get(target.sub("Valency_counterion"),Valency_counterion);
		ufout->put(target.sub("Valency_counterion"),Valency_counterion);
		ufin->get(target.sub("Onsager_coeff_counterion"),Onsager_coeff_counterion);
		ufout->put(target.sub("Onsager_coeff_counterion"),Onsager_coeff_counterion);
	      }
	    }else if(str == "salt"){
	      N_spec =2;
	      target.down("salt");
	      {
		ufin->get(target.sub("Valency_positive_ion"),Valency_positive_ion);
		ufout->put(target.sub("Valency_positive_ion"),Valency_positive_ion);
		ufin->get(target.sub("Valency_negative_ion"),Valency_negative_ion);
		ufout->put(target.sub("Valency_negative_ion"),Valency_negative_ion);
		ufin->get(target.sub("Onsager_coeff_positive_ion"),Onsager_coeff_positive_ion);
		ufout->put(target.sub("Onsager_coeff_positive_ion"),Onsager_coeff_positive_ion);
		ufin->get(target.sub("Onsager_coeff_negative_ion"),Onsager_coeff_negative_ion);
		ufout->put(target.sub("Onsager_coeff_negative_ion"),Onsager_coeff_negative_ion);
		ufin->get(target.sub("Debye_length"),Debye_length);
		ufout->put(target.sub("Debye_length"),Debye_length);
	      }
	      target.up();
	    }else{
	      fprintf(stderr, "invalid Add_salt\n"); 
	      exit_job(EXIT_FAILURE);
	    }
	  }
	  {
	    Location target("constitutive_eq.Electrolyte.Electric_field");
	    ufin->get(target.sub("type"),str);
	    ufout->put(target.sub("type"),str);
	    if(str == "OFF"){
	      External_field = 0;
	      for(int d=0;d<DIM;d++){
		E_ext[d] = 0.0;
	      }
	    }else if(str == "ON"){
	      External_field = 1;
	      target.down("ON");
	      {
		ufin->get(target.sub("type"),str);
		ufout->put(target.sub("type"),str);
		if(str == "DC"){
		  target.down("DC");
		  AC = 0;
		  ufin->get(target.sub("Ex"),E_ext[0]);
		  ufin->get(target.sub("Ey"),E_ext[1]);
		  ufin->get(target.sub("Ez"),E_ext[2]);
		  ufout->put(target.sub("Ex"),E_ext[0]);
		  ufout->put(target.sub("Ey"),E_ext[1]);
		  ufout->put(target.sub("Ez"),E_ext[2]);
		}else if(str == "AC"){
		  target.down("AC");
		  AC = 1;
		  ufin->get(target.sub("Ex"),E_ext[0]);
		  ufin->get(target.sub("Ey"),E_ext[1]);
		  ufin->get(target.sub("Ez"),E_ext[2]);
		  ufin->get(target.sub("Frequency"),Frequency);
		  ufout->put(target.sub("Ex"),E_ext[0]);
		  ufout->put(target.sub("Ey"),E_ext[1]);
		  ufout->put(target.sub("Ez"),E_ext[2]);
		  ufout->put(target.sub("Frequency"),Frequency);
		}else {
		  fprintf(stderr, "invalid switch for DC or AC\n"); 
		  exit_job(EXIT_FAILURE);
		}
	      }
	      target.up();
	    }else{
	      fprintf(stderr, "invalid Electric_field\n"); 
	      exit_job(EXIT_FAILURE);
	    }
	  }
	}
      }
    }else if(str == EQ_name[Two_fluid]){
      SW_EQ = Two_fluid;
      N_spec = 1;
      {
	target.down(EQ_name[Two_fluid]);
	{
	  ufin->get(target.sub("DX"),DX);
	  ufin->get(target.sub("RHO"),RHO);
	  ufin->get(target.sub("ETA"),ETA);
	  ufin->get(target.sub("kBT"),kBT);
	  ufin->get(target.sub("Mean_Bulk_concentration")
		    ,Mean_Bulk_concentration);
	  ufin->get(target.sub("Onsager_coeff")
		    ,Onsager_solute_coeff);
	  ufout->put(target.sub("DX"),DX);
	  ufout->put(target.sub("RHO"),RHO);
	  ufout->put(target.sub("ETA"),ETA);
	  ufout->put(target.sub("kBT"),kBT);
	  ufout->put(target.sub("Mean_Bulk_concentration")
		     ,Mean_Bulk_concentration);
	  ufout->put(target.sub("Onsager_coeff")
		     ,Onsager_solute_coeff);
	  {
	    ufin->get(target.sub("External_field"),str);
	    ufout->put(target.sub("External_field"),str);
	    if(str == "OFF"){
	      External_field=0;
	    }else if(str == "ON"){
	      External_field=1;
	      target.down("ON");
	      {
		ufin->get(target.sub("Ex"),E_ext[0]);
		ufin->get(target.sub("Ey"),E_ext[1]);
		ufin->get(target.sub("Ez"),E_ext[2]);
		ufout->put(target.sub("Ex"),E_ext[0]);
		ufout->put(target.sub("Ey"),E_ext[1]);
		ufout->put(target.sub("Ez"),E_ext[2]);
	      }
	      target.up();
	    }else{
	      fprintf(stderr, "invalid External_field\n"); 
	      exit_job(EXIT_FAILURE);
	    }
	  }
	}
      }
      if(0){
	fprintf(stderr, "%s: under construction\n"
		,EQ_name[Two_fluid]); 
	exit_job(EXIT_FAILURE);
      }
    }else{
      fprintf(stderr, "invalid constitutive_eq\n"); 
      exit_job(EXIT_FAILURE);
    }
    fprintf(stderr,"#\n# %s eq. selected.\n", EQ_name[SW_EQ]);
  }
  /////// 計算条件の設定
  Component_Number= ufin->size("Particle_spec[]");
  {
    MASS_RATIOS=alloc_1d_double(Component_Number);
    Particle_Numbers=alloc_1d_int(Component_Number);
    RHO_particle=alloc_1d_double(Component_Number);
    MASS=alloc_1d_double(Component_Number);
    IMASS=alloc_1d_double(Component_Number);
    IMASS_RATIOS=alloc_1d_double(Component_Number);
    MOI=alloc_1d_double(Component_Number);
    IMOI=alloc_1d_double(Component_Number);

    S_surfaces = alloc_1d_double(Component_Number);
    W_surfaces = alloc_1d_double(Component_Number);

    Surface_charge = alloc_1d_double(Component_Number);
    Surface_charge_e = alloc_1d_double(Component_Number);
    
  }
  {
    {
      fprintf(stderr, "#\n");
      if(SW_EQ == nematic_Allen_Cahn || 
	 SW_EQ == Olmsted_Goldbart || 
	 SW_EQ == Qian_Sheng ){ 
	int d=1;
	fprintf(stderr, "#%d:species",d++);
	fprintf(stderr, " %d:number_of_particle[i]",d++);
	fprintf(stderr, " %d:mass_density_ratio[i]",d++);
	fprintf(stderr, " %d:pretilt_scalar_order[i]",d++);
	fprintf(stderr, " %d:anchoring_strength[i]",d++);
      }else if(SW_EQ == Electrolyte){
	int d=1;
	fprintf(stderr, "#%d:species",d++);
	fprintf(stderr, " %d:number_of_particle[i]",d++);
	fprintf(stderr, " %d:mass_density_ratio[i]",d++);
	fprintf(stderr, " %d:Surface_charge[i]",d++);
      }else {
	int d=1;
	fprintf(stderr, "#%d:species",d++);
	fprintf(stderr, " %d:number_of_particle[i]",d++);
	fprintf(stderr, " %d:mass_density_ratio[i]",d++);
      }
      fprintf(stderr, "\n");
    }
    for(int i=0; i<Component_Number; i++){
      char str[256];
      sprintf(str,"Particle_spec[%d]",i);
      Location target(str);
      ufin->get(target.sub("Particle_number"),Particle_Numbers[i]);
      ufout->put(target.sub("Particle_number"),Particle_Numbers[i]);
      ufin->get(target.sub("MASS_RATIO"),MASS_RATIOS[i]);
      ufout->put(target.sub("MASS_RATIO"),MASS_RATIOS[i]);
      ufin->get(target.sub("Pretilt_scalar_order"),S_surfaces[i]);
      ufout->put(target.sub("Pretilt_scalar_order"),S_surfaces[i]);
      ufin->get(target.sub("Anchoring_coeff"),W_surfaces[i]);
      ufout->put(target.sub("Anchoring_coeff"),W_surfaces[i]);
      ufin->get(target.sub("Surface_charge"),Surface_charge[i]);
      ufout->put(target.sub("Surface_charge"),Surface_charge[i]);
      if(SW_EQ == nematic_Allen_Cahn || 
	 SW_EQ == Olmsted_Goldbart || 
	 SW_EQ == Qian_Sheng ){ 
	fprintf(stderr, "#%d %d %g %g %g\n"
		,i
		,Particle_Numbers[i]
		,MASS_RATIOS[i]
		,S_surfaces[i]
		,W_surfaces[i]
		);
      }else if(SW_EQ == Electrolyte){
	fprintf(stderr, "#%d %d %g %g\n"
		,i
		,Particle_Numbers[i]
		,MASS_RATIOS[i]
		,Surface_charge[i]
		);
      }else {
	fprintf(stderr, "#%d %d %f\n"
		,i
		,Particle_Numbers[i]
		,MASS_RATIOS[i]);
      }
    }
    fprintf(stderr, "#\n");
  }
  {
    ufin->get("A_XI",A_XI);
    ufout->put("A_XI",A_XI);
    XI= A_XI * DX; // surface thickness
    
    ufin->get("A",A);
    ufout->put("A",A);
    RADIUS= A*DX;
    SIGMA = 2.0 * RADIUS;
  }
  {
    Location target("gravity");
    ufin->get(target.sub("G"),G);
    string str;
    ufin->get(target.sub("G_direction"),str);
    if(str == "-X"){
      G_direction = 0;
    }else if(str == "-Y"){
      G_direction = 1;
    }else if(str == "-Z"){
      G_direction = 2;
    }else {
      fprintf(stderr, "invalid G_direction\n"); 
      exit_job(EXIT_FAILURE);
    }
    ufout->put(target.sub("G"),G);
    ufout->put(target.sub("G_direction"),str);
  }
  ufin->get("EPSILON",EPSILON);
  ufout->put("EPSILON",EPSILON);
  string str;
  ufin->get("LJ_powers",str);
  if(str == "12:6"){
    LJ_powers = 0;
  }else if(str == "24:12"){
    LJ_powers = 1;
  }else if(str == "36:18"){
    LJ_powers = 2;
  }else {
    fprintf(stderr, "invalid LJ_powers\n"); 
    exit_job(EXIT_FAILURE);
  }
  ufout->put("LJ_powers",str);  
  //  printf("%d\n",LJ_powers);
  {
    int np[DIM];
    Location target("mesh");
    ufin->get(target.sub("NPX"),np[0]);
    ufin->get(target.sub("NPY"),np[1]);
    ufin->get(target.sub("NPZ"),np[2]);
    ufout->put(target.sub("NPX"),np[0]);
    ufout->put(target.sub("NPY"),np[1]);
    ufout->put(target.sub("NPZ"),np[2]);
    NX = 1<<np[0];
    Ns_shear[0] = NX;
    NY = 1<<np[1];
    Ns_shear[1] = NY;
    if(SW_EQ == Shear_Navier_Stokes){
      NY_shear = 1<<(np[1]);
      LY_shear = (double)NY_shear * DX;
      NY = 1<<(np[1]+1);
      Ns_shear[1] = NY_shear;
    }
    NZ = 1<<np[2];
    Ns_shear[2] = NZ;
    Nmax = MAX(NX, MAX(NY, NZ));
    Nmin = MAX(NX, MAX(NY, NZ));
  }
  {
    Location target("time_increment");
    string str;
    ufin->get(target.sub("type"),str);
    ufout->put(target.sub("type"),str);
    if(str == "auto"){
      SW_TIME = AUTO;
      ufin->get(target.sub("auto.factor"),Axel);
      ufout->put(target.sub("auto.factor"),Axel);
    }else if(str == "manual"){
      SW_TIME = MANUAL;
      ufin->get(target.sub("manual.delta_t"),DT);
      ufout->put(target.sub("manual.delta_t"),DT);
    }else {
      fprintf(stderr, "invalid time_increment\n"); 
      exit_job(EXIT_FAILURE);
    }
  }
  {
    Location target("switch");
    string str;

    ufin->get(target.sub("ROTATION"),str);
    ufout->put(target.sub("ROTATION"),str);
    if(SW_EQ == Slippy_Navier_Stokes){
      ROTATION = 1;
    }else{
      if(str == "OFF"){
	ROTATION = 0;
      }else if(str == "ON"){
	ROTATION = 1;
      }else{
	fprintf(stderr, "invalid ROTATION\n"); 
	exit_job(EXIT_FAILURE);
      }
    }

    ufin->get(target.sub("HYDRO_int"),str);
    ufout->put(target.sub("HYDRO_int"),str);
    if(SW_EQ == Slippy_Navier_Stokes){
      HYDRO_int = 1;
    }else{
      if(str == "Correct"){
	HYDRO_int = 1;
      }else if(str == "free draining"){
	HYDRO_int = 0;
      }else if(str == "squeeze-lubrication and drain"){
	HYDRO_int = -1;
      }else{
	fprintf(stderr, "invalid HYDRO_int\n"); 
	exit_job(EXIT_FAILURE);
      }
    }

    ufin->get(target.sub("Stokes"),str);
    ufout->put(target.sub("Stokes"),str);
    if(str == "with advection"){
      STOKES  = 0;
    }else if(str == "w/o advection"){
      STOKES  = 1;
    }else{
      fprintf(stderr, "invalid switch.Stokes\n"); 
      exit_job(EXIT_FAILURE);
    }

    ufin->get(target.sub("LJ_truncate"),str);
    ufout->put(target.sub("LJ_truncate"),str);
    if(str == "ON"){
      LJ_truncate = 1;
    }else if(str == "OFF"){
      LJ_truncate = 0;
    }else if(str == "NONE"){
      LJ_truncate = -1;
    }else{
      fprintf(stderr, "invalid LJ_truncate\n"); 
      exit_job(EXIT_FAILURE);
    }
    if(LJ_truncate > 0){
      A_R_cutoff = pow(2.0,1./6.); //Lennard-Jones minimum;
    }else if(LJ_truncate == 0){
      const double max_A_R_cutoff = 2.5;
      A_R_cutoff = MIN(Nmin*DX*.5/SIGMA, max_A_R_cutoff);
    }else{
      A_R_cutoff = 0.;
    }

    {
      Location target("switch.INIT_distribution");
      string str;
      ufin->get(target.sub("type"),str);
      ufout->put(target.sub("type"),str);

      if(str == "NONE"){
	DISTRIBUTION = None;
      }else if(str == "uniform_random"){
	DISTRIBUTION = uniform_random;
      }else if(str == "random_walk"){
	DISTRIBUTION = random_walk;
	ufin->get(target.sub("random_walk.iteration"),N_iteration_init_distribution);
	ufout->put(target.sub("random_walk.iteration"),N_iteration_init_distribution);
      }else if(str == "FCC"){
	DISTRIBUTION = FCC;
      }else if(str == "BCC"){
	DISTRIBUTION = BCC;
      }else if(str == "user_specify"){
	DISTRIBUTION = user_specify;
      }else{
	cerr << str << endl;
	fprintf(stderr, "invalid DISTRIBUTION\n"); 
	exit_job(EXIT_FAILURE);
      }
    }

    {
      target.down("FIX_CELL");
      {
	char *xyz[DIM]={"x","y","z"};
	for(int d=0;d<DIM;d++){
	  ufin->get(target.sub(xyz[d]),str);
	  ufout->put(target.sub(xyz[d]),str);
	  if(str == "OFF"){
	    FIX_CELLxyz[d] = 0;
	  }else if(str == "ON"){
	    FIX_CELLxyz[d] = 1;
	  }else{
	    fprintf(stderr, "invalid FIX_CELL%s\n",xyz[d]); 
	    exit_job(EXIT_FAILURE);
	  }
	}
      }
      if(SW_EQ == Shear_Navier_Stokes){
	FIX_CELLxyz[1] = 1;
      }
      FIX_CELL = (FIX_CELLxyz[0] | FIX_CELLxyz[1] | FIX_CELLxyz[2]);
      target.up();
    }
  }
  {
    Location target("boundary_condition");
    string str;
    ufin->get(target.sub("type"),str);
    ufout->put(target.sub("type"),str);
    if(str == "full_periodic"){
      WALL = PBC;
      for(int d=0;d<DIM;d++){
	U_wall[d] = 0.0;
      }
    }else if(str == "z_dirichlet"){
      WALL = z_dirichlet;
      target.down("z_dirichlet");
      ufin->get(target.sub("wall_velocity_x"),U_wall[0]);
      ufout->put(target.sub("wall_velocity_x"),U_wall[0]);
      ufin->get(target.sub("wall_velocity_y"),U_wall[1]);
      ufout->put(target.sub("wall_velocity_y"),U_wall[1]);
      ufin->get(target.sub("wall_velocity_z"),U_wall[2]);
      ufout->put(target.sub("wall_velocity_z"),U_wall[2]);
      target.up();

      Radius_wall = A_radius_wall * DX;
    }else {
      fprintf(stderr, "invalid boundary_condition %s\n"
	      ,str.c_str()); 
      exit_job(EXIT_FAILURE);
    }
  }
  { /////// output;
    string str;
    Location target("output");
    ufin->get(target.sub("GTS"),GTS);
    ufout->put(target.sub("GTS"),GTS);
    ufin->get(target.sub("Num_snap"),Num_snap);
    ufout->put(target.sub("Num_snap"),Num_snap);
    { ////// AVS
      ufin->get(target.sub("AVS"),str);
      ufout->put(target.sub("AVS"),str);
      if(str == "ON"){
	SW_AVS = 1;
	target.down("ON");
	{
	  ufin->get(target.sub("Out_dir"),str);
	  ufout->put(target.sub("Out_dir"),str);
	  strcpy(Out_dir,str.c_str());
	  {
	    if (opendir(Out_dir)==NULL) {
	      fprintf(stderr,
		      "\tdirectory \"%s\",\"%s/avs\" does not seemed to exist.\n",
		      Out_dir,Out_dir);
	      fprintf(stderr,
		      "\texecute \"mkdir %s;mkdir %s/avs\" before run.\n",
		      Out_dir,Out_dir);
	      exit_job(EXIT_FAILURE);
	    }
	  }
	  ufin->get(target.sub("Out_name"),str);
	  ufout->put(target.sub("Out_name"),str);
	  strcpy(Out_name,str.c_str());
	  ufin->get(target.sub("FileType"),str);
	  ufout->put(target.sub("FileType"),str);
	  if(str == "BINARY"){
	    BINARY = 1;
	  }else if(str == "ASCII"){
	    BINARY = 0;
	  }else{
	    fprintf(stderr, "invalid FileType %s\n",str.c_str()); 
	    exit_job(EXIT_FAILURE);
	  }
	}
	target.up();
      }else if(str == "OFF"){
	SW_AVS = 0;
      }else{
	fprintf(stderr, "invalid switch for AVS\n"); 
	exit_job(EXIT_FAILURE);
      }
    }
    { ////// UDF
      ufin->get(target.sub("UDF"),str);
      ufout->put(target.sub("UDF"),str);
      if(str == "ON"){
	SW_UDF = 1;
      }else if(str == "OFF"){
	SW_UDF = 0;
      }else{
	fprintf(stderr, "invalid switch for UDF\n"); 
	exit_job(EXIT_FAILURE);
      }
    }
  }

  if((!RESUMED)
     &&(DISTRIBUTION != user_specify)
     ){
    delete ufin;
  }
  
  Set_global_parameters();

}

char *In_udf,*Sum_udf,*Out_udf,*Def_udf,*Ctrl_udf,*Res_udf;

//GOURMET上で与えられたファイル名を取得します

void file_get(const int argc, char *argv[]){

  const int Number_of_reuired_arguments = 5;

  if(argc < Number_of_reuired_arguments){
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "> %s -I[input UDF] -O[output UDF] -D[define UDF] -R[restart UDF]\n",
	   argv[0]);
    fprintf(stderr, "\n");
    exit_job(EXIT_FAILURE);
  }
  int R_selected = 0;
  In_udf=Sum_udf=Out_udf=Def_udf=Ctrl_udf=Def_udf=Res_udf=NULL;
  for(int i=1; i<argc; i++){
    char c=' ';
    char *p=argv[i];
    if(*p=='-' && *++p) c=*p++;
    switch(c){
    case 'I':   //インプットUDF
      In_udf=p;
      fprintf(stderr, "#using %s as input\n",p);
      break;
    case 'S':   //計算途中経過出力UDF
      Sum_udf=p;
      fprintf(stderr,"#using %s as summary\n",p);
      break;
    case 'O':   //アウトプットUDF
      Out_udf=p;
      fprintf(stderr,"#using %s as output\n",p);
      break;
    case 'D':   //定義UDF
      Def_udf=p;
      fprintf(stderr, "#using %s as definition\n",p);
      break;
    case 'M':   //制御用ファイル
      Ctrl_udf=p;
      fprintf(stderr,"#using %s as control\n",p);
      break;
    case 'R':   // リスタートUDF
      Res_udf=p;
      fprintf(stderr,"#using %s as restart\n",p);
      R_selected = 1;
      break;
    default:
      break;
    }
  }
  if((In_udf==NULL)
     ||(Out_udf==NULL)
     ||(Def_udf==NULL)
     ||(Res_udf==NULL)
     ){
    fprintf(stderr, "Program stopped because required udf file(s) is not given.\n");
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "> %s -I[input UDF] -O[output UDF] -D[define UDF] -R[restart UDF]\n",
	   argv[0]);
    fprintf(stderr, "\n");
    exit_job(EXIT_FAILURE);
  }

  if(R_selected){ // input.udf を restart.udf にコピーしとく
    FILE *fin, *fout;
    char s[256];
    if((fin=fopen(In_udf,"r"))==NULL){
      printf("cannot open file\n");
      exit(0);
    }
    if((fout=fopen(Res_udf,"w"))==NULL){
      printf("cannot open file\n");
      exit(0);
    }
    while(fgets(s,256,fin)!=NULL){
      fputs(s,fout);
    }
    fclose(fin);
    fclose(fout);
  }
  
}
