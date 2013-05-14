//
// $Id: input.cxx,v 1.3 2006/11/14 20:07:12 nakayama Exp $
//
#include "input.h"

/////////////////////
/////////////////////
int Fixed_particle = 0;
/////////////////////
/////////////////////

//////
EQ SW_EQ;
char *EQ_name[]={"Navier_Stokes"
		 ,"Shear_Navier_Stokes"
		 ,"Shear_Navier_Stokes_Lees_Edwards"
		 ,"Electrolyte"
};
//////
PT SW_PT;
char *PT_name[]={"spherical_particle"
		 ,"chain"
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
//int SW_FFT=Ooura;
int SW_FFT=IMKL_FFT;
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
int STOKES;
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
double ikBT;
double Shear_rate;
double Shear_rate_eff;
double Shear_strain_realized;
double Shear_strain;
double Shear_frequency;
double Inertia_stress;
int Shear_strain_int;
double dev_shear_stress[0];
double &dev_shear_stress_lj = dev_shear_stress[0];
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
int *Beads_Numbers;
int *Chain_Numbers;
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
/////// Two_fluid
double Mean_Bulk_concentration;
int N_spec;
double Onsager_solute_coeff;
/////// Electrolyte
int Poisson_Boltzmann;
int External_field;
int AC;
int Shear_AC;
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

//double kT_snap_v=0.59;
//double kT_snap_omega=0.61;
//double kT_snap_v=1.763;
double kT_snap_v=1.439;
//ldouble kT_snap_v=2.11;
double kT_snap_o=1.22;

double alpha_v;
double alpha_o;

//////shear_degree
double degree_oblique;

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
    if(kBT > 0.){ 
	ikBT=1./kBT;
    } 
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
    }
    WAVE_X= PI2/LX;
    WAVE_Y= PI2/LY;
    WAVE_Z= PI2/LZ;
    KMAX2 = SQ(WAVE_X * TRN_X) 
	+ SQ(WAVE_Y * TRN_Y) 
	+ SQ(WAVE_Z * TRN_Z);
    //////
    {
	if(SW_EQ == Navier_Stokes ){
	    Tdump=1./(NU * KMAX2);
	}else if ((SW_EQ == Shear_Navier_Stokes ) || (SW_EQ == Shear_Navier_Stokes_Lees_Edwards)){
	    double shear_CFL_time = DX/(Shear_rate*LY);
	    //double shear_stokes_time = RADIUS/(Shear_rate*LY*0.5);
	    double shear_stokes_time = XI/(Shear_rate*LY*0.5);
	    
	    double mass_min = DBL_MAX;
	    {
		for(int i=0; i<Component_Number; i++){
		    mass_min = MIN(mass_min, MASS[i]);
		}
	    }
	    double LJ_stokes_time = sqrt(mass_min * XI/Srate_depend_LJ_cap);
	    Tdump=1./(NU * KMAX2);
	    fprintf(stderr, "# aaaaa :vis_time 2:shearCFLtime 3:shearstokestime 4:LJstokestime\n");
	    fprintf(stderr, "# %g %g %g %g\n" , Tdump*Axel, shear_CFL_time, shear_stokes_time, LJ_stokes_time);
	    //Tdump = MIN(Tdump, shear_CFL_time);
	    //Tdump = MIN(Tdump, shear_stokes_time);
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
	}
	if(SW_TIME == AUTO){
	    DT = Axel * Tdump;
	    if(SW_EQ == Navier_Stokes ){
		if(kBT > 0){
		    if(1){
			double sdv2 = (NX*NY*NZ)/POW3(DX) * ETA * kBT;
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
    double dummy_pow;
    if((SW_EQ == Shear_Navier_Stokes ) || (SW_EQ == Shear_Navier_Stokes_Lees_Edwards)){
	if(LJ_powers == 0){
	    dummy_pow = pow(2.,1./6.);
	}	
	if(LJ_powers == 1){
	    dummy_pow = pow(2.,1./12.);
	}
	if(LJ_powers == 2){
	    dummy_pow = pow(2.,1./18.);
	}
	LJ_dia = MIN((SIGMA+XI)/dummy_pow, SIGMA);
    }else {
	LJ_dia = SIGMA;
    }
    R_cutoff = A_R_cutoff * LJ_dia;
    {
	double radius_dmy = dummy_pow*LJ_dia*.5;
	Ivolume = 1./(LX * LY * LZ);
	double dmy = (double)Particle_Number * 4./3.*M_PI * Ivolume;
	VF = dmy * POW3(RADIUS);
	VF_LJ = dmy * POW3(radius_dmy);
    }
}

UDFManager *ufin;
UDFManager *ufout;
//UDFManager *ufsum;
UDFManager *ufres;
void Gourmet_file_io(const char *infile
		     ,const char *outfile
		     ,const char *sumfile
		     ,const char *deffile
		     ,const char *ctrlfile
		     ,const char *resfile
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
	//if(file_check(deffile)) ufout= new UDFManager(outfile,deffile,true);
	if(file_check(deffile)) ufout= new UDFManager(outfile,deffile,false);
    }
    else{
	//if(file_check(deffile)) ufout= new UDFManager(outfile);
	if(file_check(deffile)) ufout= new UDFManager(outfile, 2);
    }
    //if(file_check(resfile)) ufres= new UDFManager(resfile);
    //if(file_check(deffile)) ufres= new UDFManager(resfile,2);
    if(file_check(deffile)) ufres= new UDFManager(resfile,deffile,false);
    /*
      if(file_check(deffile)) ufout= new UDFManager(outfile,deffile,true);
      if(file_check(resfile)) ufres= new UDFManager(resfile);
    */
    // --------------------------------------------------------
    
    /////// resumed or not
    {
	ufin->get("resume.CONTINUE.Saved_Data.jikan.ts",last_ts);
	//ufout->put("resume.CONTINUE.Saved_Data.jikan.ts",last_ts);
	//ufres->put("resume.CONTINUE.Saved_Data.jikan.ts",last_ts);
    }
    
    /////// select constitutive eq
    {
	Location target("constitutive_eq");
	string str;
	ufin->get(target.sub("type"),str);
	ufout->put(target.sub("type"),str);
	ufres->put(target.sub("type"),str);
	if(str == EQ_name[Navier_Stokes]){
	    SW_EQ=Navier_Stokes;
	    {
		target.down(EQ_name[SW_EQ]);
		{
		    ufin->get(target.sub("DX"),DX);
		    ufin->get(target.sub("RHO"),RHO);
		    ufin->get(target.sub("ETA"),ETA);
		    ufin->get(target.sub("kBT"),kBT);
		    ufin->get(target.sub("alpha_v"),alpha_v);
		    ufin->get(target.sub("alpha_o"),alpha_o);
		}
		{
		    ufout->put(target.sub("DX"),DX);
		    ufout->put(target.sub("RHO"),RHO);
		    ufout->put(target.sub("ETA"),ETA);
		    ufout->put(target.sub("kBT"),kBT);
		    ufout->put(target.sub("alpha_v"),alpha_v);
		    ufout->put(target.sub("alpha_o"),alpha_o);
		}
		{
		    ufres->put(target.sub("DX"),DX);
		    ufres->put(target.sub("RHO"),RHO);
		    ufres->put(target.sub("ETA"),ETA);
		    ufres->put(target.sub("kBT"),kBT);
		    ufres->put(target.sub("alpha_v"),alpha_v);
		    ufres->put(target.sub("alpha_o"),alpha_o);
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
		    ufin->get(target.sub("alpha_v"),alpha_v);
		    ufin->get(target.sub("alpha_o"),alpha_o);
		}
		Shear_strain_realized = 0.0;
		Shear_strain = 0.0;
		Shear_strain_int = 0;
		{
		    Srate_depend_LJ_cap = DBL_MAX;
		}
		{
		    ufout->put(target.sub("DX"),DX);
		    ufout->put(target.sub("RHO"),RHO);
		    ufout->put(target.sub("ETA"),ETA);
		    ufout->put(target.sub("kBT"),kBT);
		    ufout->put(target.sub("alpha_v"),alpha_v);
		    ufout->put(target.sub("alpha_o"),alpha_o);
		}
		{
		    ufres->put(target.sub("DX"),DX);
		    ufres->put(target.sub("RHO"),RHO);
		    ufres->put(target.sub("ETA"),ETA);
		    ufres->put(target.sub("kBT"),kBT);
		    ufres->put(target.sub("alpha_v"),alpha_v);
		    ufres->put(target.sub("alpha_o"),alpha_o);
		}
		{
		    Location target("constitutive_eq.Shear_Navier_Stokes.External_field");
		    ufin->get(target.sub("type"),str);
		    ufout->put(target.sub("type"),str);
		    ufres->put(target.sub("type"),str);
		    if(str == "DC"){
			target.down("DC");
			Shear_AC = 0;
			ufin->get(target.sub("Shear_rate"),Shear_rate);
			ufout->put(target.sub("Shear_rate"),Shear_rate);
			ufres->put(target.sub("Shear_rate"),Shear_rate);
			fprintf(stderr,"# DC steady shear: shear rate %f \n", Shear_rate);
		    }
		    if(str == "AC"){
			target.down("AC");
			Shear_AC = 1;
			ufin->get(target.sub("Shear_rate"),Shear_rate);
			ufout->put(target.sub("Shear_rate"),Shear_rate);
			ufres->put(target.sub("Shear_rate"),Shear_rate);
			
			ufin->get(target.sub("Frequency"),Shear_frequency);
			ufout->put(target.sub("Frequency"),Shear_frequency);
			ufres->put(target.sub("Frequency"),Shear_frequency);
			
			fprintf(stderr,"# AC oscillatory shear: (shear rate, frequency, the maximum amp of strain)= %f %f %f\n",Shear_rate, Shear_frequency, Shear_rate/Shear_frequency);
		    }
		}
	    }
	}else if(str == EQ_name[Shear_Navier_Stokes_Lees_Edwards]){
	    SW_EQ=Shear_Navier_Stokes_Lees_Edwards;
	    {
		target.down(EQ_name[SW_EQ]);
		{
		    ufin->get(target.sub("DX"),DX);
		    ufin->get(target.sub("RHO"),RHO);
		    ufin->get(target.sub("ETA"),ETA);
		    ufin->get(target.sub("kBT"),kBT);
		    ufin->get(target.sub("alpha_v"),alpha_v);
		    ufin->get(target.sub("alpha_o"),alpha_o);
		}
		Shear_strain_realized = 0.0;
		Shear_strain = 0.0;
		Shear_strain_int = 0;
		{
		    Srate_depend_LJ_cap = DBL_MAX;
		}
		{
		    ufout->put(target.sub("DX"),DX);
		    ufout->put(target.sub("RHO"),RHO);
		    ufout->put(target.sub("ETA"),ETA);
		    ufout->put(target.sub("kBT"),kBT);
		    ufout->put(target.sub("alpha_v"),alpha_v);
		    ufout->put(target.sub("alpha_o"),alpha_o);
		}
		{
		    ufres->put(target.sub("DX"),DX);
		    ufres->put(target.sub("RHO"),RHO);
		    ufres->put(target.sub("ETA"),ETA);
		    ufres->put(target.sub("kBT"),kBT);
		    ufres->put(target.sub("alpha_v"),alpha_v);
		    ufres->put(target.sub("alpha_o"),alpha_o);
		}
		{
		    Location target("constitutive_eq.Shear_Navier_Stokes_Lees_Edwards.External_field");
		    ufin->get(target.sub("type"),str);
		    ufout->put(target.sub("type"),str);
		    ufres->put(target.sub("type"),str);
		    if(str == "DC"){
			target.down("DC");
			Shear_AC = 0;
			ufin->get(target.sub("Shear_rate"),Shear_rate);
			ufout->put(target.sub("Shear_rate"),Shear_rate);
			ufres->put(target.sub("Shear_rate"),Shear_rate);
			fprintf(stderr,"# DC steady shear: shear rate %f \n", Shear_rate);
		    }
		    if(str == "AC"){// in near future, someone will extend this section.
			target.down("AC");
			Shear_AC = 0;
			ufin->get(target.sub("Shear_rate"),Shear_rate);
			ufout->put(target.sub("Shear_rate"),Shear_rate);
			ufres->put(target.sub("Shear_rate"),Shear_rate);
			fprintf(stderr,"# DC steady shear: shear rate %f \n", Shear_rate);
		    }
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
		    
		    {
			ufout->put(target.sub("DX"),DX);
			ufout->put(target.sub("RHO"),RHO);
			ufout->put(target.sub("ETA"),ETA);
			ufout->put(target.sub("kBT"),kBT);
			ufout->put(target.sub("Dielectric_cst"),Dielectric_cst);
		    }
		    {
			ufres->put(target.sub("DX"),DX);
			ufres->put(target.sub("RHO"),RHO);
			ufres->put(target.sub("ETA"),ETA);
			ufres->put(target.sub("kBT"),kBT);
			ufres->put(target.sub("Dielectric_cst"),Dielectric_cst);
		    }
		    {
			ufin->get(target.sub("INIT_profile"),str);
			ufout->put(target.sub("INIT_profile"),str);
			ufres->put(target.sub("INIT_profile"),str);
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
			ufres->put(target.sub("type"),str);
			if(str == "saltfree"){
			    N_spec =1;
			    target.down("saltfree");
			    {
				ufin->get(target.sub("Valency_counterion"),Valency_counterion);
				ufout->put(target.sub("Valency_counterion"),Valency_counterion);
				ufres->put(target.sub("Valency_counterion"),Valency_counterion);
				ufin->get(target.sub("Onsager_coeff_counterion"),Onsager_coeff_counterion);
				ufout->put(target.sub("Onsager_coeff_counterion"),Onsager_coeff_counterion);
				ufres->put(target.sub("Onsager_coeff_counterion"),Onsager_coeff_counterion);
			    }
			}else if(str == "salt"){
			    N_spec =2;
			    target.down("salt");
			    {
				{
				    ufin->get(target.sub("Valency_positive_ion"),Valency_positive_ion);
				    ufin->get(target.sub("Valency_negative_ion"),Valency_negative_ion);
				    ufin->get(target.sub("Onsager_coeff_positive_ion"),Onsager_coeff_positive_ion);
				    ufin->get(target.sub("Onsager_coeff_negative_ion"),Onsager_coeff_negative_ion);
				    ufin->get(target.sub("Debye_length"),Debye_length);
				}
				{
				    ufout->put(target.sub("Valency_positive_ion"),Valency_positive_ion);
				    ufout->put(target.sub("Valency_negative_ion"),Valency_negative_ion);
				    ufout->put(target.sub("Onsager_coeff_positive_ion"),Onsager_coeff_positive_ion);
				    ufout->put(target.sub("Onsager_coeff_negative_ion"),Onsager_coeff_negative_ion);
				    ufout->put(target.sub("Debye_length"),Debye_length);
				}
				{
				    ufres->put(target.sub("Valency_positive_ion"),Valency_positive_ion);
				    ufres->put(target.sub("Valency_negative_ion"),Valency_negative_ion);
				    ufres->put(target.sub("Onsager_coeff_positive_ion"),Onsager_coeff_positive_ion);
				    ufres->put(target.sub("Onsager_coeff_negative_ion"),Onsager_coeff_negative_ion);
				    ufres->put(target.sub("Debye_length"),Debye_length);
				}
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
			ufres->put(target.sub("type"),str);
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
				ufres->put(target.sub("type"),str);
				if(str == "DC"){
				    target.down("DC");
				    AC = 0;
				    {
					ufin->get(target.sub("Ex"),E_ext[0]);
					ufin->get(target.sub("Ey"),E_ext[1]);
					ufin->get(target.sub("Ez"),E_ext[2]);
				    }
				    {
					ufout->put(target.sub("Ex"),E_ext[0]);
					ufout->put(target.sub("Ey"),E_ext[1]);
					ufout->put(target.sub("Ez"),E_ext[2]);
				    }
				    {
					ufres->put(target.sub("Ex"),E_ext[0]);
					ufres->put(target.sub("Ey"),E_ext[1]);
					ufres->put(target.sub("Ez"),E_ext[2]);
				    }
				}else if(str == "AC"){
				    target.down("AC");
				    AC = 1;
				    {
					ufin->get(target.sub("Ex"),E_ext[0]);
					ufin->get(target.sub("Ey"),E_ext[1]);
					ufin->get(target.sub("Ez"),E_ext[2]);
					ufin->get(target.sub("Frequency"),Frequency);
				    }
				    {
					ufout->put(target.sub("Ex"),E_ext[0]);
					ufout->put(target.sub("Ey"),E_ext[1]);
					ufout->put(target.sub("Ez"),E_ext[2]);
					ufout->put(target.sub("Frequency"),Frequency);
				    }
				    {
					ufres->put(target.sub("Ex"),E_ext[0]);
					ufres->put(target.sub("Ey"),E_ext[1]);
					ufres->put(target.sub("Ez"),E_ext[2]);
					ufres->put(target.sub("Frequency"),Frequency);
				    }
				    
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
	}else{
	    fprintf(stderr, "invalid constitutive_eq\n"); 
	    exit_job(EXIT_FAILURE);
	}
	fprintf(stderr,"#\n# %s eq. selected.\n", EQ_name[SW_EQ]);
    }
    
    /////// 計算条件の設定
    {
	Location target("object_type");
	string str;
	ufin->get(target.sub("type"),str);
	ufout->put(target.sub("type"),str);
	ufres->put(target.sub("type"),str);
	if(str == PT_name[spherical_particle]){
	    SW_PT=spherical_particle;
	    {
		Component_Number= ufin->size("object_type.spherical_particle.Particle_spec[]");
		//fprintf(stderr, "# %s %d\n",str, Component_Number);
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
	    }
	}else if(str == PT_name[chain]){
	    SW_PT=chain;
	    Component_Number= ufin->size("object_type.chain.Chain_spec[]");
	    {
		MASS_RATIOS=alloc_1d_double(Component_Number);
		Particle_Numbers=alloc_1d_int(Component_Number);
		Beads_Numbers=alloc_1d_int(Component_Number);
		Chain_Numbers=alloc_1d_int(Component_Number);
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
	}
    }
    
    {
	{
	    fprintf(stderr, "#\n");
	    if(SW_PT == spherical_particle){
		if(SW_EQ == Electrolyte){
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
	    }else if(SW_PT == chain){
		if(SW_EQ == Electrolyte){
		    int d=1;
		    fprintf(stderr, "#%d:species",d++);
		    fprintf(stderr, " %d:total_number_of_particle[i]",d++);
		    fprintf(stderr, " %d:number_of_beads[i]",d++);
		    fprintf(stderr, " %d:number_of_chain[i]",d++);
		    fprintf(stderr, " %d:mass_density_ratio[i]",d++);
		    fprintf(stderr, " %d:Surface_charge[i]",d++);
		}else {
		    int d=1;
		    fprintf(stderr, "#%d:species",d++);
		    fprintf(stderr, " %d:total_number_of_particle[i]",d++);
		    fprintf(stderr, " %d:number_of_beads[i]",d++);
		    fprintf(stderr, " %d:number_of_chain[i]",d++);
		    fprintf(stderr, " %d:mass_density_ratio[i]",d++);
		}
	    }
	    fprintf(stderr, "\n");
	}
	
	if(SW_PT == spherical_particle){
	    for(int i=0; i<Component_Number; i++){
		char str[256];
		sprintf(str,"object_type.spherical_particle.Particle_spec[%d]",i);
		Location target(str);
		{
		    ufin->get(target.sub("Particle_number"),Particle_Numbers[i]);
		    ufin->get(target.sub("MASS_RATIO"),MASS_RATIOS[i]);
		    ufin->get(target.sub("Surface_charge"),Surface_charge[i]);
		}
		{
		    ufout->put(target.sub("Particle_number"),Particle_Numbers[i]);
		    ufout->put(target.sub("MASS_RATIO"),MASS_RATIOS[i]);
		    ufout->put(target.sub("Surface_charge"),Surface_charge[i]);
		}
		{
		    ufres->put(target.sub("Particle_number"),Particle_Numbers[i]);
		    ufres->put(target.sub("MASS_RATIO"),MASS_RATIOS[i]);
		    ufres->put(target.sub("Surface_charge"),Surface_charge[i]);
		}
		if(SW_EQ == Electrolyte){
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
		
		fprintf(stderr, "#\n");
		fprintf(stderr, "# Spherical Particles selected.\n");
	    }
	}else if(SW_PT == chain){
	    for(int i=0; i<Component_Number; i++){
		char str[256];
		sprintf(str,"object_type.chain.Chain_spec[%d]",i);
		Location target(str);
		{
		    ufin->get(target.sub("Beads_number"),Beads_Numbers[i]);
		    ufin->get(target.sub("Chain_number"),Chain_Numbers[i]);
		    ufin->get(target.sub("MASS_RATIO"),MASS_RATIOS[i]);
		    ufin->get(target.sub("Surface_charge"),Surface_charge[i]);
		}
		{
		    ufout->put(target.sub("Beads_number"),Beads_Numbers[i]);
		    ufout->put(target.sub("Chain_number"),Chain_Numbers[i]);
		    ufout->put(target.sub("MASS_RATIO"),MASS_RATIOS[i]);
		    ufout->put(target.sub("Surface_charge"),Surface_charge[i]);
		}
		{
		    ufres->put(target.sub("Beads_number"),Beads_Numbers[i]);
		    ufres->put(target.sub("Chain_number"),Chain_Numbers[i]);
		    ufres->put(target.sub("MASS_RATIO"),MASS_RATIOS[i]);
		    ufres->put(target.sub("Surface_charge"),Surface_charge[i]);
		}
		
		Particle_Numbers[i] = Beads_Numbers[i]*Chain_Numbers[i];
		
		if(SW_EQ == Electrolyte){
		    fprintf(stderr, "#%d %d %d %d %g %g\n"
			    ,i
			    ,Particle_Numbers[i]
			    ,Beads_Numbers[i]
			    ,Chain_Numbers[i]
			    ,MASS_RATIOS[i]
			    ,Surface_charge[i]
			);
		}else {
		    fprintf(stderr, "#%d %d %d %d %f\n"
			    ,i
			    ,Particle_Numbers[i]
			    ,Beads_Numbers[i]
			    ,Chain_Numbers[i]
			    ,MASS_RATIOS[i]);
		}
	    }
	    fprintf(stderr, "#\n");
	    fprintf(stderr, "# Flexible chains selected.\n");
	}
    }
    fprintf(stderr, "#\n");
    {
	ufin->get("A_XI",A_XI);
	ufout->put("A_XI",A_XI);
	ufres->put("A_XI",A_XI);
	XI= A_XI * DX; // surface thickness
	ufin->get("A",A);
	ufout->put("A",A);
	ufres->put("A",A);
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
	
	ufres->put(target.sub("G"),G);
	ufres->put(target.sub("G_direction"),str);
    }
    ufin->get("EPSILON",EPSILON);
    ufout->put("EPSILON",EPSILON);
    ufres->put("EPSILON",EPSILON);
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
    ufres->put("LJ_powers",str);  
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
	
	ufres->put(target.sub("NPX"),np[0]);
	ufres->put(target.sub("NPY"),np[1]);
	ufres->put(target.sub("NPZ"),np[2]);
	
	NX = 1<<np[0];
	Ns_shear[0] = NX;
	NY = 1<<np[1];
	Ns_shear[1] = NY;
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
	ufres->put(target.sub("type"),str);
	if(str == "auto"){
	    SW_TIME = AUTO;
	    ufin->get(target.sub("auto.factor"),Axel);
	    ufout->put(target.sub("auto.factor"),Axel);
	    ufres->put(target.sub("auto.factor"),Axel);
	}else if(str == "manual"){
	    SW_TIME = MANUAL;
	    ufin->get(target.sub("manual.delta_t"),DT);
	    ufout->put(target.sub("manual.delta_t"),DT);
	    ufres->put(target.sub("manual.delta_t"),DT);
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
	ufres->put(target.sub("ROTATION"),str);
	{
	    if(str == "OFF"){
		ROTATION = 0;
	    }else if(str == "ON"){
		ROTATION = 1;
	    }else{
		fprintf(stderr, "invalid ROTATION\n"); 
		exit_job(EXIT_FAILURE);
	    }
	}
	
	ufin->get(target.sub("Stokes"),str);
	ufout->put(target.sub("Stokes"),str);
	ufres->put(target.sub("Stokes"),str);
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
	ufres->put(target.sub("LJ_truncate"),str);
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
	    // A_R_cutoff = pow(2.0,1./6.); //Lennard-Jones minimum;
	    if(LJ_powers == 0){
		A_R_cutoff = pow(2.,1./6.);
	    }	
	    if(LJ_powers == 1){
		A_R_cutoff = pow(2.,1./12.);
		fprintf(stderr,"# A_R_cutoff %f\n", A_R_cutoff);
	    }
	    if(LJ_powers == 2){
		A_R_cutoff = pow(2.,1./18.);
		fprintf(stderr,"# A_R_cutoff %f\n", A_R_cutoff); 
	    }
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
	    ufres->put(target.sub("type"),str);
	    
	    if(str == "NONE"){
		DISTRIBUTION = None;
	    }else if(str == "uniform_random"){
		DISTRIBUTION = uniform_random;
	    }else if(str == "random_walk"){
		DISTRIBUTION = random_walk;
		ufin->get(target.sub("random_walk.iteration"),N_iteration_init_distribution);
		ufout->put(target.sub("random_walk.iteration"),N_iteration_init_distribution);
		ufres->put(target.sub("random_walk.iteration"),N_iteration_init_distribution);
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
		    ufres->put(target.sub(xyz[d]),str);
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
	    if((SW_EQ == Shear_Navier_Stokes) || (SW_EQ == Shear_Navier_Stokes_Lees_Edwards)){
		FIX_CELLxyz[1] = 1;
	    }
	    FIX_CELL = (FIX_CELLxyz[0] | FIX_CELLxyz[1] | FIX_CELLxyz[2]);
	    target.up();
	}
    }
    
    { /////// output;
	string str;
	Location target("output");
	ufin->get(target.sub("GTS"),GTS);
	ufin->get(target.sub("Num_snap"),Num_snap);
	
	ufout->put(target.sub("GTS"),GTS);
	ufout->put(target.sub("Num_snap"),Num_snap);
	
	ufres->put(target.sub("GTS"),GTS);
	ufres->put(target.sub("Num_snap"),Num_snap);
	{ ////// AVS
	    ufin->get(target.sub("AVS"),str);
	    ufout->put(target.sub("AVS"),str);
	    ufres->put(target.sub("AVS"),str);
	    if(str == "ON"){
		SW_AVS = 1;
		target.down("ON");
		{
		    ufin->get(target.sub("Out_dir"),str);
		    ufout->put(target.sub("Out_dir"),str);
		    ufres->put(target.sub("Out_dir"),str);
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
		    ufres->put(target.sub("Out_name"),str);
		    strcpy(Out_name,str.c_str());
		    ufin->get(target.sub("FileType"),str);
		    ufout->put(target.sub("FileType"),str);
		    ufres->put(target.sub("FileType"),str);
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
	    ufres->put(target.sub("UDF"),str);
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
    
    if(0){ // input.udf を restart.udf にコピーしとく
	//if(R_selected){ // input.udf を restart.udf にコピーしとく
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
