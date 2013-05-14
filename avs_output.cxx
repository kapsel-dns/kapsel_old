//
// $Id: avs_output.cxx,v 1.34 2006/05/15 09:44:04 kin Exp $
//

#include "avs_output.h"

const int Veclen = 5+5; 
const char *Label="ux uy uz phi pressure tau_xy tau_yz tau_zx tau_xx tau_yy"; // avs 出力ラベル用
const int Veclen_two_fluid = 5; 
const char *Label_two_fluid="ux uy uz phi concentration"; // avs 出力ラベル用
const int Veclen_QS = 4+5; 
//const char *Label_QS="ux uy uz phi q11 q12 q13 q22 q23"; // avs 出力ラベル用
const char *Label_QS="ux uy uz phi s1 s2 nx ny nz"; // avs 出力ラベル用
const int Veclen_charge = 4+3; 
const char *Label_charge="ux uy uz phi surface_charge rho e_potential"; // avs 出力ラベル用

//const AVS_Field Field = irregular;
const AVS_Field Field = uniform;
AVS_parameters Avs_parameters;
void Init_avs(const AVS_parameters &Avs_parameters){
  FILE *fout;
  fout=filecheckopen(Avs_parameters.fld_file,"w");
  fprintf(fout,"# AVS field file\n");
  fprintf(fout,"ndim=%d\n",DIM);
  fprintf(fout,"dim1=%d\n", Avs_parameters.nx);
  fprintf(fout,"dim2=%d\n", Avs_parameters.ny);
  fprintf(fout,"dim3=%d\n", Avs_parameters.nz);
  fprintf(fout,"nspace=%d\n", DIM);
  if(SW_EQ == Navier_Stokes 
     || SW_EQ == Slippy_Navier_Stokes
     || SW_EQ == Shear_Navier_Stokes 
     ){
    fprintf(fout,"veclen=%d\n", Veclen);
  }else if( SW_EQ == Two_fluid ){
    fprintf(fout,"veclen=%d\n", Veclen_two_fluid);
  }else if(SW_EQ==Qian_Sheng|| SW_EQ == nematic_Allen_Cahn || SW_EQ == Olmsted_Goldbart){
    fprintf(fout,"veclen=%d\n", Veclen_QS);
  }else if(SW_EQ==Electrolyte){
    fprintf(fout,"veclen=%d\n", Veclen_charge);
  }
  fprintf(fout,"data=float\n");
  if(Field == irregular ){
    fprintf(fout,"field=irregular\n");
  }else if(Field == uniform ){
    fprintf(fout,"field=uniform\n");
  }else {
    fprintf(stderr, "invalid Field\n"); 
    exit_job(EXIT_FAILURE);
  }
  fprintf(fout,"nstep=%d\n",Avs_parameters.nstep);
  if(SW_EQ == Navier_Stokes 
     || SW_EQ == Slippy_Navier_Stokes
     || SW_EQ == Shear_Navier_Stokes
     ){
    fprintf(fout,"label = %s\n", Label);
  }else if(SW_EQ == Two_fluid){
    fprintf(fout,"label = %s\n", Label_two_fluid);
  }else if(SW_EQ==Qian_Sheng|| SW_EQ == nematic_Allen_Cahn || SW_EQ == Olmsted_Goldbart){
    fprintf(fout,"label = %s\n", Label_QS);
  }else if(SW_EQ==Electrolyte){
    fprintf(fout,"label = %s\n", Label_charge);
  }
  
  fclose(fout);
  
  if(Field == irregular){
    fout=filecheckopen(Avs_parameters.cod_file,"wb");
    if(BINARY){
      for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
	for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend ; j++){
	  for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend ; k++){
	    float dmy= (float)(i*DX);
	    fwrite(&dmy,sizeof(float),1,fout);
	  }
	}
      }
      for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
	for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend ; j++){
	  for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend ; k++){
	    float dmy= (float)(j*DX);
	    fwrite(&dmy,sizeof(float),1,fout);
	  }
	}
      }
      for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
	for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend ; j++){
	  for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend ; k++){
	    float dmy= (float)(k*DX);
	    fwrite(&dmy,sizeof(float),1,fout);
	  }
	}
      }
    }else{
      fprintf(fout,"X Y Z\n");
      for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
	for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend ; j++){
	  for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend ; k++){
	  fprintf(fout,"%g %g %g\n",
		  (double)i*DX, (double)j*DX, (double)k*DX);
	  }
	}
      }
    }
    fclose(fout);
  }else if(Field == uniform){
    fout=filecheckopen(Avs_parameters.cod_file,"wb");
    {
      fprintf(fout,"%g\n%g\n%g\n%g\n%g\n%g\n"
	      ,(double)(Avs_parameters.istart)*DX
	      ,(double)(Avs_parameters.iend)*DX
	      ,(double)(Avs_parameters.jstart)*DX
	      ,(double)(Avs_parameters.jend)*DX
	      ,(double)(Avs_parameters.kstart)*DX
	      ,(double)(Avs_parameters.kend)*DX
	      );
    }
    fclose(fout);
  }else {
    fprintf(stderr, "invalid Field\n"); 
    exit_job(EXIT_FAILURE);
  }
}


void Set_avs_parameters(AVS_parameters &Avs_parameters){
  {
    Avs_parameters.nx=NX;
    if(SW_EQ == Shear_Navier_Stokes){
      Avs_parameters.ny=HNY;
    }else {
      Avs_parameters.ny=NY;
    }
    Avs_parameters.nz=NZ;
  }
  
  {
    Avs_parameters.istart = 0;
    Avs_parameters.iend = NX-1;
    Avs_parameters.jstart = 0;
    if(SW_EQ == Shear_Navier_Stokes){
      Avs_parameters.jend = HNY-1;
    }else {
      Avs_parameters.jend = NY-1;
    }
    Avs_parameters.kstart = 0;
    Avs_parameters.kend = NZ-1;
  }
  {
    sprintf(Avs_parameters.out_fld, "%s.fld", Out_name);
    {
      sprintf(Avs_parameters.out_cod, "%s.cod", Out_name);
      sprintf(Avs_parameters.cod_file, "%s/%s",
	      Out_dir, Avs_parameters.out_cod);
    }
    {
      sprintf(Avs_parameters.out_pfx, "avs/%s_", Out_name);
      sprintf(Avs_parameters.fld_file, "%s/%s",
	      Out_dir, Avs_parameters.out_fld);
    }

    sprintf(Avs_parameters.out_pfld, "%sp.fld", Out_name);
    //sprintf(Avs_parameters.out_pcod, "%sp.cod", Out_name);
    sprintf(Avs_parameters.out_ppfx, "avs/%sp_", Out_name);

    sprintf(Avs_parameters.pfld_file, "%s/%s",
	    Out_dir, Avs_parameters.out_pfld);
    //sprintf(Avs_parameters.pcod_file, "%s/%s",
    //    Out_dir, Avs_parameters.out_pcod);
  }
  {
    Avs_parameters.nstep=(Num_snap+1);
  }

}
inline void Binary_write(FILE *fout
			 ,AVS_parameters &Avs_parameters
			 ,Value &a
			 ){
  for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend; k++){
    for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend; j++){
      for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
	float dmy= (float)a[i][j][k];
	fwrite(&dmy,sizeof(float),1,fout);
      }
    }
  }
}
inline void Add_field_description(AVS_parameters &Avs_parameters
				  ,const CTime &time
				  ,const int &veclen
				  ){
  FILE *fout;
  char line[512];
  fout=filecheckopen(Avs_parameters.fld_file,"a");
  fprintf(fout,"time value = \"step%dtime%g\"\n"
	  ,time.ts, time.time);
  if(BINARY){
    static const int data_size=sizeof(float)*
      Avs_parameters.nx * Avs_parameters.ny * Avs_parameters.nz;
    if(Field == irregular){
      for(int n=0; n < DIM; n++){
	fprintf(fout,
		"coord %d file = %s filetype = binary skip = %d\n",
		n+1,
		Avs_parameters.out_cod, n*data_size);
      }
    }else if(Field == uniform){
      {
	for(int n=0; n < DIM; n++){
	  fprintf(fout,
		  "coord %d file = %s filetype = ascii skip = %d\n",
		  n+1, Avs_parameters.out_cod, n*2);
	}
      }
    }
    for(int n=0; n < veclen; n++){
      fprintf(fout,
	      "variable %d file = %s%d.dat filetype = binary skip = %d\n",
	      n+1,
	      Avs_parameters.out_pfx, time.ts, n * data_size);
    }
  }else{
    if(Field == irregular){
      for(int n=0; n < DIM; n++){
	fprintf(fout,
		"coord %d file = %s filetype = ascii skip = 1 offset = %d stride =%d\n",
		n+1, Avs_parameters.out_cod, n, DIM);
      }
    }else if(Field == uniform){
      for(int n=0; n < DIM; n++){
	fprintf(fout,
		"coord %d file = %s filetype = ascii skip = %d\n",
		n+1, Avs_parameters.out_cod, n*2);
      }
    }else {
      fprintf(stderr, "invalid Field\n"); 
      exit_job(EXIT_FAILURE);
    }
    for(int n=0; n < veclen; n++){
      fprintf(fout,
	      "variable %d file = %s%d.dat filetype = ascii skip = 2 offset = %d stride = %d\n",
	      n+1, Avs_parameters.out_pfx, time.ts, n, veclen);
    }
  }
  fprintf(fout,"EOT\n");
  fclose(fout);
}
void Output_avs(AVS_parameters &Avs_parameters
		,Value *zeta
		,double *uk_dc
		,Particle *p
		,const CTime &time){

  Value strain[QDIM]={f_particle[0]
		      ,f_particle[1]
		      ,f_particle[2]
		      ,f_ns0[0]
		      ,f_ns0[1]
  };
  {// strain tensor
    for(int d=0;d<DIM-1;d++){
      for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	  for(int k=0; k<NZ_; k++){
	    u[d][i][j][k] = ETA * zeta[d][i][j][k];
	  }
	}
      }
    }
    Zeta_k2Strain_k(u, strain);
    for(int d=0;d<QDIM;d++){
      A_k2a(strain[d]);
    }
  }

  Zeta_k2u(zeta, uk_dc, u);
  A_k2a(Pressure);
  {
    Reset_phi(phi);
    Make_phi_particle(phi, p);
  }

  Add_field_description(Avs_parameters,time, Veclen);

  FILE *fout;
  char line[512];
  sprintf(line,"timesteps=%d time=%f\n%s",time.ts,time.time, Label);
  sprintf(Avs_parameters.data_file,"%s/%s%d.dat",
	  Out_dir, Avs_parameters.out_pfx, time.ts);
  fout=filecheckopen(Avs_parameters.data_file,"wb");
  
  if(BINARY){
    Binary_write(fout, Avs_parameters, u[0]);
    Binary_write(fout, Avs_parameters, u[1]);
    Binary_write(fout, Avs_parameters, u[2]);
    Binary_write(fout, Avs_parameters, phi);
    Binary_write(fout, Avs_parameters, Pressure);
    {

      Binary_write(fout, Avs_parameters, strain[1]); // 12
      Binary_write(fout, Avs_parameters, strain[4]); // 23
      Binary_write(fout, Avs_parameters, strain[2]); // 13

      Binary_write(fout, Avs_parameters, strain[0]); // 11
      Binary_write(fout, Avs_parameters, strain[3]); // 22

    }
  }else{
    fprintf(fout,"%s\n", line);
    for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend; k++){
      for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend; j++){
	for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
	  fprintf(fout,"%.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g\n"
		  ,u[0][i][j][k]
		  ,u[1][i][j][k]
		  ,u[2][i][j][k]
		  ,phi[i][j][k]
		  ,Pressure[i][j][k]
		  ,strain[1][i][j][k]
		  ,strain[4][i][j][k]
		  ,strain[2][i][j][k]
		  ,strain[0][i][j][k]
		  ,strain[3][i][j][k]
		  );
	}
      }
    }
  }
  fclose(fout);

  if(0){ // ux(y) at x=Lx/2, z= Lz/2
    for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend ; j++){
      int i=HNX;
      int k=HNZ;
      fprintf(stdout,"%g %g\n",(double)j*DX, u[0][i][j][k]);
    }
    fprintf(stdout,"\n\n");
  }
}
void Output_avs_two_fluid(AVS_parameters &Avs_parameters
			  ,Value zeta[DIM-1]
			  ,double uk_dc[DIM]
			  ,Value u_r[DIM] // working memory
			  ,Value *conc_k
			  ,Value *conc_r // working memory
			  ,Particle *p
			  ,Value &phi // working memory
			  ,const CTime &time){
  Zeta_k2u(zeta, uk_dc, u_r);
  
  {
    Reset_phi(phi);
    Make_phi_particle(phi, p);
  }
  {
    A_k2a_out(conc_k[0], conc_r[0]);
    if(0){
      for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
	for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend ; j++){
	  for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend ; k++){
	    conc_r[0][i][j][k] *= (1.-phi[i][j][k]);
	  }
	}
      }
    }
  }

  Add_field_description(Avs_parameters,time, Veclen_two_fluid);

  FILE *fout;
  char line[512];
  sprintf(line,"timesteps=%d time=%f\n%s"
	  ,time.ts,time.time, Label_two_fluid);
  sprintf(Avs_parameters.data_file,"%s/%s%d.dat",
	  Out_dir, Avs_parameters.out_pfx, time.ts);
  fout=filecheckopen(Avs_parameters.data_file,"wb");
  
  if(BINARY){
    Binary_write(fout, Avs_parameters, u[0]);
    Binary_write(fout, Avs_parameters, u[1]);
    Binary_write(fout, Avs_parameters, u[2]);
    Binary_write(fout, Avs_parameters, phi);
    Binary_write(fout, Avs_parameters, conc_r[0]);
  }else{
    fprintf(fout,"%s\n", line);
    for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend; k++){
      for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend; j++){
	for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
	  fprintf(fout,"%.3g %.3g %.3g %.3g %.3g\n"
		  ,u[0][i][j][k]
		  ,u[1][i][j][k]
		  ,u[2][i][j][k]
		  ,phi[i][j][k]
		  ,conc_r[0][i][j][k]
		  );
	}
      }
    }
  }
  fclose(fout);

  if(1){
    fprintf(stderr, "%d\n",time.ts);

    int i=HNX;
    int j=HNY;
    //for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
    //for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend ; j++){
	for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend ; k++){
	  fprintf(stdout, "%g %g\n", conc_r[0][i][j][k]
		  ,phi[i][j][k]);
	}
	fprintf(stdout, "\n\n");
  }
}
void Output_udf(UDFManager *ufout
		,AVS_parameters &Avs_parameters
		,Value *zeta
		,double *uk_dc
		,const Particle *p
		,const CTime &time
		){
  ufout->newRecord();
  ufout->put("E",1.0);
  ufout->put("t",time.ts);
  for(int j=0;j<Particle_Number;j++){
    char str[256];
    sprintf(str,"Particles[%d]",j);
    Location target(str);
    ufout->put(target.sub("R.x"),p[j].x[0]);
    ufout->put(target.sub("R.y"),p[j].x[1]);
    ufout->put(target.sub("R.z"),p[j].x[2]);
    ufout->put(target.sub("v.x"),p[j].v[0]);
    ufout->put(target.sub("v.y"),p[j].v[1]);
    ufout->put(target.sub("v.z"),p[j].v[2]);
  }
}
void Output_avs_QS(AVS_parameters &Avs_parameters
		   ,Value *zeta
		   ,double *uk_dc
		   ,Value *tensor_order
		   ,Particle *p
		   ,const CTime &time
		   ){
  static double scale = POW3(.5);
  static double scale2 = SQ(scale);

  static const double s1_coeff = sqrt(1.5);
  static const double s2_coeff = sqrt(0.5);

  Zeta_k2u(zeta, uk_dc, u);
  
  {
    Reset_phi(phi);
    Make_phi_particle(phi, p);
  }
  {
    if(1){
      Q_k2q_out(tensor_order, tmp1_QDIM);
    }else{
      for(int d=0;d<QDIM;d++){
	Truncate_four_seconds_rule(tensor_order[d]);
	for(int i=0; i<NX; i++){
	  for(int j=0; j<NY; j++){
	    for(int k=0; k<NZ_; k++){
	      tmp1_QDIM[d][i][j][k]=tensor_order[d][i][j][k];
	    }
	  }
	}
	K_space_extension(tmp1_QDIM[d]);
	A_k2a_extended(tmp1_QDIM[d]);
      }
    }


    {
      double qq_trace_max = 0.0;
      double qq_trace_min = DBL_MAX;
      //for(int i=0; i<N2X; i++){
      //for(int j=0; j<N2Y; j++){
      //  for(int k=0; k<N2Z; k++){
      for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	  for(int k=0; k<NZ; k++){
	    double q_full[DIM][DIM];
	    Q2q_full(tmp1_QDIM, i, j, k, q_full);
	    double qqij[DIM][DIM];
	    AdotB(q_full, q_full, qqij);
	    double qq_trace= Q2qq(q_full);
	    double qqqij[DIM][DIM];
	    AdotB(q_full, qqij, qqqij);
	    double qqq_trace= 0.0;
	    for(int m=0;m<DIM;m++){
	      qqq_trace += qqqij[m][m];
	    }


	    double intensity = sqrt(qq_trace * scale2);
	    double chi;
	    if(intensity > 0.0){
	      chi = One_third * acos(Root_six * qqq_trace/POW3(intensity));
	    }else{
	      chi = 0.0;
	    }
	    
	    qq_trace_min = MIN(qq_trace_min,qq_trace);
	    qq_trace_max = MAX(qq_trace_max,qq_trace);

	    
	    double phi = atan2(q_full[1][2], q_full[0][2]);
	    double theta = (atan2(q_full[0][1], sin(phi)* q_full[0][2]));
	    
	    {
	      double sdmy = sqrt(qq_trace * scale2 * Q2S_prefactor)/A_eq;
	      tmp1_QDIM[0][i][j][k] = sdmy;
	      //tmp1_QDIM[0][i][j][k] = intensity * cos(chi) * s1_coeff/A_eq;
	      tmp1_QDIM[1][i][j][k] = intensity * cos(chi) * s1_coeff/A_eq;
	      //tmp1_QDIM[1][i][j][k] = intensity * sin(chi) * s2_coeff;
	      //tmp1_QDIM[2][i][j][k] = theta;
	      //tmp1_QDIM[3][i][j][k] = phi;
	      //tmp1_QDIM[4][i][j][k] = qqq_trace;
	      {
		double sin_theta = sin(theta);
		double xdmy = sin_theta * cos(phi) * sdmy;
		double ydmy = sin_theta * sin(phi) * sdmy;
		double zdmy = cos(theta) * sdmy;
		if(xdmy < 0.0){
		  tmp1_QDIM[2][i][j][k] = -xdmy;
		  tmp1_QDIM[3][i][j][k] = -ydmy;
		  tmp1_QDIM[4][i][j][k] = -zdmy;
		}else {
		  tmp1_QDIM[2][i][j][k] = xdmy;
		  tmp1_QDIM[3][i][j][k] = ydmy;
		  tmp1_QDIM[4][i][j][k] = zdmy;
		}
	      }
	    }
	  }
	}
      }
      if(0){
	fprintf(stdout, "#Output_avs_QS: (%g %g) (%g %g) (%g %g)\n"
		,qq_trace_max
		,qq_trace_min
		,qq_trace_max * scale2
		,qq_trace_min * scale2
		, sqrt(qq_trace_max * scale2 *Q2S_prefactor)/A_eq
		, sqrt(qq_trace_min * scale2 *Q2S_prefactor)/A_eq
		);
      }
    }
  }

  Add_field_description(Avs_parameters,time, Veclen_QS);

  FILE *fout;
  char line[512];
  sprintf(line,"timesteps=%d time=%f\n%s",time.ts,time.time, Label_QS);
  sprintf(Avs_parameters.data_file,"%s/%s%d.dat",
	  Out_dir, Avs_parameters.out_pfx, time.ts);
  fout=filecheckopen(Avs_parameters.data_file,"wb");
  
  if(BINARY){
    Binary_write(fout, Avs_parameters, u[0]);
    Binary_write(fout, Avs_parameters, u[1]);
    Binary_write(fout, Avs_parameters, u[2]);
    Binary_write(fout, Avs_parameters, phi);
    Binary_write(fout, Avs_parameters, tmp1_QDIM[0]);
    Binary_write(fout, Avs_parameters, tmp1_QDIM[1]);
    Binary_write(fout, Avs_parameters, tmp1_QDIM[2]);
    Binary_write(fout, Avs_parameters, tmp1_QDIM[3]);
    Binary_write(fout, Avs_parameters, tmp1_QDIM[4]);
  }else{
    fprintf(fout,"%s\n", line);
    for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend; k++){
      for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend; j++){
	for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
	  fprintf(fout,"%.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g\n"
 		  ,u[0][i][j][k],u[1][i][j][k],u[2][i][j][k]
 		  ,phi[i][j][k]
 		  ,tmp1_QDIM[0][i][j][k]
 		  ,tmp1_QDIM[1][i][j][k]
 		  ,tmp1_QDIM[2][i][j][k]
 		  ,tmp1_QDIM[3][i][j][k]
 		  ,tmp1_QDIM[4][i][j][k]
		  );
	}
      }
    }
  }
  fclose(fout);
}

void Output_avs_charge(AVS_parameters &Avs_parameters
		       ,Value *zeta
		       ,double *uk_dc
		       ,Value *Concentration
		       ,Particle *p
		       ,const CTime &time
		       ){
  Zeta_k2u(zeta, uk_dc, u);
  
  Value potential = f_particle[0];
  Value dmy_value0 = f_particle[1];
  {
    Conc_k2charge_field(p, Concentration, potential, phi, dmy_value0);
    A2a_k(potential);
    Charge_field_k2Coulomb_potential_k_PBC(potential);
    A_k2a(potential);
    for(int n=0;n<N_spec;n++){
      A_k2a(Concentration[n]);
    }
  }
  {
    Reset_phi(phi);
    Reset_phi(up[0]);
    Make_phi_qq_particle(phi, up[0], p);
  }

  Add_field_description(Avs_parameters,time, Veclen_charge);

  FILE *fout;
  char line[512];
  sprintf(line,"timesteps=%d time=%f\n%s",time.ts,time.time, Label_charge);
  sprintf(Avs_parameters.data_file,"%s/%s%d.dat",
	  Out_dir, Avs_parameters.out_pfx, time.ts);
  fout=filecheckopen(Avs_parameters.data_file,"wb");
  
  double dmy_surface_area = PI4 * RADIUS * RADIUS;
  if(BINARY){
    Binary_write(fout, Avs_parameters, u[0]);
    Binary_write(fout, Avs_parameters, u[1]);
    Binary_write(fout, Avs_parameters, u[2]);
    Binary_write(fout, Avs_parameters, phi);
    for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend ; k++){
      for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend ; j++){
	for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
	  float dmy = up[0][i][j][k] * dmy_surface_area;
	  fwrite(&dmy,sizeof(float),1,fout);
	}
      }
    }
    for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend; k++){
      for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend; j++){
	for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
	  float dmy=0.;
	  for(int n=0; n<N_spec; n++){
	    dmy += (float)(Elementary_charge*Valency[n]*Concentration[n][i][j][k]);
	  }
	  dmy *= (1.-phi[i][j][k]);
	  fwrite(&dmy,sizeof(float),1,fout);
	}
      }
    }
    Binary_write(fout, Avs_parameters, potential);
  }else{
    fprintf(fout,"%s\n", line);
    for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend; k++){
      for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend; j++){
	for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
	  double dmy = 0.;
	  for(int n=0; n<N_spec; n++){
	    dmy += Elementary_charge*Valency[n]*Concentration[n][i][j][k];
	  }
	  fprintf(fout,"%.3g %.3g %.3g %.3g %.3g %.3g %.3g\n"
		  ,u[0][i][j][k],u[1][i][j][k],u[2][i][j][k]
		  ,phi[i][j][k]
		  ,up[0][i][j][k]*dmy_surface_area
		  ,(1.-phi[i][j][k])*dmy
		  ,potential[i][j][k]
		  );
	}
      }
    }
  }
  fclose(fout);

  for(int n=0; n<N_spec; n++){
    A2a_k(Concentration[n]);
  }
}
