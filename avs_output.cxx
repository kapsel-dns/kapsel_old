/*!
  \file avs_output.cxx
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief Output routines for field data in AVS/Express format
 */

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
     || SW_EQ == Shear_Navier_Stokes || SW_EQ == Shear_Navier_Stokes_Lees_Edwards
     ){
    fprintf(fout,"veclen=%d\n", Veclen);
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
     || SW_EQ == Shear_Navier_Stokes || SW_EQ == Shear_Navier_Stokes_Lees_Edwards
     ){
    fprintf(fout,"label = %s\n", Label);
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
    Avs_parameters.ny=NY;
    Avs_parameters.nz=NZ;
  }
  
  {
    Avs_parameters.istart = 0;
    Avs_parameters.iend = NX-1;
    Avs_parameters.jstart = 0;
    Avs_parameters.jend = NY-1;
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
			 ,double *a
			 ){
  int im;
  for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend; k++){
    for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend; j++){
      for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
	im = (i * NY * NZ_) + (j * NZ_) + k;
	float dmy= (float)a[im];
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
		,double **zeta
		,double *uk_dc
		,Particle *p
		,const CTime &time){

  double *strain[QDIM]={f_particle[0]
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
		int im=(i*NY*NZ_)+(j*NZ_)+k;
	    u[d][im] = ETA * zeta[d][im];
	  }
	}
      }
    }
    Zeta_k2Strain_k(u, strain);
    for(int d=0;d<QDIM;d++){
      A_k2a(strain[d]);
    }
  }

  if (SW_EQ == Shear_Navier_Stokes_Lees_Edwards) {
      Zeta_k2u_k_OBL(zeta, uk_dc, u);
      U_k2u(u);
      U_oblique2u(u);
  } else {
      Zeta_k2u(zeta, uk_dc, u);
  }

  A_k2a(Pressure);
  {
    Reset_phi(phi);
    if (SW_EQ == Shear_Navier_Stokes_Lees_Edwards) {
	  Make_phi_particle_OBL(phi, p);
	  if(SW_JANUS){
	    Make_phi_janus_particle_OBL(phi, work_v1, p); // +1/-1 janus polarity
	  }
    }else {
	  Make_phi_particle(phi, p);
	  if(SW_JANUS){
	    Make_phi_janus_particle(phi, work_v1, p); // +1/-1 janus polarity
	  }
    }
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
		int im=(i*NY*NZ_)+(j*NZ_)+k;
	  fprintf(fout,"%.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g\n"
		  ,u[0][im]
		  ,u[1][im]
		  ,u[2][im]
		  ,phi[im]
		  ,Pressure[im]
		  ,strain[1][im]
		  ,strain[4][im]
		  ,strain[2][im]
		  ,strain[0][im]
		  ,strain[3][im]
		  );
	}
      }
    }
  }
  fclose(fout);

}

void Output_udf(UDFManager *ufout
                , AVS_parameters &Avs_parameters
                , double **zeta
                , double *uk_dc
                , const Particle *p
                , const CTime &time
		)
{
  ufout->newRecord();
  ufout->put("E", 1.0);
  ufout->put("t", time.ts);
  for(int j = 0; j < Particle_Number; j++) {
    char str[256];
    sprintf(str, "Particles[%d]", j);
    Location target(str);
    ufout->put(target.sub("R.x"), p[j].x[0]);
    ufout->put(target.sub("R.y"), p[j].x[1]);
    ufout->put(target.sub("R.z"), p[j].x[2]);
    ufout->put(target.sub("R_raw.x"), p[j].x_nopbc[0]);
    ufout->put(target.sub("R_raw.y"), p[j].x_nopbc[1]);
    ufout->put(target.sub("R_raw.z"), p[j].x_nopbc[2]);
    ufout->put(target.sub("v.x"), p[j].v[0]);
    ufout->put(target.sub("v.y"), p[j].v[1]);
    ufout->put(target.sub("v.z"), p[j].v[2]);

    qtn_isnormal(p[j].q);
    ufout->put(target.sub("q.q0"), qtn_q0(p[j].q));
    ufout->put(target.sub("q.q1"), qtn_q1(p[j].q));
    ufout->put(target.sub("q.q2"), qtn_q2(p[j].q));
    ufout->put(target.sub("q.q3"), qtn_q3(p[j].q));
    ufout->put(target.sub("omega.x"), p[j].omega[0]);
    ufout->put(target.sub("omega.y"), p[j].omega[1]);
    ufout->put(target.sub("omega.z"), p[j].omega[2]);

    ufout->put(target.sub("f_hydro.x"), p[j].f_hydro_previous[0]);
    ufout->put(target.sub("f_hydro.y"), p[j].f_hydro_previous[1]);
    ufout->put(target.sub("f_hydro.z"), p[j].f_hydro_previous[2]);
    ufout->put(target.sub("torque_hydro.x"), p[j].torque_hydro_previous[0]);
    ufout->put(target.sub("torque_hydro.y"), p[j].torque_hydro_previous[1]);
    ufout->put(target.sub("torque_hydro.z"), p[j].torque_hydro_previous[2]);

    ufout->put(target.sub("f_r.x"), p[j].fr_previous[0]);
    ufout->put(target.sub("f_r.y"), p[j].fr_previous[1]);
    ufout->put(target.sub("f_r.z"), p[j].fr_previous[2]);
    ufout->put(target.sub("torque_r.x"), 0.0);
    ufout->put(target.sub("torque_r.y"), 0.0);
    ufout->put(target.sub("torque_r.z"), 0.0);


    ufout->put(target.sub("f_slip.x"), p[j].f_slip_previous[0]);
    ufout->put(target.sub("f_slip.y"), p[j].f_slip_previous[1]);
    ufout->put(target.sub("f_slip.z"), p[j].f_slip_previous[2]);
    ufout->put(target.sub("torque_slip.x"), p[j].torque_slip_previous[0]);
    ufout->put(target.sub("torque_slip.y"), p[j].torque_slip_previous[1]);
    ufout->put(target.sub("torque_slip.z"), p[j].torque_slip_previous[2]);
  }
  if(SW_PT == rigid){
    for(int rigidID = 0; rigidID < Rigid_Number; rigidID++){
      char str[256];
      sprintf(str, "RigidParticles[%d]", rigidID);
      Location target(str);

      int rigid_first_n = Rigid_Particle_Cumul[rigidID];
      quaternion qGs;
      qtn_init(qGs, p[rigid_first_n].q);

      ufout->put(target.sub("R.x"), xGs[rigidID][0]);
      ufout->put(target.sub("R.y"), xGs[rigidID][1]);
      ufout->put(target.sub("R.z"), xGs[rigidID][2]);
      ufout->put(target.sub("R_raw.x"), xGs_nopbc[rigidID][0]);
      ufout->put(target.sub("R_raw.y"), xGs_nopbc[rigidID][1]);
      ufout->put(target.sub("R_raw.z"), xGs_nopbc[rigidID][2]);
      ufout->put(target.sub("v.x"), velocityGs[rigidID][0]);
      ufout->put(target.sub("v.y"), velocityGs[rigidID][1]);
      ufout->put(target.sub("v.z"), velocityGs[rigidID][2]);

      ufout->put(target.sub("q.q0"), qtn_q0(qGs));
      ufout->put(target.sub("q.q1"), qtn_q1(qGs));
      ufout->put(target.sub("q.q2"), qtn_q2(qGs));
      ufout->put(target.sub("q.q3"), qtn_q3(qGs));
      ufout->put(target.sub("omega.x"), omegaGs[rigidID][0]);
      ufout->put(target.sub("omega.y"), omegaGs[rigidID][1]);
      ufout->put(target.sub("omega.z"), omegaGs[rigidID][2]);

      ufout->put(target.sub("f_hydro.x"), forceGs_previous[rigidID][0]);
      ufout->put(target.sub("f_hydro.y"), forceGs_previous[rigidID][1]);
      ufout->put(target.sub("f_hydro.z"), forceGs_previous[rigidID][2]);
      ufout->put(target.sub("torque_hydro.x"), torqueGs_previous[rigidID][0]);
      ufout->put(target.sub("torque_hydro.y"), torqueGs_previous[rigidID][1]);
      ufout->put(target.sub("torque_hydro.z"), torqueGs_previous[rigidID][2]);

      ufout->put(target.sub("f_r.x"), forceGrs_previous[rigidID][0]);
      ufout->put(target.sub("f_r.y"), forceGrs_previous[rigidID][1]);
      ufout->put(target.sub("f_r.z"), forceGrs_previous[rigidID][2]);
      ufout->put(target.sub("torque_r.x"), torqueGrs_previous[rigidID][0]);
      ufout->put(target.sub("torque_r.y"), torqueGrs_previous[rigidID][1]);
      ufout->put(target.sub("torque_r.z"), torqueGrs_previous[rigidID][2]);

      ufout->put(target.sub("f_slip.x"), 0.0);
      ufout->put(target.sub("f_slip.y"), 0.0);
      ufout->put(target.sub("f_slip.z"), 0.0);
      ufout->put(target.sub("torque_slip.x"), 0.0);
      ufout->put(target.sub("torque_slip.y"), 0.0);
      ufout->put(target.sub("torque_slip.z"), 0.0);
    }
  }
}

void Output_avs_charge(AVS_parameters &Avs_parameters
		       ,double **zeta
		       ,double *uk_dc
		       ,double **Concentration
		       ,Particle *p
		       ,const CTime &time
		       ){
  Zeta_k2u(zeta, uk_dc, u);
  
  double *potential = f_particle[0];
  double *dmy_value0 = f_particle[1];
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
	  float dmy = up[0][(i*NY*NZ_)+(j*NZ_)+k] * dmy_surface_area;
	  fwrite(&dmy,sizeof(float),1,fout);
	}
      }
    }
    for(int k=Avs_parameters.kstart; k<= Avs_parameters.kend; k++){
      for(int j=Avs_parameters.jstart; j<= Avs_parameters.jend; j++){
	for(int i=Avs_parameters.istart; i<=Avs_parameters.iend; i++){
	  float dmy=0.;
	  for(int n=0; n<N_spec; n++){
	    dmy += (float)(Elementary_charge*Valency[n]*Concentration[n][(i*NY*NZ_)+(j*NZ_)+k]);
	  }
	  dmy *= (1.-phi[(i*NY*NZ_)+(j*NZ_)+k]);
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
	    dmy += Elementary_charge*Valency[n]*Concentration[n][(i*NY*NZ_)+(j*NZ_)+k];
	  }
	  fprintf(fout,"%.3g %.3g %.3g %.3g %.3g %.3g %.3g\n"
		  ,u[0][(i*NY*NZ_)+(j*NZ_)+k],u[1][(i*NY*NZ_)+(j*NZ_)+k],u[2][(i*NY*NZ_)+(j*NZ_)+k]
		  ,phi[(i*NY*NZ_)+(j*NZ_)+k]
		  ,up[0][(i*NY*NZ_)+(j*NZ_)+k]*dmy_surface_area
		  ,(1.-phi[(i*NY*NZ_)+(j*NZ_)+k])*dmy
		  ,potential[(i*NY*NZ_)+(j*NZ_)+k]
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

