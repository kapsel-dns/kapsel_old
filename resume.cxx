//
// $Id: resume.cxx,v 1.1 2006/05/15 09:19:02 kin Exp $
//
#include "resume.h"

void Save_Restart_udf(UDFManager *ufout
		      ,Value *zeta
		      ,double *uk_dc
		      ,const Particle *p
		      ,const CTime &time
		      ,Value *conc_k
		      ){
  {
    {
      char sw_eq[128];
      ufres->put("E",1.0);
      ufres->put("t",time.ts);
      for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	  for(int k=0; k<NZ_; k++){
	    char str[256];
	    {
	      sprintf(str,"resume.CONTINUE.Saved_Data.Zeta[%d][%d][%d]",i,j,k);
	      Location target(str);
	      ufres->put(target.sub("zeta0"),zeta[0][i][j][k]);
	      ufres->put(target.sub("zeta1"),zeta[1][i][j][k]);
	    }
	    if( SW_EQ == Two_fluid || Electrolyte ){
	      for(int n=0;n<N_spec;n++){ // Two_fluid では N_spec =1
	      sprintf(str,"resume.CONTINUE.Saved_Data.Concentration[%d][%d][%d][%d]",n,i,j,k);
	      Location target(str);
	      ufres->put(target.sub("ck"),conc_k[n][i][j][k]);
	      }
	    }
	  }
	}
      }
     
      for(int j=0;j<Particle_Number;j++){
	char str[256];
	sprintf(str,"resume.CONTINUE.Saved_Data.Particles[%d]",j);
	Location target(str);
	ufres->put(target.sub("R.x"),p[j].x[0]);
	ufres->put(target.sub("R.y"),p[j].x[1]);
	ufres->put(target.sub("R.z"),p[j].x[2]);

	ufres->put(target.sub("v.x"),p[j].v[0]);
	ufres->put(target.sub("v.y"),p[j].v[1]);
	ufres->put(target.sub("v.z"),p[j].v[2]);
	ufres->put(target.sub("v_old.x"),p[j].v_old[0]);
	ufres->put(target.sub("v_old.y"),p[j].v_old[1]);
	ufres->put(target.sub("v_old.z"),p[j].v_old[2]);

	ufres->put(target.sub("f_hydro.x"),p[j].f_hydro[0]);
	ufres->put(target.sub("f_hydro.y"),p[j].f_hydro[1]);
	ufres->put(target.sub("f_hydro.z"),p[j].f_hydro[2]);
	ufres->put(target.sub("f_hydro1.x"),p[j].f_hydro1[0]);
	ufres->put(target.sub("f_hydro1.y"),p[j].f_hydro1[1]);
	ufres->put(target.sub("f_hydro1.z"),p[j].f_hydro1[2]);

	ufres->put(target.sub("fr.x"),p[j].fr[0]);
	ufres->put(target.sub("fr.y"),p[j].fr[1]);
	ufres->put(target.sub("fr.z"),p[j].fr[2]);
	ufres->put(target.sub("fr_previous.x"),p[j].fr_previous[0]);
	ufres->put(target.sub("fr_previous.y"),p[j].fr_previous[1]);
	ufres->put(target.sub("fr_previous.z"),p[j].fr_previous[2]);

	ufres->put(target.sub("fv.x"),p[j].fv[0]);
	ufres->put(target.sub("fv.y"),p[j].fv[1]);
	ufres->put(target.sub("fv.z"),p[j].fv[2]);
	ufres->put(target.sub("fv_previous.x"),p[j].fv_previous[0]);
	ufres->put(target.sub("fv_previous.y"),p[j].fv_previous[1]);
	ufres->put(target.sub("fv_previous.z"),p[j].fv_previous[2]);

	ufres->put(target.sub("f_collison.x"),p[j].f_collison[0]);
	ufres->put(target.sub("f_collison.y"),p[j].f_collison[1]);
	ufres->put(target.sub("f_collison.z"),p[j].f_collison[2]);
	ufres->put(target.sub("f_collison_previous.x"),p[j].f_collison_previous[0]);
	ufres->put(target.sub("f_collison_previous.y"),p[j].f_collison_previous[1]);
	ufres->put(target.sub("f_collison_previous.z"),p[j].f_collison_previous[2]);

	ufres->put(target.sub("omega.x"),p[j].omega[0]);
	ufres->put(target.sub("omega.y"),p[j].omega[1]);
	ufres->put(target.sub("omega.z"),p[j].omega[2]);
	ufres->put(target.sub("omega_old.x"),p[j].omega_old[0]);
	ufres->put(target.sub("omega_old.y"),p[j].omega_old[1]);
	ufres->put(target.sub("omega_old.z"),p[j].omega_old[2]);

	ufres->put(target.sub("torque_hydro.x"),p[j].torque_hydro[0]);
	ufres->put(target.sub("torque_hydro.y"),p[j].torque_hydro[1]);
	ufres->put(target.sub("torque_hydro.z"),p[j].torque_hydro[2]);
	ufres->put(target.sub("torque_hydro1.x"),p[j].torque_hydro1[0]);
	ufres->put(target.sub("torque_hydro1.y"),p[j].torque_hydro1[1]);
	ufres->put(target.sub("torque_hydro1.z"),p[j].torque_hydro1[2]);

	ufres->put(target.sub("torquer.x"),p[j].torquer[0]);
	ufres->put(target.sub("torquer.y"),p[j].torquer[1]);
	ufres->put(target.sub("torquer.z"),p[j].torquer[2]);
	ufres->put(target.sub("torquer_previous.x"),p[j].torquer_previous[0]);
	ufres->put(target.sub("torquer_previous.y"),p[j].torquer_previous[1]);
	ufres->put(target.sub("torquer_previous.z"),p[j].torquer_previous[2]);

	ufres->put(target.sub("torquev.x"),p[j].torquev[0]);
	ufres->put(target.sub("torquev.y"),p[j].torquev[1]);
	ufres->put(target.sub("torquev.z"),p[j].torquev[2]);
	ufres->put(target.sub("torquev_previous.x"),p[j].torquev_previous[0]);
	ufres->put(target.sub("torquev_previous.y"),p[j].torquev_previous[1]);
	ufres->put(target.sub("torquev_previous.z"),p[j].torquev_previous[2]);

      }
      {
	char str[256];
	sprintf(str,"resume.CONTINUE.Saved_Data.uk_dc");
	Location target(str);
	ufres->put(target.sub("x"),uk_dc[0]);
	ufres->put(target.sub("y"),uk_dc[1]);
	ufres->put(target.sub("z"),uk_dc[2]);
      }
      {
	char str[256];
	sprintf(str,"resume.CONTINUE.Saved_Data.jikan");
	Location target(str);
	ufres->put(target.sub("ts"),time.ts);
	ufres->put(target.sub("time"),time.time);
      }
      
    }
  }
}

void Force_restore_parameters(UDFManager *ufout
			      ,Value *zeta
			      ,double *uk_dc
			      ,Particle *p
			      ,CTime &time
			      ,Value *conc_k
			      ){
  {
    {
      char sw_eq[128];
      for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	  for(int k=0; k<NZ_; k++){
	    char str[256];
	    {
	      sprintf(str,"resume.CONTINUE.Saved_Data.Zeta[%d][%d][%d]",i,j,k);
	      Location target(str);
	      ufin->get(target.sub("zeta0"),zeta[0][i][j][k]);
	      ufin->get(target.sub("zeta1"),zeta[1][i][j][k]);
	    }
	    if( SW_EQ == Two_fluid || Electrolyte ){
	      for(int n=0;n<N_spec;n++){ // Two_fluid では N_spec =1
	      sprintf(str,"resume.CONTINUE.Saved_Data.Concentration[%d][%d][%d][%d]",n,i,j,k);
	      Location target(str);
	      ufin->get(target.sub("ck"),conc_k[n][i][j][k]);
	      }
	    }
	  }
	}
      }


      for(int j=0;j<Particle_Number;j++){
	char str[256];
	sprintf(str,"resume.CONTINUE.Saved_Data.Particles[%d]",j);
	Location target(str);
	ufin->get(target.sub("R.x"),p[j].x[0]);
	ufin->get(target.sub("R.y"),p[j].x[1]);
	ufin->get(target.sub("R.z"),p[j].x[2]);
	ufin->get(target.sub("v.x"),p[j].v[0]);
	ufin->get(target.sub("v.y"),p[j].v[1]);
	ufin->get(target.sub("v.z"),p[j].v[2]);
	
	ufin->get(target.sub("v_old.x"),p[j].v_old[0]);
	ufin->get(target.sub("v_old.y"),p[j].v_old[1]);
	ufin->get(target.sub("v_old.z"),p[j].v_old[2]);

	ufin->get(target.sub("f_hydro.x"),p[j].f_hydro[0]);
	ufin->get(target.sub("f_hydro.y"),p[j].f_hydro[1]);
	ufin->get(target.sub("f_hydro.z"),p[j].f_hydro[2]);
	ufin->get(target.sub("f_hydro1.x"),p[j].f_hydro1[0]);
	ufin->get(target.sub("f_hydro1.y"),p[j].f_hydro1[1]);
	ufin->get(target.sub("f_hydro1.z"),p[j].f_hydro1[2]);

	ufin->get(target.sub("fr.x"),p[j].fr[0]);
	ufin->get(target.sub("fr.y"),p[j].fr[1]);
	ufin->get(target.sub("fr.z"),p[j].fr[2]);
	ufin->get(target.sub("fr_previous.x"),p[j].fr_previous[0]);
	ufin->get(target.sub("fr_previous.y"),p[j].fr_previous[1]);
	ufin->get(target.sub("fr_previous.z"),p[j].fr_previous[2]);

	ufin->get(target.sub("fv.x"),p[j].fv[0]);
	ufin->get(target.sub("fv.y"),p[j].fv[1]);
	ufin->get(target.sub("fv.z"),p[j].fv[2]);
	ufin->get(target.sub("fv_previous.x"),p[j].fv_previous[0]);
	ufin->get(target.sub("fv_previous.y"),p[j].fv_previous[1]);
	ufin->get(target.sub("fv_previous.z"),p[j].fv_previous[2]);

	ufin->get(target.sub("f_collison.x"),p[j].f_collison[0]);
	ufin->get(target.sub("f_collison.y"),p[j].f_collison[1]);
	ufin->get(target.sub("f_collison.z"),p[j].f_collison[2]);
	ufin->get(target.sub("f_collison_previous.x"),p[j].f_collison_previous[0]);
	ufin->get(target.sub("f_collison_previous.y"),p[j].f_collison_previous[1]);
	ufin->get(target.sub("f_collison_previous.z"),p[j].f_collison_previous[2]);
	
	ufin->get(target.sub("omega.x"),p[j].omega[0]);
	ufin->get(target.sub("omega.y"),p[j].omega[1]);
	ufin->get(target.sub("omega.z"),p[j].omega[2]);
	
	ufin->get(target.sub("omega_old.x"),p[j].omega_old[0]);
	ufin->get(target.sub("omega_old.y"),p[j].omega_old[1]);
	ufin->get(target.sub("omega_old.z"),p[j].omega_old[2]);

	ufin->get(target.sub("torque_hydro.x"),p[j].torque_hydro[0]);
	ufin->get(target.sub("torque_hydro.y"),p[j].torque_hydro[1]);
	ufin->get(target.sub("torque_hydro.z"),p[j].torque_hydro[2]);
	ufin->get(target.sub("torque_hydro1.x"),p[j].torque_hydro1[0]);
	ufin->get(target.sub("torque_hydro1.y"),p[j].torque_hydro1[1]);
	ufin->get(target.sub("torque_hydro1.z"),p[j].torque_hydro1[2]);

	ufin->get(target.sub("torquer.x"),p[j].torquer[0]);
	ufin->get(target.sub("torquer.y"),p[j].torquer[1]);
	ufin->get(target.sub("torquer.z"),p[j].torquer[2]);
	ufin->get(target.sub("torquer_previous.x"),p[j].torquer_previous[0]);
	ufin->get(target.sub("torquer_previous.y"),p[j].torquer_previous[1]);
	ufin->get(target.sub("torquer_previous.z"),p[j].torquer_previous[2]);

	ufin->get(target.sub("torquev.x"),p[j].torquev[0]);
	ufin->get(target.sub("torquev.y"),p[j].torquev[1]);
	ufin->get(target.sub("torquev.z"),p[j].torquev[2]);
	ufin->get(target.sub("torquev_previous.x"),p[j].torquev_previous[0]);
	ufin->get(target.sub("torquev_previous.y"),p[j].torquev_previous[1]);
	ufin->get(target.sub("torquev_previous.z"),p[j].torquev_previous[2]);
	
      }
      {
	char str[256];
	sprintf(str,"resume.CONTINUE.Saved_Data.uk_dc");
	Location target(str);
	
	ufin->get(target.sub("x"),uk_dc[0]);
	ufin->get(target.sub("y"),uk_dc[1]);
	ufin->get(target.sub("z"),uk_dc[2]);
	
      }
      {
	char str[256];
	sprintf(str,"resume.CONTINUE.Saved_Data.jikan");
	Location target(str);
	ufin->get(target.sub("ts"),time.ts);
	ufin->get(target.sub("time"),time.time);
      }      
    }
  }
}
