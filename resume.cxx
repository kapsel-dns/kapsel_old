/*!
  \file resume.cxx
  \brief Routines to read/write restart file
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
*/
#include "resume.h"

void Save_Restart_udf(
		      double **zeta
		      ,double *uk_dc
		      ,const Particle *p
		      ,const CTime &time
		      ,double **conc_k
		      ){
  ufres->put("resume.Calculation","CONTINUE");
  {
    {
      for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	  for(int k=0; k<NZ_; k++){
	    char str[256];
	    {
	      sprintf(str,"resume.CONTINUE.Saved_Data.Zeta[%d][%d][%d]",i,j,k);
	      Location target(str);
	      ufres->put(target.sub("zeta0"),zeta[0][(i*NY*NZ_)+(j*NZ_)+k]);
	      ufres->put(target.sub("zeta1"),zeta[1][(i*NY*NZ_)+(j*NZ_)+k]);
	    }
	    if(Electrolyte){
	      for(int n=0;n<N_spec;n++){ // Two_fluid では N_spec =1
	      sprintf(str,"resume.CONTINUE.Saved_Data.Concentration[%d][%d][%d][%d]",n,i,j,k);
	      Location target(str);
	      ufres->put(target.sub("ck"),conc_k[n][(i*NY*NZ_)+(j*NZ_)+k]);
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

	qtn_isnormal(p[j].q);
	ufres->put(target.sub("q.q0"), qtn_q0(p[j].q));
	ufres->put(target.sub("q.q1"), qtn_q1(p[j].q));
	ufres->put(target.sub("q.q2"), qtn_q2(p[j].q));
	ufres->put(target.sub("q.q3"), qtn_q3(p[j].q));
	qtn_isnormal(p[j].q_old);
	ufres->put(target.sub("q_old.q0"), qtn_q0(p[j].q_old));
	ufres->put(target.sub("q_old.q1"), qtn_q1(p[j].q_old));
	ufres->put(target.sub("q_old.q2"), qtn_q2(p[j].q_old));
	ufres->put(target.sub("q_old.q3"), qtn_q3(p[j].q_old));

	ufres->put(target.sub("f_hydro.x"),p[j].f_hydro[0]);
	ufres->put(target.sub("f_hydro.y"),p[j].f_hydro[1]);
	ufres->put(target.sub("f_hydro.z"),p[j].f_hydro[2]);
	ufres->put(target.sub("f_hydro_previous.x"),p[j].f_hydro_previous[0]);
	ufres->put(target.sub("f_hydro_previous.y"),p[j].f_hydro_previous[1]);
	ufres->put(target.sub("f_hydro_previous.z"),p[j].f_hydro_previous[2]);
	ufres->put(target.sub("f_hydro1.x"),p[j].f_hydro1[0]);
	ufres->put(target.sub("f_hydro1.y"),p[j].f_hydro1[1]);
	ufres->put(target.sub("f_hydro1.z"),p[j].f_hydro1[2]);
	ufres->put(target.sub("f_slip.x"), p[j].f_slip[0]);
	ufres->put(target.sub("f_slip.y"), p[j].f_slip[1]);
	ufres->put(target.sub("f_slip.z"), p[j].f_slip[2]);
	ufres->put(target.sub("f_slip_previous.x"), p[j].f_slip_previous[0]);
	ufres->put(target.sub("f_slip_previous.y"), p[j].f_slip_previous[1]);
	ufres->put(target.sub("f_slip_previous.z"), p[j].f_slip_previous[2]);
	

	ufres->put(target.sub("fr.x"),p[j].fr[0]);
	ufres->put(target.sub("fr.y"),p[j].fr[1]);
	ufres->put(target.sub("fr.z"),p[j].fr[2]);
	ufres->put(target.sub("fr_previous.x"),p[j].fr_previous[0]);
	ufres->put(target.sub("fr_previous.y"),p[j].fr_previous[1]);
	ufres->put(target.sub("fr_previous.z"),p[j].fr_previous[2]);

	ufres->put(target.sub("omega.x"),p[j].omega[0]);
	ufres->put(target.sub("omega.y"),p[j].omega[1]);
	ufres->put(target.sub("omega.z"),p[j].omega[2]);
	ufres->put(target.sub("omega_old.x"),p[j].omega_old[0]);
	ufres->put(target.sub("omega_old.y"),p[j].omega_old[1]);
	ufres->put(target.sub("omega_old.z"),p[j].omega_old[2]);

	ufres->put(target.sub("torque_hydro.x"),p[j].torque_hydro[0]);
	ufres->put(target.sub("torque_hydro.y"),p[j].torque_hydro[1]);
	ufres->put(target.sub("torque_hydro.z"),p[j].torque_hydro[2]);
	ufres->put(target.sub("torque_hydro_previous.x"),p[j].torque_hydro_previous[0]);
	ufres->put(target.sub("torque_hydro_previous.y"),p[j].torque_hydro_previous[1]);
	ufres->put(target.sub("torque_hydro_previous.z"),p[j].torque_hydro_previous[2]);
	ufres->put(target.sub("torque_hydro1.x"),p[j].torque_hydro1[0]);
	ufres->put(target.sub("torque_hydro1.y"),p[j].torque_hydro1[1]);
	ufres->put(target.sub("torque_hydro1.z"),p[j].torque_hydro1[2]);
	ufres->put(target.sub("torque_slip.x"), p[j].torque_slip[0]);
	ufres->put(target.sub("torque_slip.y"), p[j].torque_slip[1]);
	ufres->put(target.sub("torque_slip.z"), p[j].torque_slip[2]);
	ufres->put(target.sub("torque_slip_previous.x"), p[j].torque_slip_previous[0]);
	ufres->put(target.sub("torque_slip_previous.y"), p[j].torque_slip_previous[1]);
	ufres->put(target.sub("torque_slip_previous.z"), p[j].torque_slip_previous[2]);

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
      {
	  char str[256];
	  sprintf(str,"resume.CONTINUE.Saved_Data.oblique");
	  Location target(str);
	  ufres->put(target.sub("degree_oblique"),degree_oblique);
      }
    }
  }
}
void Force_restore_parameters(double **zeta
			      ,double *uk_dc
			      ,Particle *p
			      ,CTime &time
			      ,double **conc_k
    ){
    {
	{
	    for(int i=0; i<NX; i++){
		for(int j=0; j<NY; j++){
		    for(int k=0; k<NZ_; k++){
			char str[256];
			{
			    sprintf(str,"resume.CONTINUE.Saved_Data.Zeta[%d][%d][%d]",i,j,k);
			    Location target(str);
			    ufin->get(target.sub("zeta0"),zeta[0][(i*NY*NZ_)+(j*NZ_)+k]);
			    ufin->get(target.sub("zeta1"),zeta[1][(i*NY*NZ_)+(j*NZ_)+k]);
			}
			if(Electrolyte ){
			    for(int n=0;n<N_spec;n++){ // Two_fluid では N_spec =1
				sprintf(str,"resume.CONTINUE.Saved_Data.Concentration[%d][%d][%d][%d]",n,i,j,k);
				Location target(str);
				ufin->get(target.sub("ck"),conc_k[n][(i*NY*NZ_)+(j*NZ_)+k]);
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
		for(int d = 0; d < DIM; d++){
		  p[j].x_nopbc[d] = p[j].x[d];
		}

		ufin->get(target.sub("v.x"),p[j].v[0]);
		ufin->get(target.sub("v.y"),p[j].v[1]);
		ufin->get(target.sub("v.z"),p[j].v[2]);
		
		ufin->get(target.sub("v_old.x"),p[j].v_old[0]);
		ufin->get(target.sub("v_old.y"),p[j].v_old[1]);
		ufin->get(target.sub("v_old.z"),p[j].v_old[2]);
		
		{
		  double q0, q1, q2, q3;
		  ufin->get(target.sub("q.q0"), q0);
		  ufin->get(target.sub("q.q1"), q1);
		  ufin->get(target.sub("q.q2"), q2);
		  ufin->get(target.sub("q.q3"), q3);
		  qtn_init(p[j].q, q0, q1, q2, q3);
		  qtn_isnormal(p[j].q);

		  ufin->get(target.sub("q_old.q0"), q0);
		  ufin->get(target.sub("q_old.q1"), q1);
		  ufin->get(target.sub("q_old.q2"), q2);
		  ufin->get(target.sub("q_old.q3"), q3);
		  qtn_init(p[j].q_old, q0, q1, q2, q3);
		  qtn_isnormal(p[j].q_old);
		}
		
		ufin->get(target.sub("f_hydro.x"),p[j].f_hydro[0]);
		ufin->get(target.sub("f_hydro.y"),p[j].f_hydro[1]);
		ufin->get(target.sub("f_hydro.z"),p[j].f_hydro[2]);
		ufin->get(target.sub("f_hydro_previous.x"),p[j].f_hydro_previous[0]);
		ufin->get(target.sub("f_hydro_previous.y"),p[j].f_hydro_previous[1]);
		ufin->get(target.sub("f_hydro_previous.z"),p[j].f_hydro_previous[2]);
		ufin->get(target.sub("f_hydro1.x"),p[j].f_hydro1[0]);
		ufin->get(target.sub("f_hydro1.y"),p[j].f_hydro1[1]);
		ufin->get(target.sub("f_hydro1.z"),p[j].f_hydro1[2]);
		ufin->get(target.sub("f_slip.x"), p[j].f_slip[0]);
		ufin->get(target.sub("f_slip.y"), p[j].f_slip[1]);
		ufin->get(target.sub("f_slip.z"), p[j].f_slip[2]);
		ufin->get(target.sub("f_slip_previous.x"), p[j].f_slip_previous[0]);
		ufin->get(target.sub("f_slip_previous.y"), p[j].f_slip_previous[1]);
		ufin->get(target.sub("f_slip_previous.z"), p[j].f_slip_previous[2]);
		
		ufin->get(target.sub("fr.x"),p[j].fr[0]);
		ufin->get(target.sub("fr.y"),p[j].fr[1]);
		ufin->get(target.sub("fr.z"),p[j].fr[2]);
		ufin->get(target.sub("fr_previous.x"),p[j].fr_previous[0]);
		ufin->get(target.sub("fr_previous.y"),p[j].fr_previous[1]);
		ufin->get(target.sub("fr_previous.z"),p[j].fr_previous[2]);
		
		ufin->get(target.sub("omega.x"),p[j].omega[0]);
		ufin->get(target.sub("omega.y"),p[j].omega[1]);
		ufin->get(target.sub("omega.z"),p[j].omega[2]);
		
		ufin->get(target.sub("omega_old.x"),p[j].omega_old[0]);
		ufin->get(target.sub("omega_old.y"),p[j].omega_old[1]);
		ufin->get(target.sub("omega_old.z"),p[j].omega_old[2]);
		
		ufin->get(target.sub("torque_hydro.x"),p[j].torque_hydro[0]);
		ufin->get(target.sub("torque_hydro.y"),p[j].torque_hydro[1]);
		ufin->get(target.sub("torque_hydro.z"),p[j].torque_hydro[2]);
		ufin->get(target.sub("torque_hydro_previous.x"),p[j].torque_hydro_previous[0]);
		ufin->get(target.sub("torque_hydro_previous.y"),p[j].torque_hydro_previous[1]);
		ufin->get(target.sub("torque_hydro_previous.z"),p[j].torque_hydro_previous[2]);
		ufin->get(target.sub("torque_hydro1.x"),p[j].torque_hydro1[0]);
		ufin->get(target.sub("torque_hydro1.y"),p[j].torque_hydro1[1]);
		ufin->get(target.sub("torque_hydro1.z"),p[j].torque_hydro1[2]);
		ufin->get(target.sub("torque_slip.x"), p[j].torque_slip[0]);
		ufin->get(target.sub("torque_slip.y"), p[j].torque_slip[1]);
		ufin->get(target.sub("torque_slip.z"), p[j].torque_slip[2]);
		ufin->get(target.sub("torque_slip_previous.x"), p[j].torque_slip_previous[0]);
		ufin->get(target.sub("torque_slip_previous.y"), p[j].torque_slip_previous[1]);
		ufin->get(target.sub("torque_slip_previous.z"), p[j].torque_slip_previous[2]);
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
	    {
		char str[256];
		sprintf(str,"resume.CONTINUE.Saved_Data.oblique");
		Location target(str);
		ufin->get(target.sub("degree_oblique"),degree_oblique);
	    }      
	}
    }
}
