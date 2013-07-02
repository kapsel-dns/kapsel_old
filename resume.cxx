/*!
  \file resume.cxx
  \brief Routines to read/write restart file
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
*/
#include "resume.h"

void Save_Restart_udf(double **zeta,
                      double *uk_dc,
                      const Particle *p,
                      const CTime &time,
                      double **conc_k){
  ufres->put("resume.Calculation", "CONTINUE");
  {//fluid data
    int im;
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
        for(int k=0; k<NZ_; k++){
          im = (i*NY*NZ_) + (j*NZ_) + k;
          char str[256];
          {
            sprintf(str, "resume.CONTINUE.Saved_Data.Zeta[%d][%d][%d]",i,j,k);
            Location target(str);
            ufres->put(target.sub("zeta0"),zeta[0][im]);
            ufres->put(target.sub("zeta1"),zeta[1][im]);
          }

          if(Electrolyte){
            for(int n=0;n<N_spec;n++){
              sprintf(str, "resume.CONTINUE.Saved_Data.Concentration[%d][%d][%d][%d]",n,i,j,k);
              Location target(str);
              ufres->put(target.sub("ck"),conc_k[n][im]);
            }
          }//Electrolyte

        }//k
      }//j
    }//i
  }
  {//uk_dc
    char str[256];
    sprintf(str,"resume.CONTINUE.Saved_Data.uk_dc");
    Location target(str);
    ufres->put(target.sub("x"), uk_dc[0]);
    ufres->put(target.sub("y"), uk_dc[1]);
    ufres->put(target.sub("z"), uk_dc[2]);
  }

  //particle data
  if(SW_PT != rigid){
    Save_Particle_udf(p, Particle_Number);
  }else{
    Particle *rigid_p = new Particle [Rigid_Number];
    Get_Rigid_Particle_Data(rigid_p, p);

    Save_Particle_udf(rigid_p, Rigid_Number);
    Save_Rigid_Particle_udf();

    delete [] rigid_p;
  }
  
  {//time
    char str[256];
    sprintf(str,"resume.CONTINUE.Saved_Data.jikan");
    Location target(str);
    ufres->put(target.sub("ts"),time.ts);
    ufres->put(target.sub("time"),time.time);
  }

  {//strain
    char str[256];
    sprintf(str,"resume.CONTINUE.Saved_Data.oblique");
    Location target(str);
    ufres->put(target.sub("degree_oblique"),degree_oblique);
  }
}

void Get_Rigid_Particle_Data( Particle *rigid_p, const Particle *p){
#pragma omp parallel for schedule(dynamic, 1)
  for(int rigidID = 0; rigidID < Rigid_Number; rigidID++){
    int rigid_first_n = Rigid_Particle_Cumul[rigidID];

    //orientations
    qtn_init(rigid_p[rigidID].q, p[rigid_first_n].q);
    qtn_init(rigid_p[rigidID].q_old, p[rigid_first_n].q_old);
    for(int d = 0; d < DIM; d++){
      for(int l = 0; l < DIM; l++){
        rigid_p[rigidID].QR[d][l] = 0.0;
        rigid_p[rigidID].QR_old[d][l] = 0.0;
      }
    }

    //vector data
    for(int d = 0; d < DIM; d++){
      //positions
      rigid_p[rigidID].x[d] = xGs[rigidID][d];
      rigid_p[rigidID].x_previous[d] = xGs_previous[rigidID][d];
      rigid_p[rigidID].x_nopbc[d] = xGs_nopbc[rigidID][d];

      //velocities
      rigid_p[rigidID].v[d] = velocityGs[rigidID][d];
      rigid_p[rigidID].v_old[d] = velocityGs_old[rigidID][d];
      rigid_p[rigidID].v_slip[d] = 0.0;

      //angular velocities
      rigid_p[rigidID].omega[d] = omegaGs[rigidID][d];
      rigid_p[rigidID].omega_old[d] = omegaGs_old[rigidID][d];
      rigid_p[rigidID].omega_slip[d] = 0.0;
      
      //hydrodynamic forces
      rigid_p[rigidID].f_hydro[d] = forceGs[rigidID][d];
      rigid_p[rigidID].f_hydro_previous[d] = forceGs_previous[rigidID][d];
      rigid_p[rigidID].f_hydro1[d] = 0.0;
      rigid_p[rigidID].f_slip[d] = 0.0;
      rigid_p[rigidID].f_slip_previous[d] = 0.0;

      //ext+lj+random forces
      rigid_p[rigidID].fr[d] = forceGrs[rigidID][d];
      rigid_p[rigidID].fr_previous[d] = forceGrs_previous[rigidID][d];

      //hydrodynamic torques
      rigid_p[rigidID].torque_hydro[d] = torqueGs[rigidID][d];
      rigid_p[rigidID].torque_hydro_previous[d] = torqueGs_previous[rigidID][d];
      rigid_p[rigidID].torque_hydro1[d] = 0.0;
      rigid_p[rigidID].torque_slip[d] = 0.0;
      rigid_p[rigidID].torque_slip_previous[d] = 0.0;

      //momentum_depend_fr
      rigid_p[rigidID].momentum_depend_fr[d] = 0.0;
    }

    //temp data (ignored)
    rigid_p[rigidID].mass = 0.0;
    rigid_p[rigidID].surface_mass = 0.0;
    for(int d = 0; d < DIM; d++){
      rigid_p[rigidID].mass_center[d] = 0.0;
      rigid_p[rigidID].surface_mass_center[d] = 0.0;
      rigid_p[rigidID].surface_dv[d] = 0.0;
      rigid_p[rigidID].surface_dw[d] = 0.0;
    }
    for(int d = 0; d < DIM; d++){
      for(int l = 0; l < DIM; l++){
        rigid_p[rigidID].inertia[d][l] = 0.0;
        rigid_p[rigidID].surface_inertia[d][l] = 0.0;
      }
    }
  }//rigid
}

void Save_Rigid_Particle_udf(){
  if(SW_PT == rigid){
    for(int j = 0; j < Particle_Number; j++){
      char str[256];
      sprintf(str, "resume.CONTINUE.Saved_Data.GR_body[%d]",j);
      Location target(str);

      ufres->put(target.sub("x"), GRvecs_body[j][0]);
      ufres->put(target.sub("y"), GRvecs_body[j][1]);
      ufres->put(target.sub("z"), GRvecs_body[j][2]);
    }
  }
}

void Save_Particle_udf(const Particle *p, const int &n_out_particles){
  
  for(int j=0;j<n_out_particles;j++){
    char str[256];
    sprintf(str,"resume.CONTINUE.Saved_Data.Particles[%d]",j);
    Location target(str);

    //positions
    ufres->put(target.sub("R.x"),p[j].x[0]);
    ufres->put(target.sub("R.y"),p[j].x[1]);
    ufres->put(target.sub("R.z"),p[j].x[2]);
    
    //velocities
    ufres->put(target.sub("v.x"),p[j].v[0]);
    ufres->put(target.sub("v.y"),p[j].v[1]);
    ufres->put(target.sub("v.z"),p[j].v[2]);

    //old velocities
    ufres->put(target.sub("v_old.x"),p[j].v_old[0]);
    ufres->put(target.sub("v_old.y"),p[j].v_old[1]);
    ufres->put(target.sub("v_old.z"),p[j].v_old[2]);
    
    //orientations
    qtn_isnormal(p[j].q);
    ufres->put(target.sub("q.q0"), qtn_q0(p[j].q));
    ufres->put(target.sub("q.q1"), qtn_q1(p[j].q));
    ufres->put(target.sub("q.q2"), qtn_q2(p[j].q));
    ufres->put(target.sub("q.q3"), qtn_q3(p[j].q));

    //old orientations
    qtn_isnormal(p[j].q_old);
    ufres->put(target.sub("q_old.q0"), qtn_q0(p[j].q_old));
    ufres->put(target.sub("q_old.q1"), qtn_q1(p[j].q_old));
    ufres->put(target.sub("q_old.q2"), qtn_q2(p[j].q_old));
    ufres->put(target.sub("q_old.q3"), qtn_q3(p[j].q_old));

    //hydrodynamic force
    ufres->put(target.sub("f_hydro.x"),p[j].f_hydro[0]);
    ufres->put(target.sub("f_hydro.y"),p[j].f_hydro[1]);
    ufres->put(target.sub("f_hydro.z"),p[j].f_hydro[2]);

    //old hydrodynamic force
    ufres->put(target.sub("f_hydro_previous.x"),p[j].f_hydro_previous[0]);
    ufres->put(target.sub("f_hydro_previous.y"),p[j].f_hydro_previous[1]);
    ufres->put(target.sub("f_hydro_previous.z"),p[j].f_hydro_previous[2]);

    //electrolyte force
    ufres->put(target.sub("f_hydro1.x"),p[j].f_hydro1[0]);
    ufres->put(target.sub("f_hydro1.y"),p[j].f_hydro1[1]);
    ufres->put(target.sub("f_hydro1.z"),p[j].f_hydro1[2]);

    //slip force
    ufres->put(target.sub("f_slip.x"), p[j].f_slip[0]);
    ufres->put(target.sub("f_slip.y"), p[j].f_slip[1]);
    ufres->put(target.sub("f_slip.z"), p[j].f_slip[2]);

    //old slip force
    ufres->put(target.sub("f_slip_previous.x"), p[j].f_slip_previous[0]);
    ufres->put(target.sub("f_slip_previous.y"), p[j].f_slip_previous[1]);
    ufres->put(target.sub("f_slip_previous.z"), p[j].f_slip_previous[2]);
    
    //random+lj+ext force
    ufres->put(target.sub("fr.x"),p[j].fr[0]);
    ufres->put(target.sub("fr.y"),p[j].fr[1]);
    ufres->put(target.sub("fr.z"),p[j].fr[2]);

    //old random+lj+ext force
    ufres->put(target.sub("fr_previous.x"),p[j].fr_previous[0]);
    ufres->put(target.sub("fr_previous.y"),p[j].fr_previous[1]);
    ufres->put(target.sub("fr_previous.z"),p[j].fr_previous[2]);

    //angular velocity
    ufres->put(target.sub("omega.x"),p[j].omega[0]);
    ufres->put(target.sub("omega.y"),p[j].omega[1]);
    ufres->put(target.sub("omega.z"),p[j].omega[2]);

    //old angular velocity
    ufres->put(target.sub("omega_old.x"),p[j].omega_old[0]);
    ufres->put(target.sub("omega_old.y"),p[j].omega_old[1]);
    ufres->put(target.sub("omega_old.z"),p[j].omega_old[2]);

    //hydrodynamic torque
    ufres->put(target.sub("torque_hydro.x"),p[j].torque_hydro[0]);
    ufres->put(target.sub("torque_hydro.y"),p[j].torque_hydro[1]);
    ufres->put(target.sub("torque_hydro.z"),p[j].torque_hydro[2]);

    //old hydrodynamic torque
    ufres->put(target.sub("torque_hydro_previous.x"),p[j].torque_hydro_previous[0]);
    ufres->put(target.sub("torque_hydro_previous.y"),p[j].torque_hydro_previous[1]);
    ufres->put(target.sub("torque_hydro_previous.z"),p[j].torque_hydro_previous[2]);

    //electrolyte torque
    ufres->put(target.sub("torque_hydro1.x"),p[j].torque_hydro1[0]);
    ufres->put(target.sub("torque_hydro1.y"),p[j].torque_hydro1[1]);
    ufres->put(target.sub("torque_hydro1.z"),p[j].torque_hydro1[2]);

    //slip torque
    ufres->put(target.sub("torque_slip.x"), p[j].torque_slip[0]);
    ufres->put(target.sub("torque_slip.y"), p[j].torque_slip[1]);
    ufres->put(target.sub("torque_slip.z"), p[j].torque_slip[2]);

    //old slip torque
    ufres->put(target.sub("torque_slip_previous.x"), p[j].torque_slip_previous[0]);
    ufres->put(target.sub("torque_slip_previous.y"), p[j].torque_slip_previous[1]);
    ufres->put(target.sub("torque_slip_previous.z"), p[j].torque_slip_previous[2]);
  }
}

void Force_restore_parameters(double **zeta
			      ,double *uk_dc
			      ,Particle *p
			      ,CTime &time
			      ,double **conc_k
                              ){
  {//fluid data
    int im;
    for(int i=0; i<NX; i++){
      for(int j=0; j<NY; j++){
        for(int k=0; k<NZ_; k++){
          im = (i*NY*NZ_) + (j*NZ_) + k;
          char str[256];
          {
            sprintf(str,"resume.CONTINUE.Saved_Data.Zeta[%d][%d][%d]",i,j,k);
            Location target(str);
            ufin->get(target.sub("zeta0"),zeta[0][im]);
            ufin->get(target.sub("zeta1"),zeta[1][im]);
          }
          if(Electrolyte ){
            for(int n=0;n<N_spec;n++){ // Two_fluid では N_spec =1
              sprintf(str,"resume.CONTINUE.Saved_Data.Concentration[%d][%d][%d][%d]",n,i,j,k);
              Location target(str);
              ufin->get(target.sub("ck"),conc_k[n][im]);
            }
          }
        }
      }
    }
  }
  
  {//uk_dc
    
    char str[256];
    sprintf(str,"resume.CONTINUE.Saved_Data.uk_dc");
    Location target(str);
    
    ufin->get(target.sub("x"),uk_dc[0]);
    ufin->get(target.sub("y"),uk_dc[1]);
    ufin->get(target.sub("z"),uk_dc[2]);
  }

  //particle_data
  if(SW_PT != rigid){
    Read_Particle_udf(p, Particle_Number);
  }else{
    Particle *rigid_p = new Particle [Rigid_Number];
    Read_Particle_udf(rigid_p, Rigid_Number);
    Read_Rigid_Particle_udf();
    Set_Rigid_Particle_Data(rigid_p, p);
    delete [] rigid_p;
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
  
  if(SW_EQ == Shear_Navier_Stokes_Lees_Edwards){
    int im;
#pragma omp parallel for schedule(dynamic, 1) private(im)
    for(int i = 0; i < NX; i++){
      for(int j = 0; j < NY; j++){
        for(int k = 0; k < NZ; k++){
          im = (i*NY*NZ_) + (j*NZ_) + k;
          
          K2[im] = 
            SQ(WAVE_X*KX_int[im]) +
            SQ(WAVE_Y*KY_int[im] -
               WAVE_X*degree_oblique*KX_int[im]) +
            SQ(WAVE_Z*KZ_int[im]);
          if(K2[im] > 0.0){
            IK2[im] = 1.0/K2[im];
          }else{
            IK2[im] = 0.0;
          }
        }
      }
    }
  }

}      

void Read_Particle_udf(Particle *p, const int &n_in_particles){
  
  for(int j=0;j<n_in_particles;j++){
    char str[256];
    sprintf(str,"resume.CONTINUE.Saved_Data.Particles[%d]",j);
    Location target(str);

    //positions
    ufin->get(target.sub("R.x"),p[j].x[0]);
    ufin->get(target.sub("R.y"),p[j].x[1]);
    ufin->get(target.sub("R.z"),p[j].x[2]);

    //old/raw positions (not saved)
    for(int d = 0; d < DIM; d++){
      p[j].x_previous[d] = p[j].x[d];
      p[j].x_nopbc[d] = p[j].x[d];
    }
    
    //velocities
    ufin->get(target.sub("v.x"),p[j].v[0]);
    ufin->get(target.sub("v.y"),p[j].v[1]);
    ufin->get(target.sub("v.z"),p[j].v[2]);
    
    //old velocities
    ufin->get(target.sub("v_old.x"),p[j].v_old[0]);
    ufin->get(target.sub("v_old.y"),p[j].v_old[1]);
    ufin->get(target.sub("v_old.z"),p[j].v_old[2]);

    //slip velocities (not saved)
    for(int d = 0; d < DIM; d++){
      p[j].v_slip[d] = 0.0;
    }
    
    {//orientations
      double q0, q1, q2, q3;
      ufin->get(target.sub("q.q0"), q0);
      ufin->get(target.sub("q.q1"), q1);
      ufin->get(target.sub("q.q2"), q2);
      ufin->get(target.sub("q.q3"), q3);
      qtn_init(p[j].q, q0, q1, q2, q3);
      qtn_isnormal(p[j].q, HUGE_TOL_MP);
      qtn_normalize(p[j].q);
      
      ufin->get(target.sub("q_old.q0"), q0);
      ufin->get(target.sub("q_old.q1"), q1);
      ufin->get(target.sub("q_old.q2"), q2);
      ufin->get(target.sub("q_old.q3"), q3);
      qtn_init(p[j].q_old, q0, q1, q2, q3);
      qtn_isnormal(p[j].q_old, HUGE_TOL_MP);
      qtn_normalize(p[j].q_old);
    }
    
    //hydrodynamic force
    ufin->get(target.sub("f_hydro.x"),p[j].f_hydro[0]);
    ufin->get(target.sub("f_hydro.y"),p[j].f_hydro[1]);
    ufin->get(target.sub("f_hydro.z"),p[j].f_hydro[2]);

    //old hydrodynamic force
    ufin->get(target.sub("f_hydro_previous.x"),p[j].f_hydro_previous[0]);
    ufin->get(target.sub("f_hydro_previous.y"),p[j].f_hydro_previous[1]);
    ufin->get(target.sub("f_hydro_previous.z"),p[j].f_hydro_previous[2]);

    //electrolyte force
    ufin->get(target.sub("f_hydro1.x"),p[j].f_hydro1[0]);
    ufin->get(target.sub("f_hydro1.y"),p[j].f_hydro1[1]);
    ufin->get(target.sub("f_hydro1.z"),p[j].f_hydro1[2]);

    //slip force
    ufin->get(target.sub("f_slip.x"), p[j].f_slip[0]);
    ufin->get(target.sub("f_slip.y"), p[j].f_slip[1]);
    ufin->get(target.sub("f_slip.z"), p[j].f_slip[2]);

    //old slip force
    ufin->get(target.sub("f_slip_previous.x"), p[j].f_slip_previous[0]);
    ufin->get(target.sub("f_slip_previous.y"), p[j].f_slip_previous[1]);
    ufin->get(target.sub("f_slip_previous.z"), p[j].f_slip_previous[2]);
    
    //random+lj+ext force
    ufin->get(target.sub("fr.x"),p[j].fr[0]);
    ufin->get(target.sub("fr.y"),p[j].fr[1]);
    ufin->get(target.sub("fr.z"),p[j].fr[2]);

    //old random+lj+ext force
    ufin->get(target.sub("fr_previous.x"),p[j].fr_previous[0]);
    ufin->get(target.sub("fr_previous.y"),p[j].fr_previous[1]);
    ufin->get(target.sub("fr_previous.z"),p[j].fr_previous[2]);
    
    //angular velocity
    ufin->get(target.sub("omega.x"),p[j].omega[0]);
    ufin->get(target.sub("omega.y"),p[j].omega[1]);
    ufin->get(target.sub("omega.z"),p[j].omega[2]);

    //old angular velocity
    ufin->get(target.sub("omega_old.x"),p[j].omega_old[0]);
    ufin->get(target.sub("omega_old.y"),p[j].omega_old[1]);
    ufin->get(target.sub("omega_old.z"),p[j].omega_old[2]);

    //slip angular velocity (not saved)
    for(int d = 0; d < DIM; d++){
      p[j].omega_slip[d] = 0.0;
    }

    // hydrodynamic torque
    ufin->get(target.sub("torque_hydro.x"),p[j].torque_hydro[0]);
    ufin->get(target.sub("torque_hydro.y"),p[j].torque_hydro[1]);
    ufin->get(target.sub("torque_hydro.z"),p[j].torque_hydro[2]);

    //old hydrodynamic torque
    ufin->get(target.sub("torque_hydro_previous.x"),p[j].torque_hydro_previous[0]);
    ufin->get(target.sub("torque_hydro_previous.y"),p[j].torque_hydro_previous[1]);
    ufin->get(target.sub("torque_hydro_previous.z"),p[j].torque_hydro_previous[2]);

    //electrolyte torque
    ufin->get(target.sub("torque_hydro1.x"),p[j].torque_hydro1[0]);
    ufin->get(target.sub("torque_hydro1.y"),p[j].torque_hydro1[1]);
    ufin->get(target.sub("torque_hydro1.z"),p[j].torque_hydro1[2]);

    //slip torque
    ufin->get(target.sub("torque_slip.x"), p[j].torque_slip[0]);
    ufin->get(target.sub("torque_slip.y"), p[j].torque_slip[1]);
    ufin->get(target.sub("torque_slip.z"), p[j].torque_slip[2]);
    
    //old slip torque
    ufin->get(target.sub("torque_slip_previous.x"), p[j].torque_slip_previous[0]);
    ufin->get(target.sub("torque_slip_previous.y"), p[j].torque_slip_previous[1]);
    ufin->get(target.sub("torque_slip_previous.z"), p[j].torque_slip_previous[2]);


    {    //temp data (not saved)
      p[j].mass = 0.0;
      p[j].surface_mass = 0.0;
      for(int d = 0; d < DIM; d++){
        p[j].mass_center[d] = 0.0;
        p[j].surface_mass_center[d] = 0.0;
        p[j].surface_dv[d] = 0.0;
        p[j].surface_dw[d] = 0.0;
      }
      for(int d = 0; d < DIM; d++){
        for(int l = 0; l < DIM; l++){
          p[j].inertia[d][l] = 0.0;
          p[j].surface_inertia[d][l] = 0.0;
        }
      }
    }
  }
}

void Read_Rigid_Particle_udf(){
  if(SW_PT == rigid){
    for(int j = 0; j < Particle_Number; j++){
      char str[256];
      sprintf(str, "resume.CONTINUE.Saved_Data.GR_body[%d]", j);
      Location target(str);

      ufin->get(target.sub("x"), GRvecs_body[j][0]);
      ufin->get(target.sub("y"), GRvecs_body[j][1]);
      ufin->get(target.sub("z"), GRvecs_body[j][2]);

    }
  }
}

void Set_Rigid_Particle_Data(Particle *rigid_p, Particle *p){

  //rigid particles
#pragma omp parallel for schedule(dynamic, 1)
  for(int rigidID = 0; rigidID < Rigid_Number; rigidID++){
    int rigid_first_n = Rigid_Particle_Cumul[rigidID];
    int rigid_last_n = Rigid_Particle_Cumul[rigidID];

    for(int d = 0; d < DIM; d++){
      //positions
      xGs[rigidID][d] = rigid_p[rigidID].x[d];
      xGs_previous[rigidID][d] = rigid_p[rigidID].x_previous[d];
      xGs_nopbc[rigidID][d] = rigid_p[rigidID].x_nopbc[d];

      //velocities
      velocityGs[rigidID][d] = rigid_p[rigidID].v[d];
      velocityGs_old[rigidID][d] = rigid_p[rigidID].v_old[d];

      //angular velocities
      omegaGs[rigidID][d] = rigid_p[rigidID].omega[d];
      omegaGs_old[rigidID][d] = rigid_p[rigidID].omega_old[d];

      //hydrodynamic forces
      forceGs[rigidID][d] = rigid_p[rigidID].f_hydro[d];
      forceGs_previous[rigidID][d] = rigid_p[rigidID].f_hydro_previous[d];

      //ext+lj+random forces
      forceGrs[rigidID][d] = rigid_p[rigidID].fr[d];
      forceGrs_previous[rigidID][d] = rigid_p[rigidID].fr_previous[d];

      //hydrodynamic torque
      torqueGs[rigidID][d] = rigid_p[rigidID].torque_hydro[d];
      torqueGs_previous[rigidID][d] = rigid_p[rigidID].torque_hydro_previous[d];
    }
  }

  //bead particle data
#pragma omp parallel for schedule(dynamic, 1)
  for(int n = 0; n < Particle_Number; n++){
    int rigidID = Particle_RigidID[n];

    //orientations
    qtn_init(p[n].q, rigid_p[rigidID].q);
    qtn_init(p[n].q_old, rigid_p[rigidID].q_old);

    //positions & velocities
    rigid_body_rotation(GRvecs[n], GRvecs_body[n], p[n].q, BODY2SPACE);
    for(int d = 0; d < DIM; d++){
      p[n].x[d] = xGs[rigidID][d] + GRvecs[n][d];
      p[n].omega[d] = omegaGs[rigidID][d];

    }
    rigid_Velocity(p[n].v, GRvecs[n], velocityGs[rigidID], omegaGs[rigidID]);

    //periodic boundary conditions
    int sign;
    double delta_vx;
    if(SW_EQ != Shear_Navier_Stokes_Lees_Edwards){
      PBC(p[n].x);
    }else{
      PBC_OBL(p[n].x, delta_vx);
      p[n].v[0] += delta_vx;
    }

    //old data
    for(int d = 0; d < DIM; d++){
      p[n].x_previous[d] = p[n].x[d];
      p[n].x_nopbc[d] = p[n].x[d];

      p[n].v_old[d] = p[n].v[d];
      p[n].omega_old[d] = p[n].omega[d];
    }

    //unused data
    for(int d = 0; d < DIM; d++){
      {
        p[n].v_slip[d] = 0.0;
        
        p[n].f_hydro[d] = 0.0;
        p[n].f_hydro_previous[d] = 0.0;
        p[n].f_hydro1[d] = 0.0;
        p[n].f_slip[d] = 0.0;
        p[n].f_slip_previous[d] = 0.0;
        
        p[n].fr[d] = 0.0;
        p[n].fr_previous[d] = 0.0;
      }
      {
        p[n].omega_slip[d] = 0.0;
        
        p[n].torque_hydro[d] = 0.0;
        p[n].torque_hydro_previous[d] = 0.0;
        p[n].torque_hydro1[d] = 0.0;
        p[n].torque_slip[d] = 0.0;
        p[n].torque_slip_previous[d] = 0.0;

        p[n].momentum_depend_fr[d] = 0.0;
      }
    }
    for(int d = 0; d < DIM; d++){
      for(int l = 0; l < DIM; l++){
        p[n].QR[d][l] = 0.0;
        p[n].QR_old[d][l] = 0.0;
      }
    }

    //temp data
    p[n].mass = 0.0;
    p[n].surface_mass = 0.0;
    for(int d = 0; d < DIM; d++){
      p[n].mass_center[d] = 0.0;
      p[n].surface_mass_center[d] = 0.0;
      p[n].surface_dv[d] = 0.0;
      p[n].surface_dw[d] = 0.0;
    }
    for(int d = 0; d < DIM; d++){
      for(int l = 0; l < DIM; l++){
        p[n].inertia[d][l] = 0.0;
        p[n].surface_inertia[d][l] = 0.0;
      }
    }
  }
}
