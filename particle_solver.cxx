/*!
  \file particle_solver.cxx
  \brief Solver routines for particle position and velocity
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
 */

#include "particle_solver.h"

inline void reset_Forces(Particle *p){
#pragma omp parallel for schedule(dynamic, 1)
  for(int n = 0; n < Particle_Number; n++){
    for(int d = 0; d < DIM; d++){
      {
        p[n].fr_previous[d] = p[n].fr[d];
        p[n].fr[d] = 0.0;
        
        p[n].f_hydro_previous[d] = p[n].f_hydro[d];
        p[n].f_hydro[d] = 0.0;
        
        p[n].f_slip_previous[d] = p[n].f_slip[d];
        p[n].f_slip[d] = 0.0;
      }
      {
        p[n].torque_hydro_previous[d] = p[n].torque_hydro[d];
        p[n].torque_hydro[d] = 0.0;
        
        p[n].torque_slip_previous[d] = p[n].torque_slip[d];
        p[n].torque_slip[d] = 0.0;
      }
    }
  }
}

void MD_solver_position_Euler(Particle *p, const CTime &jikan)
{
  if(SW_PT != rigid){
    double delta_x;
#pragma omp parallel for schedule(dynamic, 1) private(delta_x)
    for(int n = 0; n < Particle_Number; n++) {
      
      if(janus_propulsion[p[n].spec] != obstacle){
        for(int d = 0; d < DIM; d++) {
          p[n].x_previous[d] = p[n].x[d];
          
          delta_x = jikan.dt_md * p[n].v[d];
          p[n].x_nopbc[d] += delta_x;
          p[n].x[d] += delta_x;
        }
        PBC(p[n].x);
        
        if(ROTATION)
          MD_solver_orientation_Euler(p[n], jikan.dt_md);
      }
    }
  }else{
    solver_Rigid_Position(p, jikan, "Euler");
    update_Particle_Configuration(p);
  }
}

void MD_solver_position_AB2(Particle *p, const CTime &jikan)
{
  if(SW_PT != rigid){
    double delta_x;
#pragma omp parallel for schedule(dynamic, 1) private(delta_x)
    for(int n = 0; n < Particle_Number; n++) {
      
      if(janus_propulsion[p[n].spec] != obstacle){
        for(int d = 0; d < DIM; d++) {
          p[n].x_previous[d] = p[n].x[d];
          
          delta_x = jikan.hdt_md * (3.0 * p[n].v[d] - p[n].v_old[d]);
          p[n].x_nopbc[d] += delta_x;
          p[n].x[d] += delta_x;
        }
        PBC(p[n].x);
        
        if(ROTATION)
          MD_solver_orientation_AB2(p[n], jikan.hdt_md);
      }
    }
  }else{
    solver_Rigid_Position(p, jikan, "AB2");
    update_Particle_Configuration(p);
  }
}

// self-force and torque in space coordinates 
inline void self_propulsion(Particle &p, double *force, double *torque){
  for(int d = 0; d < DIM; d++){
    force[d] = 0.0;
    torque[d] = 0.0;
  }
  
  if(SW_JANUS_MOTOR){
    int &spec = p.spec;
    if(janus_propulsion[spec] == motor){
      rigid_body_rotation(force, janus_force[spec], p.q, BODY2SPACE);
      rigid_body_rotation(torque, janus_torque[spec], p.q, BODY2SPACE);
    }
  }
}

void MD_solver_velocity_Euler(Particle *p, const CTime &jikan)
{
  Force(p);
  
  if(SW_PT != rigid){
    double dmy;
    double dmy_rot;
    double self_force[DIM];
    double self_torque[DIM];
    
#pragma omp parallel for schedule(dynamic,1) private(dmy, dmy_rot, self_force, self_torque)
    for(int n = 0; n < Particle_Number; n++){
      dmy = jikan.dt_md * IMASS[p[n].spec];
      dmy_rot = jikan.dt_md * IMOI[p[n].spec];
      
      if(janus_propulsion[p[n].spec] != obstacle){
        self_propulsion(p[n], self_force, self_torque);
        for(int d = 0; d < DIM; d++) {
          
          p[n].v_old[d] = p[n].v[d];
          p[n].omega_old[d] = p[n].omega[d];
          
          p[n].v[d] += ( dmy * (p[n].f_hydro[d] + p[n].fr[d] + self_force[d]) );
          p[n].omega[d] += ( dmy_rot * (p[n].torque_hydro[d] + self_torque[d]));
        }      
      }else{
        for(int d = 0; d < DIM; d++){
          p[n].v_old[d] = p[n].v[d];
          p[n].v[d] = 0.0;
          
          p[n].omega_old[d] = p[n].omega[d];
          p[n].omega[d] = 0.0;
        }
      }             
    }// Particle_Number
  }else{
    calc_Rigid_VOGs(p, jikan, "Euler");
    set_Particle_Velocities(p);
  }
  
  reset_Forces(p);
}

void MD_solver_velocity_slip_Euler(Particle *p, const CTime &jikan){
  double dmy;
  double dmy_rot;
  double self_force[DIM];
  double self_torque[DIM];
#pragma omp parallel for schedule(dynamic, 1) private(dmy, dmy_rot, self_force, self_torque)
  for(int n = 0; n < Particle_Number; n++){
    dmy = jikan.dt_md * IMASS[p[n].spec];
    dmy_rot = jikan.dt_md * IMOI[p[n].spec];
    
    if(janus_propulsion[p[n].spec] != obstacle){
      self_propulsion(p[n], self_force, self_torque);
      for(int d = 0; d < DIM; d++){
        p[n].v_old[d] = p[n].v[d];
        p[n].v[d] += (dmy * (p[n].f_hydro[d] + p[n].f_slip[d] + p[n].fr[d] + self_force[d]));
        
        p[n].omega_old[d] = p[n].omega[d];
        p[n].omega[d] += (dmy_rot * (p[n].torque_hydro[d] + p[n].torque_slip[d] + self_torque[d]));
      }
    }else{
      for(int d = 0; d < DIM; d++){
        p[n].v_old[d] = p[n].v[d];
        p[n].v[d] = 0.0;
        
        p[n].omega_old[d] = p[n].omega[d];
        p[n].omega[d] = 0.0;
      }
    }
  }//Particle_Number
}

void MD_solver_velocity_slip_AB2(Particle *p, const CTime &jikan){
  double dmy;
  double dmy_rot;
  double self_force[DIM];
  double self_torque[DIM];
#pragma omp parallel for schedule(dynamic, 1) private(dmy, dmy_rot, self_force, self_torque)
  for(int n = 0; n < Particle_Number; n++){
    dmy = jikan.hdt_md * IMASS[p[n].spec];
    dmy_rot = jikan.hdt_md * IMOI[p[n].spec];
    
    if(janus_propulsion[p[n].spec] != obstacle){
      self_propulsion(p[n], self_force, self_torque);
      for(int d = 0; d < DIM; d++){
        p[n].v_old[d] = p[n].v[d];
        p[n].v[d] += dmy * (2.0 * (p[n].f_hydro[d] + p[n].f_slip[d] + self_force[d]) 
                            + p[n].fr[d] + p[n].fr_previous[d]);

        
        p[n].omega_old[d] = p[n].omega[d];
        p[n].omega[d] += dmy_rot * (2.0 * (p[n].torque_hydro[d] + p[n].torque_slip[d] 
                                           + self_torque[d]));
      }
    }else{
      for(int d = 0; d < DIM; d++){
        p[n].v_old[d] = p[n].v[d];
        p[n].v[d] = 0.0;
        
        p[n].omega_old[d] = p[n].omega[d];
        p[n].omega[d] = 0.0;
      }
    }
  }//Particle_Number
}

void MD_solver_velocity_slip_iter(Particle *p, const CTime &jikan, 
                                  const ITER &iter_flag){
  if(iter_flag != start_iter && iter_flag != new_iter && 
     iter_flag != end_iter && iter_flag != reset_iter){
    fprintf(stderr, "Error: wrong Euler iter_flag\n");
    exit_job(EXIT_FAILURE);
  }
  if(iter_flag == start_iter){//only compute external forces once
    Force(p);
  }
  
  if(iter_flag == start_iter || iter_flag == new_iter){
    
    if(jikan.ts != 0){
      MD_solver_velocity_slip_AB2(p, jikan);
    }else{
      MD_solver_velocity_slip_Euler(p, jikan);
    }
    
  }else if(iter_flag == reset_iter){
#pragma omp parallel for schedule(dynamic, 1)
    for(int n = 0; n < Particle_Number; n++){
      for(int d = 0; d < DIM; d++){
	p[n].v[d] = p[n].v_old[d];
	p[n].omega[d] = p[n].omega_old[d];
      }
    }
    
  }else if(iter_flag == end_iter){
    reset_Forces(p);
  }//end_iter
}

void MD_solver_velocity_AB2_hydro(Particle *p, const CTime &jikan){
  Force(p);
  
  if(SW_PT != rigid){
    double dmy;
    double dmy_rot;
    double self_force[DIM];
    double self_torque[DIM];
    
#pragma omp parallel for schedule(dynamic,1) private(dmy, dmy_rot, self_force, self_torque)
    for(int n = 0; n < Particle_Number; n++) {
      
      dmy = jikan.hdt_md * IMASS[p[n].spec];
      dmy_rot = jikan.hdt_md * IMOI[p[n].spec];
      
      if(janus_propulsion[p[n].spec] != obstacle){
        self_propulsion(p[n], self_force, self_torque);
        for(int d = 0; d < DIM; d++){
          p[n].v_old[d] = p[n].v[d];
          p[n].omega_old[d] = p[n].omega[d];
          
          p[n].v[d] += dmy * ( 2.0 * (p[n].f_hydro[d] + self_force[d]) 
                               + p[n].fr[d] + p[n].fr_previous[d]); // CN
          p[n].omega[d] += dmy_rot * 2.0 *(p[n].torque_hydro[d] + self_torque[d]);
        }
      }else{
        for(int d = 0; d < DIM; d++){
          p[n].v_old[d] = p[n].v[d];
          p[n].v[d] = 0.0;
          
          p[n].omega_old[d] = p[n].omega[d];
          p[n].omega[d] = 0.0;
        }
      }
      
    }// Particle_Number
  }else{
    calc_Rigid_VOGs(p, jikan, "AB2");
    set_Particle_Velocities(p);
  }
  
  reset_Forces(p);
}
	

void MD_solver_position_Euler_OBL(Particle *p, const CTime &jikan){
  if(SW_PT != rigid){
    double delta_x, delta_vx;
#pragma omp parallel for schedule(dynamic, 1) private(delta_x, delta_vx)
    for(int n = 0; n < Particle_Number; n++) {
      if(janus_propulsion[p[n].spec] != obstacle){
        for(int d = 0; d < DIM; d++) {
          p[n].x_previous[d] = p[n].x[d];
          
          delta_x = jikan.dt_md * p[n].v[d];
          p[n].x_nopbc[d] += delta_x;
          p[n].x[d] += delta_x;
          
        }
        
        int sign = PBC_OBL(p[n].x, delta_vx);
        p[n].v[0] += delta_vx;
        p[n].v_old[0] += delta_vx;
        
        if(ROTATION)
          MD_solver_orientation_Euler(p[n], jikan.dt_md);
      }
    }
  }else{
    solver_Rigid_Position_OBL(p, jikan, "Euler");
    update_Particle_Configuration_OBL(p);
  }
}
void MD_solver_position_AB2_OBL(Particle *p, const CTime &jikan){
  if(SW_PT != rigid){
    double delta_x, delta_vx;
#pragma omp parallel for schedule(dynamic, 1) private(delta_x, delta_vx)
    for(int n = 0; n < Particle_Number; n++) {
      if(janus_propulsion[p[n].spec] != obstacle){
        for(int d = 0; d < DIM; d++) {
          p[n].x_previous[d] = p[n].x[d];
          
          delta_x = jikan.hdt_md * (3.0 * p[n].v[d] - p[n].v_old[d]);
          p[n].x_nopbc[d] += delta_x;
          p[n].x[d] += delta_x;
        }
        
        int sign = PBC_OBL(p[n].x, delta_vx);
        p[n].v[0] += delta_vx;
        p[n].v_old[0] += delta_vx;
        
        if(ROTATION)
          MD_solver_orientation_AB2(p[n], jikan.hdt_md);
      }
    }
  }else{
    solver_Rigid_Position_OBL(p, jikan, "AB2");
    update_Particle_Configuration_OBL(p);
  }
}

void MD_solver_velocity_Euler_OBL(Particle *p, const CTime &jikan){
  Force_OBL(p);

  if(SW_PT != rigid){
    double dmy;
    double dmy_rot;
    double self_force[DIM];
    double self_torque[DIM];
    
#pragma omp parallel for schedule(dynamic,1) private(dmy, dmy_rot, self_force, self_torque)
    for(int n = 0; n < Particle_Number; n++) {
      if(janus_propulsion[p[n].spec] != obstacle){
        dmy = jikan.dt_md * IMASS[p[n].spec];
        dmy_rot = jikan.dt_md * IMOI[p[n].spec];
        self_propulsion(p[n], self_force, self_torque);
        
        for(int d = 0; d < DIM; d++) {
          {
            p[n].v_old[d] = p[n].v[d];
            p[n].omega_old[d] = p[n].omega[d];
            
            p[n].v[d] += ( dmy * ( p[n].f_hydro[d] + p[n].fr[d] + self_force[d]) );
            p[n].omega[d] += ( dmy_rot * (p[n].torque_hydro[d] + self_torque[d]) );
            
            //hydro_stress calculations
            p[n].momentum_depend_fr[d] = jikan.dt_md * p[n].fr[d];
          }
        }
      }else{
        for(int d = 0; d < DIM; d++){
          p[n].v_old[d] = p[n].v[d];
          p[n].v[d] = 0.0;
          
          p[n].omega_old[d] = p[n].omega[d];
          p[n].omega[d] = 0.0;

          p[n].momentum_depend_fr[d] = jikan.dt_md * p[n].fr[d];
        }
      }
    }
  }else{
    calc_Rigid_VOGs(p, jikan, "Euler");    
    set_Particle_Velocities_OBL(p);
  }

  reset_Forces(p);
}

void MD_solver_velocity_AB2_hydro_OBL(Particle *p, const CTime &jikan){
  Force_OBL(p);

  if(SW_PT != rigid){
    double dmy;
    double dmy_rot;
    double self_force[DIM];
    double self_torque[DIM];
    
#pragma omp parallel for schedule(dynamic,1) private(dmy, dmy_rot, self_force, self_torque)
    for(int n = 0; n < Particle_Number; n++) {
      if(janus_propulsion[p[n].spec] != obstacle){
        dmy = jikan.hdt_md * IMASS[p[n].spec];
        dmy_rot = jikan.hdt_md * IMOI[p[n].spec];
        self_propulsion(p[n], self_force, self_torque);
        
        for(int d = 0; d < DIM; d++) {
          {
            p[n].v_old[d] = p[n].v[d];
            p[n].omega_old[d] = p[n].omega[d];
            
            p[n].v[d] += dmy * (2.0*(p[n].f_hydro[d] + self_force[d]) 
                                + p[n].fr[d] + p[n].fr_previous[d]); // CN
            p[n].omega[d] += dmy_rot * 2.0 *(p[n].torque_hydro[d] + self_torque[d]);
            
            //hydro_stress calculations
            p[n].momentum_depend_fr[d] = jikan.hdt_md * (p[n].fr[d] + p[n].fr_previous[d]);
          }
        }
      }else{
        for(int d = 0; d < DIM; d++){
          p[n].v_old[d] = p[n].v[d];
          p[n].v[d] = 0.0;

          p[n].omega_old[d] = p[n].omega[d];
          p[n].omega[d] = 0.0;

          p[n].momentum_depend_fr[d] = jikan.hdt_md * (p[n].fr[d] + p[n].fr_previous[d]);
        }
      }
    }
  }else{
    calc_Rigid_VOGs(p, jikan, "AB2");
    set_Particle_Velocities_OBL(p);
  }

  reset_Forces(p);
}
