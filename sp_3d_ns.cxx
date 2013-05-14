//
// $Id: sp_3d_ns.cxx,v 1.64 2006/06/11 16:15:08 nakayama Exp $
//
#include "sp_3d_ns.h"

void (*Time_evolution)(Value zeta[DIM-1], double uk_dc[DIM], Value f[DIM], Particle *p, CTime &jikan);

void Time_interval(const Particle *p, CTime &jikan){
  static double vmax;
  vmax = 0.0;
  for(int n=0; n<Particle_Number; n++){
    for(int d=0; d<DIM; d++){  
	vmax = MAX(ABS(p[n].v[d]), vmax);
    }
  }
  if(vmax == 0.0 ){
    jikan.dt_md = jikan.dt_fluid = DT;
    jikan.hdt_md = jikan.hdt_fluid = jikan.dt_fluid *.5;
    return;
  }
  static double Interface_Stokes_time;
  Interface_Stokes_time = XI/vmax;
  if(0){
    static double factor = 1./20.;
    if(Interface_Stokes_time * factor < DT){
      jikan.dt_fluid = Interface_Stokes_time * factor;
    }else{
      jikan.dt_fluid = DT;
    }
    if(0){
      jikan.dt_fluid = Interface_Stokes_time * factor;
    }
    jikan.dt_md = jikan.dt_fluid;
    jikan.hdt_md = jikan.hdt_fluid = jikan.dt_fluid *.5;
  }
  if(1){
    fprintf(stderr, "%g %g %g %g %g %g %g\n"
	    ,p[0].v[0]
	    ,p[0].v[1]
	    ,p[0].v[2]
	    ,vmax
	    ,Interface_Stokes_time
	    ,jikan.dt_fluid
	    ,Interface_Stokes_time/jikan.dt_fluid);
  }
}

void Time_evolution_slippy_NS(Value zeta[DIM-1], double uk_dc[DIM], Value f[DIM], Particle *p, CTime &jikan){
  Time_interval(p, jikan);

  {// for two-thirds rule 
    const Index_range ijk_range[] = {
      {0,TRN_X-1, 0,TRN_Y-1, 0,2*TRN_Z-1}
      ,{0,TRN_X-1, NY-TRN_Y+1,NY-1,  0,2*TRN_Z-1}
      ,{NX-TRN_X+1,NX-1,  0,TRN_Y-1, 0,2*TRN_Z-1}
      ,{NX-TRN_X+1,NX-1,  NY-TRN_Y+1,NY-1,  0,2*TRN_Z-1}
    };
    const int n_ijk_range = sizeof(ijk_range)/sizeof(Index_range);
    Slippy_NS_solver_Euler(zeta, jikan, uk_dc, ijk_range, n_ijk_range, p);
  }

  if(FIX_CELL){ // time-dependent average pressure gradient
    for(int d=0;d<DIM;d++){
      if(FIX_CELLxyz[d]){
	uk_dc[d] = 0.0;
      }
    }
  }
  
  {// construct particle velocity
    for(int n=0; n<Particle_Number; n++){
      for(int d=0; d<DIM; d++){  
	p[n].v_old[d] = p[n].v[d];
	p[n].v[d] = 0.0;
      }
    }
    Reset_phi(phi);
    Make_phi_u_particle(phi, up, p);
    Zeta_k2u(zeta, uk_dc, u);
    Calc_particle_velocity(phi, u, p);
  }
}
void Time_evolution_nohydro(Value zeta[DIM-1], double uk_dc[DIM], Value f[DIM], Particle *p, CTime &jikan){
  {
    Calc_f_hydro_draining(p, f, jikan);
    if(HYDRO_int < 0){
      Calc_f_hydro_lubrication(p);
    }
  }
  if(jikan.ts == 0){
    MD_solver_velocity_Euler(p, jikan);
  }else{
    MD_solver_velocity_AB2(p, jikan);
    //MD_solver_velocity_sAB2_nohydro(p, jikan);
  }
  MD_solver_position_CN(p, jikan);
}
void Time_evolution_nohydro_shear(Value zeta[DIM-1], double uk_dc[DIM], Value f[DIM], Particle *p, CTime &jikan){
  {
    Calc_f_hydro_draining_shear(p, f, jikan);
  }
  Shear_strain += -Shear_rate * LY_shear * jikan.dt_fluid; 
  Shear_strain -= floor(Shear_strain*iL_particle[0])*L_particle[0];

  if(jikan.ts == 0){
    MD_solver_velocity_Euler(p, jikan);
  }else{
    //MD_solver_velocity_Euler(p, jikan);
    MD_solver_velocity_AB2(p, jikan);
    //MD_solver_velocity_sAB2_nohydro(p, jikan);
  }
  MD_solver_position_CN_shear(p, jikan);
}
void Time_evolution_nohydro_brown(Value zeta[DIM-1], double uk_dc[DIM], Value f[DIM], Particle *p, CTime &jikan){
  //BD_solver_position_Euler(p, jikan);
  BD_solver_momentum(p, jikan);
}

void Time_evolution_noparticle(Value zeta[DIM-1], double uk_dc[DIM], Value f[DIM], Particle *p, CTime &jikan){
  if(SW_EQ == Navier_Stokes){
    //for two-thirds rule 
    const Index_range ijk_range[] = {
      {0,TRN_X-1, 0,TRN_Y-1, 0,2*TRN_Z-1}
      ,{0,TRN_X-1, NY-TRN_Y+1,NY-1,  0,2*TRN_Z-1}
      ,{NX-TRN_X+1,NX-1,  0,TRN_Y-1, 0,2*TRN_Z-1}
      ,{NX-TRN_X+1,NX-1,  NY-TRN_Y+1,NY-1,  0,2*TRN_Z-1}
    };
    const int n_ijk_range = sizeof(ijk_range)/sizeof(Index_range);
    //Stokes_solver(zeta, jikan, uk_dc);
    if(kBT > 0.0){
      //NS_solver_Euler2(zeta, jikan, uk_dc, ijk_range, n_ijk_range, p);
      //NS_solver_Heun(zeta, jikan, uk_dc, ijk_range, n_ijk_range, p);
      NS_solver_slavedEuler(zeta, jikan, uk_dc, ijk_range, n_ijk_range, p);
      //NS_solver_slavedHeun(zeta, jikan, uk_dc, ijk_range, n_ijk_range, p);
    }else {
      //NS_solver_Euler(zeta, jikan, uk_dc, ijk_range, n_ijk_range, p);
      //NS_solver_Heun(zeta, jikan, uk_dc, ijk_range, n_ijk_range, p);
      NS_solver_slavedEuler(zeta, jikan, uk_dc, ijk_range, n_ijk_range, p);
    }
  }else if(SW_EQ == Qian_Sheng){
    // for four-fifths rule 
    const Index_range ijk_range[] = {
      {0,TRN_QS_X-1, 0,TRN_QS_Y-1, 0,2*TRN_QS_Z-1}
      ,{0,TRN_QS_X-1, NY-TRN_QS_Y+1,NY-1,  0,2*TRN_QS_Z-1}
      ,{NX-TRN_QS_X+1,NX-1,  0,TRN_QS_Y-1, 0,2*TRN_QS_Z-1}
      ,{NX-TRN_QS_X+1,NX-1,  NY-TRN_QS_Y+1,NY-1,  0,2*TRN_QS_Z-1}
    };
    const int n_ijk_range = sizeof(ijk_range)/sizeof(Index_range);
    QS_solver_Euler(zeta, jikan, uk_dc, ijk_range, n_ijk_range);
  }else if(SW_EQ == nematic_Allen_Cahn){
    // for four-seconds rule 
    const Index_range ijk_range[] = {
      {0,NX-1, 0,NY-1, 0,NZ-1}
//       {0,HNX-1, 0,HNY-1, 0,NZ-1}
//       ,{0,HNX-1, HNY+1,NY-1, 0,NZ-1}
//       ,{HNX+1,NX-1, 0,HNY-1, 0,NZ-1}
//       ,{HNX+1,NX-1, HNY+1,NY-1, 0,NZ-1}
    };
    const int n_ijk_range = sizeof(ijk_range)/sizeof(Index_range);
    AC_solver_Euler(jikan, p, ijk_range, n_ijk_range);
  }else if(SW_EQ == Olmsted_Goldbart){
    // for four-fifths rule 
    const Index_range ijk_range[] = {
      {0,TRN_QS_X-1, 0,TRN_QS_Y-1, 0,2*TRN_QS_Z-1}
      ,{0,TRN_QS_X-1, NY-TRN_QS_Y+1,NY-1,  0,2*TRN_QS_Z-1}
      ,{NX-TRN_QS_X+1,NX-1,  0,TRN_QS_Y-1, 0,2*TRN_QS_Z-1}
      ,{NX-TRN_QS_X+1,NX-1,  NY-TRN_QS_Y+1,NY-1,  0,2*TRN_QS_Z-1}
    };
    const int n_ijk_range = sizeof(ijk_range)/sizeof(Index_range);
    OG_solver_Euler(zeta, jikan, uk_dc, ijk_range, n_ijk_range);
  }else if(SW_EQ == Shear_Navier_Stokes){
    //for two-thirds rule 
    const Index_range ijk_range[] = {
      {0,TRN_X-1, 0,TRN_Y-1, 0,2*TRN_Z-1}
      ,{0,TRN_X-1, NY-TRN_Y+1,NY-1,  0,2*TRN_Z-1}
      ,{NX-TRN_X+1,NX-1,  0,TRN_Y-1, 0,2*TRN_Z-1}
      ,{NX-TRN_X+1,NX-1,  NY-TRN_Y+1,NY-1,  0,2*TRN_Z-1}
    };
    const int n_ijk_range = sizeof(ijk_range)/sizeof(Index_range);
    //Shear_NS_solver_Euler(zeta, jikan, uk_dc, ijk_range, n_ijk_range, p, Shear_force);
    Shear_NS_solver_Heun(zeta, jikan, uk_dc, ijk_range, n_ijk_range, p, Shear_force);

    {
      static const double iLx = 1./L[0];
      Shear_strain -= Shear_rate * LY_shear * jikan.dt_fluid; 
      double dmy = fmod(Shear_strain, L[0]);
      Shear_strain_int += (int)((Shear_strain - dmy) * iLx);
      Shear_strain = dmy;
    }
  }else if(SW_EQ == Electrolyte){
    //for two-thirds rule 
    const Index_range ijk_range[] = {
      {0,TRN_X-1, 0,TRN_Y-1, 0,2*TRN_Z-1}
      ,{0,TRN_X-1, NY-TRN_Y+1,NY-1,  0,2*TRN_Z-1}
      ,{NX-TRN_X+1,NX-1,  0,TRN_Y-1, 0,2*TRN_Z-1}
      ,{NX-TRN_X+1,NX-1,  NY-TRN_Y+1,NY-1,  0,2*TRN_Z-1}
    };
    const int n_ijk_range = sizeof(ijk_range)/sizeof(Index_range);
    //NS_solver_Euler(zeta, jikan, uk_dc, ijk_range, n_ijk_range, p);
    //////////////////
    //NSsolute_solver_Euler(zeta, jikan, uk_dc, Concentration,p,ijk_range, n_ijk_range);
    NSsolute_solver_Heun(zeta, jikan, uk_dc, Concentration,p,ijk_range, n_ijk_range);
    //////////////////
  }else if(SW_EQ == Two_fluid){
    //for two-thirds rule 
    const Index_range ijk_range[] = {
      {0,TRN_X-1, 0,TRN_Y-1, 0,2*TRN_Z-1}
      ,{0,TRN_X-1, NY-TRN_Y+1,NY-1,  0,2*TRN_Z-1}
      ,{NX-TRN_X+1,NX-1,  0,TRN_Y-1, 0,2*TRN_Z-1}
      ,{NX-TRN_X+1,NX-1,  NY-TRN_Y+1,NY-1,  0,2*TRN_Z-1}
    };
    const int n_ijk_range = sizeof(ijk_range)/sizeof(Index_range);
    {
      //NSsolute_solver_Euler(zeta, jikan, uk_dc, Concentration,p
      //  ,ijk_range, n_ijk_range);
      NSsolute_solver_Heun(zeta, jikan, uk_dc, Concentration,p
			   ,ijk_range, n_ijk_range);
    }
  }
}
void Time_evolution_hydro(Value zeta[DIM-1], double uk_dc[DIM], Value f[DIM], Particle *p, CTime &jikan){
  Time_evolution_noparticle(zeta, uk_dc, f, p, jikan);

  if(Particle_Number > 0){
    if(FIX_CELL){ // time-dependent average pressure gradient
      for(int d=0;d<DIM;d++){
	if(FIX_CELLxyz[d]){
	  uk_dc[d] = 0.0;
	}
      }
    }
    Zeta_k2u(zeta, uk_dc, u);
    
    if(SW_EQ == Navier_Stokes
       || SW_EQ == Slippy_Navier_Stokes 
       || SW_EQ == Shear_Navier_Stokes
       ){
      Reset_phi(Pressure);
      Value uu[DIM*2] = {f_particle[0]
			 ,f_particle[1]
			 ,f_particle[2]
			 ,up[0]
			 ,up[1]
			 ,up[2]};
      Add_advection_pressurek(Pressure, u, uu);
    }

    {
      {
	if(!Fixed_particle){
	  if(jikan.ts == 0){
	    MD_solver_position_Euler(p, jikan);
	  }else{
	    //	    MD_solver_position_Euler(p, jikan);
	    MD_solver_position_AB2(p, jikan);
	  }
	}
      }
      if(SW_EQ == Two_fluid
	 || SW_EQ == Electrolyte 
	 ){
	//double rescale_factor[N_spec];
	double * rescale_factor = new double[N_spec];
	Rescale_solute(rescale_factor
		       ,Total_solute
		       ,Concentration, p, phi, up[0]);
	delete [] rescale_factor;
      }

      if(SW_EQ == Electrolyte){
	{
	  Reset_phi_u(phi, up);
	  Make_phi_u_particle(phi, up, p);
	  Make_f_particle_dt_nonsole(f, u, up, phi);
	  Calc_f_hydro_correct(p, f, jikan);
	  for(int n=0;n<Particle_Number;n++){
	    for(int d=0;d<DIM;d++){
	      p[n].f_hydro1[d] = p[n].f_hydro[d];
	      p[n].torque_hydro1[d] = p[n].torque_hydro[d];
	    }
	  }
	}
	Make_Coulomb_force_x_on_fluid(f, p, Concentration, up[0], up[1], jikan);
	double dmy = jikan.dt_fluid * IRHO;
	for(int d=0;d<DIM;d++){
	  for(int i=0;i<NX;i++){
	    for(int j=0;j<NY;j++){
	      for(int k=0;k<NZ;k++){
		u[d][i][j][k] += (f[d][i][j][k] * dmy);
	      }
	    }
	  }
	}
	Solenoidal_u(u);
      }
      {
	Reset_phi_u(phi, up);
	Make_phi_u_particle(phi, up, p);
	Make_f_particle_dt_nonsole(f, u, up, phi);
	Calc_f_hydro_correct(p, f, jikan);
      }

      if(!Fixed_particle){
	if(jikan.ts == 0){
	  MD_solver_velocity_Euler(p, jikan);
	}else{
	  //	  MD_solver_velocity_Euler_hydro(p, jikan);
	  MD_solver_velocity_AB2_hydro(p, jikan);
	}
      }
    }
    {
      Reset_phi_u(phi, up);
      Make_phi_u_particle(phi, up, p);
      Make_f_particle_dt_nonsole(f, u, up, phi);
      U2u_k(f);
      if(SW_EQ == Navier_Stokes
	 || SW_EQ == Slippy_Navier_Stokes 
	 || SW_EQ == Shear_Navier_Stokes
	 ){
	U_k2divergence_k(f, phi);
	const double dmy = RHO/jikan.dt_fluid;
	for(int i=0; i<NX; i++){
	  for(int j=0; j<NY; j++){
	    for(int k=0; k<NZ_; k++){
	      phi[i][j][k] *= dmy;
	    }
	  }
	}
	Add_pressurek(Pressure, phi);
      }
      Solenoidal_uk(f);
      U_k2u(f);
      Add_f_particle(u, f);
    }
    U2zeta_k(zeta, uk_dc, u);
    if(SW_EQ == Shear_Navier_Stokes){
      Symmetrize_zetak(zeta);
    }
  }
} 
enum SW_NOISE {
  Before_NS
  ,After_R
  ,After_V
};
const SW_NOISE SW_noise = After_V;
//const SW_NOISE SW_noise = After_R;
//const SW_NOISE SW_noise = Before_NS;
void Time_evolution_noisyhydro(Value zeta[DIM-1]
			       ,double uk_dc[DIM]
			       ,Value f[DIM]
			       ,Particle *p
			       ,CTime &jikan){
  if(SW_noise == Before_NS){    //////////
    Zeta_k2u(zeta, uk_dc, u);
    {
      const Index_range ijk_range[] = {
	{0,TRN_X-1, 0,TRN_Y-1, 0,2*TRN_Z-1}
	,{0,TRN_X-1, NY-TRN_Y+1,NY-1,  0,2*TRN_Z-1}
	,{NX-TRN_X+1,NX-1,  0,TRN_Y-1, 0,2*TRN_Z-1}
	,{NX-TRN_X+1,NX-1,  NY-TRN_Y+1,NY-1,  0,2*TRN_Z-1}
      };
      const int n_ijk_range = sizeof(ijk_range)/sizeof(Index_range);
      const double truncate_factor= 1.0;
      Add_random_stress_x_slavedNS(u, p, ijk_range, n_ijk_range, jikan, truncate_factor, f_particle, up);
    }
    U2zeta_k(zeta, uk_dc, u);
  }
  Time_evolution_noparticle(zeta, uk_dc, f, p, jikan);

  if(Particle_Number > 0){
    if(FIX_CELL){ // time-dependent average pressure gradient
      for(int d=0;d<DIM;d++){
	if(FIX_CELLxyz[d]){
	  uk_dc[d] = 0.0;
	}
      }
    }
    Zeta_k2u(zeta, uk_dc, u);

    {
      if(1){
	if(jikan.ts == 0){
	  MD_solver_position_Euler(p, jikan);
	}else{
	  MD_solver_position_Euler(p, jikan);
	  //MD_solver_position_AB2(p, jikan);
	}
      }
      if(SW_noise == After_R){ //////////
	const Index_range ijk_range[] = {
	  {0,TRN_X-1, 0,TRN_Y-1, 0,2*TRN_Z-1}
	  ,{0,TRN_X-1, NY-TRN_Y+1,NY-1,  0,2*TRN_Z-1}
	  ,{NX-TRN_X+1,NX-1,  0,TRN_Y-1, 0,2*TRN_Z-1}
	  ,{NX-TRN_X+1,NX-1,  NY-TRN_Y+1,NY-1,  0,2*TRN_Z-1}
	};
	const int n_ijk_range = sizeof(ijk_range)/sizeof(Index_range);
	const double truncate_factor= 1.0;
	Add_random_stress_x_slavedNS(u, p, ijk_range, n_ijk_range, jikan, truncate_factor, f_particle, up);
      }
      {
	Reset_phi_u(phi, up);
	Make_phi_u_particle(phi, up, p);
	Make_f_particle_dt_nonsole(f, u, up, phi);
	Calc_f_hydro_correct(p, f, jikan);
      }
      if(jikan.ts == 0){
	MD_solver_velocity_Euler(p, jikan);
      }else{
	//BD_solver_velocity_hydro(p, jikan);
	MD_solver_velocity_Euler_hydro(p, jikan);
      }
      //Add_random_force(p, jikan);
    }
    {
      if(SW_noise == After_V
	 && SW_EQ == Navier_Stokes
	 ){
	const Index_range ijk_range[] = {
	  {0,TRN_X-1, 0,TRN_Y-1, 0,2*TRN_Z-1}
	  ,{0,TRN_X-1, NY-TRN_Y+1,NY-1,  0,2*TRN_Z-1}
	  ,{NX-TRN_X+1,NX-1,  0,TRN_Y-1, 0,2*TRN_Z-1}
	  ,{NX-TRN_X+1,NX-1,  NY-TRN_Y+1,NY-1,  0,2*TRN_Z-1}
	};
	// 	  const Index_range ijk_range[] = {
	// 	    {0,NX-1, 0,NY-1, 0,NZ_-1}
	// 	  };
	const int n_ijk_range = sizeof(ijk_range)/sizeof(Index_range);
	const double truncate_factor= 1.0;
	if(1){
	  Add_random_stress_x_up_slavedNS(up, p, ijk_range, n_ijk_range, jikan, truncate_factor, f_particle);
	}else {
	  Add_random_stress_x_slavedNS(u, p, ijk_range, n_ijk_range, jikan, truncate_factor, f_particle, up);
	  Reset_phi_u(phi, up);
	  Make_phi_u_particle(phi, up, p);
	}
      }else {
	Reset_phi_u(phi, up);
	Make_phi_u_particle(phi, up, p);
      }
      Make_f_particle_dt_sole(f, u, up, phi);
      Add_f_particle(u, f);
    }
    U2zeta_k(zeta, uk_dc, u);
  }
} 
inline void Mem_alloc_var(Value *zeta){
  if(SW_EQ == Navier_Stokes){
    Mem_alloc_NS_solver();
  }else if(SW_EQ == Shear_Navier_Stokes){
    Mem_alloc_NS_solver();
    for(int d=0;d<DIM;d++){
      Shear_force[d] = alloc_3d_double(NX, NY, NZ_);
    }
  }else if(SW_EQ==Qian_Sheng || SW_EQ == nematic_Allen_Cahn || SW_EQ == Olmsted_Goldbart){
    Mem_alloc_qij();
    Mem_alloc_QS();
  }else if(SW_EQ==Slippy_Navier_Stokes){
    Mem_alloc_NS_solver();
    Mem_alloc_QS();
  }else if(SW_EQ==Electrolyte){
    Mem_alloc_NS_solver();
    Mem_alloc_charge();
  }else if(SW_EQ == Two_fluid){
    Mem_alloc_NS_solver();
    Mem_alloc_solute();
  }
  Mem_alloc_f_particle();
  for(int d=0;d<DIM-1;d++){
    zeta[d] = alloc_3d_double(NX, NY, NZ_);
  }
  for(int d=0;d<DIM;d++){
    u[d] = alloc_3d_double(NX, NY, NZ_);
    up[d] = alloc_3d_double(NX, NY, NZ_);
  }
  phi = alloc_3d_double(NX, NY, NZ_);
}

int main(int argc, char *argv[]){

  if(argc> 0){
    file_get(argc, argv);
    Gourmet_file_io(In_udf,Out_udf,Sum_udf,Def_udf,Ctrl_udf,Res_udf);
  }
  if(HYDRO_int > 0){
    if(SW_EQ == Slippy_Navier_Stokes){
      Time_evolution = Time_evolution_slippy_NS;
    }else if(SW_EQ == nematic_Allen_Cahn){
      Time_evolution = Time_evolution_noparticle;
    }else{
      if(kBT > 0.0){
	if(SW_EQ == Electrolyte
	   || SW_EQ == Two_fluid
	   ){
	  Time_evolution = Time_evolution_hydro;
	}else {
	  Time_evolution = Time_evolution_noisyhydro;
	}
      }else{
	Time_evolution = Time_evolution_hydro;
      }
    }
  }else{
    if(kBT > 0.0){
      if(SW_EQ == Electrolyte){
	fprintf(stderr, "no implementatin of Electrolyte with HYDRO_int <= 0.\n");
	exit_job(EXIT_FAILURE);
      }else if(SW_EQ == Two_fluid){
	fprintf(stderr, "no implementatin of Two_fluid with HYDRO_int <= 0.\n");
	exit_job(EXIT_FAILURE);
      }else {
	Time_evolution = Time_evolution_nohydro_brown;
      }
    }else {
      if(SW_EQ == Shear_Navier_Stokes){
	Time_evolution = Time_evolution_nohydro_shear;
      }else {
	Time_evolution = Time_evolution_nohydro;
      }
    }
  }

  //MT_seed(GIVEN_SEED,0);
  MT_seed(RANDOM_SEED,0);
  Init_fft();
  double uk_dc[DIM];
  Value zeta[DIM-1];
  Mem_alloc_var(zeta);

  static CTime jikan={0, 0.0, DT, DT*0.5, DT, DT*0.5};

  Set_avs_parameters(Avs_parameters);
  if(SW_AVS){
    Init_avs(Avs_parameters);
    if(Particle_Number > 0){
      Init_avs_p(Avs_parameters);
    }
  }
  //Particle particles[Particle_Number];
  Particle *particles = new Particle [Particle_Number];
  if(Particle_Number > 0){
    Init_Particle(particles);
    /*
    if(SW_EQ == Shear_Navier_Stokes){
      if(Particle_Number == 1){
	particles[0].x[0] = HLX;
	particles[0].x[1] = HLY*.5;
	particles[0].x[2] = HLZ;
      } else if(Particle_Number == 2){
	double dmy_x = 1.4;
	double dmy_y = 1.4;
	particles[0].x[0] = HLX;// + dmy_x * RADIUS;//LX-RADIUS;
	particles[0].x[1] = L_particle[1]-RADIUS*dmy_y;
	particles[0].x[2] = HLZ;

	particles[1].x[0] = HLX;// - dmy_x * RADIUS; //0.;
	particles[1].x[1] = RADIUS*dmy_y;
	particles[1].x[2] = HLZ;
      }      
    }else if(SW_EQ == Electrolyte){
      if(Particle_Number == 1){
	particles[0].x[0] = HLX;
	particles[0].x[1] = HLY;
	particles[0].x[2] = HLZ;
	particles[0].v[0] = 0.;
	particles[0].v[1] = 0.;
	particles[0].v[2] = 0.;
	particles[0].omega[0] = 0.;
	particles[0].omega[1] = 0.;
	particles[0].omega[2] = 0.;
      }else if(Particle_Number == 2){
	double dmy = 1.1;//1.2;//2.0;
	particles[0].x[0] = HLX - RADIUS * dmy;
	particles[0].x[1] = HLY;
	particles[0].x[2] = HLZ;

	particles[1].x[0] = HLX + RADIUS * dmy;
	particles[1].x[1] = HLY;
	particles[1].x[2] = HLZ;
	//particles[0].v[0] = 10.0;
      }
    }else{
      if(Particle_Number == 1){
	particles[0].x[0] = HLX;
	particles[0].x[1] = HLY;
	particles[0].x[2] = HLZ;

	particles[0].v[2] = 1.e0;
	if(0){
	  //if(kBT>0 && SW_EQ == Navier_Stokes){
	  double dmy = 10.*sqrt(kBT * IMASS[particles[0].spec]);
	  particles[0].v[0] = dmy * 2.*(genrand_real3()-.5);
	  particles[0].v[1] = dmy * 2.*(genrand_real3()-.5);
	  particles[0].v[2] = dmy * 2.*(genrand_real3()-.5);
	  dmy = 10.*sqrt(kBT * IMOI[particles[0].spec]);
	  particles[0].omega[0] = dmy * 2.*(genrand_real3()-.5);
	  particles[0].omega[1] = dmy * 2.*(genrand_real3()-.5);
	  particles[0].omega[2] = dmy * 2.*(genrand_real3()-.5);
	}
      }else if(Particle_Number == 2){
	//double dmy = HLX*.5;
	double dmy = LJ_dia*pow(2.,1./6.)*.5;
	particles[0].x[0] = HLX - dmy;
	particles[0].x[1] = HLY;
	particles[0].x[2] = HLZ;

	particles[1].x[0] = HLX + dmy;
	particles[1].x[1] = HLY;
	particles[1].x[2] = HLZ;
	//particles[0].v[0] = 10.0;
	//particles[0].v[1] = 1.0;

	for(int n=0;n< Particle_Number;n++){
	  for(int d=0;d< DIM; d++){
	    particles[n].v[d] = 0.;
	    particles[n].omega[d] = 0.;
	    particles[n].fr[d] = 0.;
	    particles[n].fv[d] = 0.;
	    particles[n].torquer[d] = 0.;
	    particles[n].torquev[d] = 0.;
	  }
	}
      }
    }
    */
    if(Lubrication_measurement 
       + Terminal_velocity_measurement 
       + DKT_measurement
       + Check_PB
       > 1){
      fprintf(stderr, "invalid compile-switch.\n");
      fprintf(stderr, "CHECK\n");
      fprintf(stderr, "\tLubrication_measurement,\n");
      fprintf(stderr, "\tTerminal_velocity_measurement,\n");
      fprintf(stderr, "\tDKT_measurement,\n");
      fprintf(stderr, "\tCheck_PB,\n");
      fprintf(stderr, "in input.cxx\n");
      exit_job(EXIT_FAILURE);
    }else {
      if(Lubrication_measurement){
	Lubrication(INIT, stderr, particles,jikan);
      }else if(Terminal_velocity_measurement){
	Terminal_velocity(INIT, stderr, particles,jikan);
      }else if(DKT_measurement){
	if(SW_EQ != Electrolyte){
	  DKT(INIT, stderr, particles,jikan);
	}else {
	  DKT(INIT, stdout, particles,jikan);
	}
      }else if(Check_PB){
	Check_Poisson_Boltzmann(INIT, stderr, particles, Concentration
				,up[0], up[1], up[2]);
      }
    }
  }
  if(0){
    fprintf(stdout,"%g %g %g %g %g %g\n"
	    ,particles[0].v[0]
	    ,particles[0].v[1]
	    ,particles[0].v[2]
	    ,particles[0].omega[0]
	    ,particles[0].omega[1]
	    ,particles[0].omega[2]);
  }
  Init_zeta_k(zeta, uk_dc);
  if(SW_EQ == Qian_Sheng 
     || SW_EQ == nematic_Allen_Cahn 
     || SW_EQ == Olmsted_Goldbart
     ){
    Init_tensor_order(Tensor_order, particles);
  }else if(SW_EQ == Two_fluid ){
    Init_solute(Concentration[0], particles);
  }
  {
    Reset_phi_u(phi,up);
    Make_phi_u_particle(phi, up, particles);
    Zeta_k2u(zeta, uk_dc, u);
    Make_f_particle_dt_sole(f_particle, u, up, phi);
    Add_f_particle(u, f_particle);
    U2zeta_k(zeta, uk_dc, u);
    if(1){
      for(int d=0;d<DIM;d++){
	uk_dc[d] = .0;
      }
    }
  }
  
  if(SW_EQ==Electrolyte){
    Init_rho_ion(Concentration, particles, jikan);
  }
  //  return EXIT_SUCCESS;
  Show_parameter(Avs_parameters, particles);
  if(SW_EQ == Navier_Stokes){
    //if(kBT > 0){
    double base_v[DIM] = {0.,0.,0.};
    Count_noisy_NS(INIT,stderr,particles,uk_dc,zeta, phi, u, base_v);
    Count_noisy_particle(INIT, stderr, base_v, particles, jikan);
  }else if(SW_EQ == Slippy_Navier_Stokes){
    Count_Slippy_NS(particles, INIT);
  }else if(SW_EQ == Shear_Navier_Stokes){
    Mean_shear_stress(INIT, stderr, NULL, particles, jikan, zeta, uk_dc, u);
  }else if(SW_EQ == Electrolyte){
    if(Lubrication_measurement
       || Terminal_velocity_measurement
       ){
      Electrolyte_free_energy(INIT,stdout,particles,Concentration,jikan);
    }else {
      Electrolyte_free_energy(INIT,stderr,particles,Concentration,jikan);
    }
  }
  if(Check_PB){
    Check_Poisson_Boltzmann(SHOW, stderr, particles, Concentration
			    ,up[0], up[1], up[2]);
  }
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  int resumed_ts = 0;
  if(RESUMED) resumed_ts = last_ts;
  for(jikan.ts = resumed_ts; jikan.ts <= MSTEP; jikan.ts++){
    int resumed_and_1st_loop = 0;
    if(RESUMED && (jikan.ts == resumed_ts)){
      resumed_and_1st_loop = 1;
    }
    if( jikan.ts % GTS == 0){
      //Show_particle(particles);
      if(!resumed_and_1st_loop){
	if(SW_AVS){
	  if(Particle_Number > 0){
	    Output_avs_p(Avs_parameters, particles, jikan);
	  }
	  if(SW_EQ == Navier_Stokes 
	     || SW_EQ == Slippy_Navier_Stokes 
	     || SW_EQ == Shear_Navier_Stokes
	     ){
	    Output_avs(Avs_parameters, zeta, uk_dc, particles, jikan);
	  }else if(SW_EQ == Two_fluid){
	    Output_avs_two_fluid(Avs_parameters, zeta, uk_dc, u
				 ,Concentration, Concentration_rhs0
				 ,particles, phi, jikan);
	  }else if(SW_EQ==Qian_Sheng
		   || SW_EQ == nematic_Allen_Cahn
		   || SW_EQ == Olmsted_Goldbart
		   ){
	    Output_avs_QS(Avs_parameters, zeta, uk_dc, Tensor_order, particles, jikan);
	  }else if(SW_EQ==Electrolyte){
	    Output_avs_charge(Avs_parameters, zeta, uk_dc, Concentration, particles, jikan);
	  }
	}
	if(SW_UDF){
	  Output_udf(ufout, Avs_parameters, zeta, uk_dc, particles, jikan);
	}
	if(SW_EQ == Electrolyte){
	  if(Lubrication_measurement
	     || Terminal_velocity_measurement
	     ){
	    Electrolyte_free_energy(SHOW,stdout,particles, Concentration,jikan);
	  }else {
	    Electrolyte_free_energy(SHOW,stderr,particles, Concentration,jikan);
	  }
	}
      }
    }

    if(0){
    //if(kBT <= 0.0){
      double dmy = Particle_Reynolds_number(particles);
      fprintf(stderr, "%g\n", dmy);
    } 
    if(Lubrication_measurement){
      Lubrication(SHOW,stderr,particles,jikan);
    }else if(Terminal_velocity_measurement){
      if(jikan.ts > MSTEP/2){
	Terminal_velocity(ADD,stderr,particles,jikan);
      }else {
	Terminal_velocity(SHOW,stderr,particles,jikan);
      }
    }else if(DKT_measurement){
      if(SW_EQ != Electrolyte){
	DKT(SHOW,stderr,particles,jikan);
      }else {
	DKT(SHOW,stdout,particles,jikan);
      }
    }
    {
      {
	if(SW_EQ == Two_fluid ){
	  if(0){
	    fprintf(stderr,"%g %g %g %g %g %g\n"
		    ,particles[0].v[0]
		    ,particles[0].v[1]
		    ,particles[0].v[2]
		    ,particles[0].omega[0]
		    ,particles[0].omega[1]
		    ,particles[0].omega[2]);
	  }
	}else if(SW_EQ == Navier_Stokes){
	  //if(kBT <= 0.){
	  //{
	  if(jikan.ts > MSTEP/2){
	    {
	      double base_v[DIM] = {0.,0.,0.};
	      Count_noisy_NS(ADD,stderr,particles,uk_dc,zeta, phi, u, base_v);    
	      if(1){
		//const double ivolume = 1./(double)(NX*NY*NZ);
		const double ivolume = Ivolume * POW3(DX);
		for(int d=0;d<DIM;d++){
		  base_v[d] = uk_dc[d] * ivolume;
		  //base_v[d] = 0.;
		}
	      }
	      if(0){
		//const double ivolume = 1./(double)(NX*NY*NZ);
		const double ivolume = Ivolume * POW3(DX);
		fprintf(stdout, "%g %g %g %g %g %g #bbb:\n"
			,base_v[0]
			,base_v[1]
			,base_v[2]
			,uk_dc[0]*ivolume
			,uk_dc[1]*ivolume
			,uk_dc[2]*ivolume
			);
	      }
	      if(Particle_Number > 0){
		Count_noisy_particle(ADD, stderr, base_v, particles, jikan);
	      }
	    }
	    //if(0){//////////////////
	    if(kBT > 0){//////////////////
	      if(Particle_Number > 0){
		Count_noisy_particle(SNAP_MEAN, stderr, NULL, particles, jikan);
	      }
	    }//////////////////
	  }
	}else if(SW_EQ == Slippy_Navier_Stokes){
	  if(jikan.ts > MSTEP/2){
	    Count_Slippy_NS(particles, ADD);
	  }
	}
      }
    }
    Time_evolution(zeta, uk_dc, f_particle, particles, jikan);
    jikan.time += jikan.dt_fluid;

    if(SW_EQ == Shear_Navier_Stokes){
      if(jikan.ts > MSTEP/4){
	Mean_shear_stress(ADD, stderr, dev_shear_stress, particles, jikan, zeta, uk_dc, u);
      }else {
	Mean_shear_stress(SHOW, stderr, dev_shear_stress, particles, jikan, zeta, uk_dc, u);
      }
    }
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //-------------------------
    if(jikan.ts == MSTEP){
      
      Save_Restart_udf(ufout
		       ,zeta
		       ,uk_dc
		       ,particles
		       ,jikan
		       ,Concentration
		       );
    }

    //-------------------------
    if(resumed_and_1st_loop){
      Force_restore_parameters(ufout
			       ,zeta
			       ,uk_dc
			       ,particles
			       ,jikan
			       ,Concentration
			       );
      delete ufin;
      fprintf(stderr,"############################ Parameters are restored.\n");
    }
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
  }

  if(SW_EQ == Shear_Navier_Stokes){
    Mean_shear_stress(MEAN, stderr, NULL, particles, jikan, zeta, uk_dc, u);
  }
  if(Terminal_velocity_measurement){
    Terminal_velocity(MEAN,stderr,particles,jikan);
  }
  if(!Lubrication_measurement 
     && SW_EQ == Navier_Stokes 
     && kBT > 0.
     ){
    if(Particle_Number > 0){
      Count_noisy_particle(MEAN, stderr, NULL, particles, jikan);
    }
    {
      double base_v[DIM] = {0.,0.,0.};
      Count_noisy_NS(MEAN,stderr, particles,uk_dc,zeta, phi, u,base_v);
    }
  }
  if(SW_EQ == Slippy_Navier_Stokes){
    Count_Slippy_NS(particles, MEAN);
  }
  if(SW_EQ == Electrolyte){
    if(Lubrication_measurement
       || Terminal_velocity_measurement
       ){
      Electrolyte_free_energy(MEAN,stdout,particles,Concentration,jikan);
    }else {
      Electrolyte_free_energy(MEAN,stderr,particles,Concentration,jikan);
    }
  }

  if(SW_UDF){
    ufout->write();
    delete ufout;
    fprintf(stderr,"#%s end.\n",Out_udf);
    ufres->write();
    delete ufres;
    fprintf(stderr,"#%s end.\n",Res_udf);
  }
  delete [] particles;
  return EXIT_SUCCESS;
}
