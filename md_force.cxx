//
// $Id: md_force.cxx,v 1.30 2005/09/16 05:05:58 nakayama Exp $
//

#include "md_force.h"

double Min_rij;
double Max_force;
double Min_rij_wall;
double Max_force_wall;

double Calc_f_Lennard_Jones_shear_cap_primitive(Particle *p
				       ,void (*distance0_func)(const double *x1,const double *x2,double &r12,double *x12)
				      ,const double cap
				       ){
  Min_rij = DBL_MAX;
  Max_force = 0.;
  int SW_minrij = 1;
  // Particle 変数の f に 
  // !! += 
  //で足す. f の初期値 が正しいと仮定している!!
  const double LJ_cutoff = A_R_cutoff * LJ_dia;

  double shear_stress=0.0;
  for(int n=0;n<Particle_Number ; n++){
    Particle *p_n = &p[n];
    for(int m=n+1; m < Particle_Number ; m++){
      double r_ij_vec[DIM];
      double r_ij;
      distance0_func( (*p_n).x, p[m].x, r_ij, r_ij_vec);
      if(SW_minrij){
	Min_rij = MIN(r_ij, Min_rij);
      }
      if(r_ij < LJ_cutoff){
	double dmy = MIN(cap/r_ij,Lennard_Jones_f( r_ij , LJ_dia));
	double dmyf = 0.;
	for(int d=0; d < DIM; d++ ){ 
	  double dmy1 = dmy * -r_ij_vec[d];
	  if(SW_minrij){
	    dmyf += SQ(dmy1);
	  }
	  (*p_n).fr[d] += dmy1;
	  p[m].fr[d] -= dmy1;
	}
	if(SW_minrij){
	  Max_force = MAX(sqrt(dmyf),Max_force);
	}
	shear_stress += ((dmy * -r_ij_vec[0])* (-r_ij_vec[1]));
      }
    }
  }
  return shear_stress;
}

void Add_f_gravity(Particle *p){
  static const double Gravity_on_fluid 
    = G*RHO * 4./3.*M_PI * SQ(RADIUS)* RADIUS;
   // Particle 変数の f に 
  // !! += 
  //で足す. f の初期値 が正しいと仮定している!!
  for(int n = 0; n < Particle_Number ; n++){
    p[n].fr[G_direction] -= Gravity_on_fluid * (MASS_RATIOS[p[n].spec] -1.0); 
  }
}

void Calc_f_hydro_draining(Particle *p, const Value *fp, const CTime &jikan){
  static const double Zeta_drag = 6.* M_PI * ETA * RADIUS;
  for(int n = 0; n < Particle_Number ; n++){
    for(int d=0; d < DIM; d++ ){ 
      p[n].f_hydro[d] = -Zeta_drag * p[n].v[d];
      //      printf("%f\n",Zeta_drag);
    }
  }
}
void Calc_f_hydro_draining_shear(Particle *p, const Value *fp, const CTime &jikan){
  static const double Zeta_drag = 6.* M_PI * ETA * RADIUS;
  const double hly_shear = LY_shear*0.5;
  for(int n = 0; n < Particle_Number ; n++){
    p[n].f_hydro[0] = -Zeta_drag
      * (p[n].v[0]-Shear_rate*(hly_shear - p[n].x[1]));
    p[n].f_hydro[1] = -Zeta_drag * p[n].v[1];
    p[n].f_hydro[2] = -Zeta_drag * p[n].v[2];
  }
}

double Collison_time(Particle *p
		     ,void (*distance0_func)(const double *x1,const double *x2,double &r12,double *x12)
		     ,const Particle_BC sw_bc
		     ){
  const double shear_max = Shear_rate * LY_shear;
  
  double t_collision = DBL_MAX;

  for(int n=0;n<Particle_Number ; n++){
    Particle *p_n = &p[n];
    for(int m=n+1; m < Particle_Number ; m++){
      double r_ij_vec[DIM];
      double r_ij;
      distance0_func( (*p_n).x, p[m].x, r_ij, r_ij_vec);
      {
	double ir_ij2 = SQ(1./r_ij);
	double amplitude=0.0;
	double v_diff[DIM]={0.,0.,0.};
	for(int d=0; d < DIM; d++ ){ 
	  v_diff[d] = p[m].v[d] - (*p_n).v[d];
	  if(d == 0 && sw_bc == Lees_Edwards ){
	    int cell1 = Nint(( p[m].x[1]-(*p_n).x[1] )/L_particle[1] );
	    v_diff[d] += cell1 * shear_max;
	  }
	  amplitude += (v_diff[d] * r_ij_vec[d]);
	}
	if(amplitude > 0.){
	  continue;
	}else {  // approaching
	  double v_sq = 0.;
	  for(int d=0; d < DIM; d++ ){ 
	    v_sq += SQ(v_diff[d]);
	  }
	  double dmy = SQ(amplitude)- v_sq * (SQ(r_ij)-SQ(SIGMA));
	  if(dmy < 0.){
	    continue;
	  }else {
	    t_collision = MIN(t_collision, ((-amplitude - sqrt(dmy))*ir_ij2));
	  }
	}
      }
    }
  }
  return t_collision;
}
double Calc_f_hydro_lubrication_shear_primitive(Particle *p
						,void (*distance0_func)(const double *x1,const double *x2,double &r12,double *x12)
						,const Particle_BC sw_bc
						){
  static const double lub_cutoff = 3.*RADIUS;//HLX;//2.*RADIUS+ SIGMA;
  static const double Zeta_lub = 6.* M_PI * ETA * SQ(RADIUS*.5);
  const double shear_max = Shear_rate * LY_shear;

  double shear_stress=0.0;
  for(int n=0;n<Particle_Number ; n++){
    Particle *p_n = &p[n];
    for(int m=n+1; m < Particle_Number ; m++){
      double r_ij_vec[DIM];
      double r_ij;
      distance0_func( (*p_n).x, p[m].x, r_ij, r_ij_vec);
      double gap;
      if(sw_bc == Lees_Edwards){
	gap = r_ij;// -SIGMA -XI;
	//gap = r_ij -SIGMA;
      }else {
	gap = r_ij - SIGMA;
      }
      if(r_ij < lub_cutoff && gap > 0.0){
	double ir_ij2 = SQ(1./r_ij);
	double igap = 1./gap;
	double amplitude=0.0;
	double v_diff[DIM]={0.,0.,0.};
	for(int d=0; d < DIM; d++ ){ 
	  v_diff[d] = p[m].v[d] - (*p_n).v[d];
	  amplitude += (v_diff[d] * r_ij_vec[d]);
	}
	if(sw_bc == Lees_Edwards){
	  int cell1 = Nint(( p[m].x[1]-(*p_n).x[1] )/L_particle[1] );
	  v_diff[0] += cell1 * shear_max;
	}
	amplitude *= (Zeta_lub * igap * ir_ij2);
	if(sw_bc != Shear_hydro){
	  for(int d=0; d < DIM; d++ ){ 
	    double dmy = amplitude * r_ij_vec[d];
	    (*p_n).fv[d] += dmy;
	    p[m].fv[d] -= dmy;
	  }
	}
	if(sw_bc != PBC_particle){
	  shear_stress += ((amplitude*r_ij_vec[0])* (-r_ij_vec[1]));
	}
      }
    }
  }
  return shear_stress;
}
inline double weight(const double &r, const double &cutoff){
  return r<cutoff?(cutoff/r-1.):0.;
}
void Calc_f_contanct_nonslip(Particle *p){
  static const double nonslip_cutoff = 1.2 * SIGMA * pow(2.,1./6.);
  static const double spring_cst = 2.5e-1;//6.* M_PI * ETA * RADIUS;
  const int NOrot = 0;

  for(int n=0;n<Particle_Number ; n++){
    Particle *p_n = &p[n];
    for(int m=n+1; m < Particle_Number ; m++){
      double r_ij_vec[DIM];
      double r_ij;
      Distance0( (*p_n).x, p[m].x, r_ij, r_ij_vec);
      if(r_ij < nonslip_cutoff){
	double ir_ij = 1./r_ij;
	double ir_ij2 = SQ(ir_ij);

	double dv_longitudinal = 0.0;
	double v_diff[DIM]={0.,0.,0.};
	double omega_diff[DIM]={0.,0.,0.};

	for(int d=0; d < DIM; d++ ){ 
	  v_diff[d] = p[m].v[d] - (*p_n).v[d];
	  omega_diff[d] = p[m].omega[d] - (*p_n).omega[d];
	  dv_longitudinal += (v_diff[d] * r_ij_vec[d]);
	}

	dv_longitudinal *= ir_ij2;

	double dmy_omega[DIM];
	dmy_omega[0] = 
	  omega_diff[1] * r_ij_vec[2] 
	  - omega_diff[2] * r_ij_vec[1];
	dmy_omega[1] = 
	  omega_diff[2] * r_ij_vec[0] 
	  - omega_diff[0] * r_ij_vec[2];
	dmy_omega[2] = 
	  omega_diff[0] * r_ij_vec[1] 
	  - omega_diff[1] * r_ij_vec[0];

	double force[DIM];
	for(int d=0; d < DIM; d++ ){ 
	  if(NOrot){
	    dmy_omega[d] = 0.;
	  }else{
	    dmy_omega[d] *= (RADIUS*ir_ij);
	  }
	  //	  fprintf(stderr, "aaa:%g %g %g %g\n"
	  //  ,spring_cst
	  //  ,spring_cst * weight(r_ij, nonslip_cutoff)
	  //  ,r_ij, nonslip_cutoff);
	  force[d] = spring_cst * weight(r_ij, nonslip_cutoff)
	    * (v_diff[d] - dv_longitudinal * r_ij_vec[d] + dmy_omega[d]);
	}
	double torque[DIM];
	{
	  torque[0] = 
	    r_ij_vec[1] * force[2]
	    - r_ij_vec[2] * force[1];
	  torque[1] = 
	    r_ij_vec[2] * force[0]
	    - r_ij_vec[0] * force[2];
	  torque[2] = 
	    r_ij_vec[0] * force[1]
	    - r_ij_vec[1] * force[0];
	}
	
	for(int d=0; d < DIM; d++ ){ 
	  (*p_n).fv[d] += force[d];
	  p[m].fv[d] -= force[d];
	  
	  torque[d] *= (RADIUS * ir_ij); 
	  (*p_n).torquev[d] += torque[d];
	  p[m].torquev[d] += torque[d];

	}
	//fprintf(stderr, "%g %g %g\n", r_ij, gap, amplitude);
      }
    }
  }
}


inline double Lennard_Jones_rep_f(const double &x, const double sigma){
  //    printf("%d\n",LJ_powers);
  // ! x== 0.0 の処理を省略
  double answer=0.0;
  if(x < A_R_cutoff * sigma ){
    int factor = LJ_powers + 1;
    // LJ_powers
    // 0: (12:6)
    // 1: (24:12)
    // 2: (36:18)
    static const double LJ_coeff1= 24. * (double)factor * EPSILON;
    double dmy0 = sigma/x;
    dmy0 =  SQ(dmy0) * SQ(dmy0) * SQ(dmy0);
    double dmy = 1.0;
    for(int n = 0;n< factor ;n++){
      dmy *= dmy0;
    }
    answer = LJ_coeff1 / SQ(x) * ( 2.0 * SQ(dmy));
  }else {
    answer = 0.0;
  }
  return answer;
}

void Calc_f_Coulomb_friction(Particle *p){
  static const double friction_cutoff = 1.2 * SIGMA * pow(2.,1./6.);
  static const double friction_cst = 1.e-1;//6.* M_PI * ETA * RADIUS;
  static const double iv_transition = DT/DX;//6.* M_PI * ETA * RADIUS;
  //fprintf(stderr, "%g %g\n", 1/iv_transition, p[0].v[0]);
  for(int n=0;n<Particle_Number ; n++){
    Particle *p_n = &p[n];
    for(int m=n+1; m < Particle_Number ; m++){
      double r_ij_vec[DIM];
      double r_ij;
      Distance0( (*p_n).x, p[m].x, r_ij, r_ij_vec);
      if(r_ij < friction_cutoff){
	double ir_ij = 1./r_ij;
	double ir_ij2 = SQ(ir_ij);

	double dv_longitudinal = 0.0;
	double v_diff[DIM]={0.,0.,0.};
	double omega_diff[DIM]={0.,0.,0.};

	for(int d=0; d < DIM; d++ ){ 
	  v_diff[d] = p[m].v[d] - (*p_n).v[d];
	  omega_diff[d] = p[m].omega[d] - (*p_n).omega[d];
	  dv_longitudinal += (v_diff[d] * r_ij_vec[d]);
	}

	dv_longitudinal *= ir_ij2;

	double dmy_omega[DIM];
	dmy_omega[0] = 
	  omega_diff[1] * r_ij_vec[2] 
	  - omega_diff[2] * r_ij_vec[1];
	dmy_omega[1] = 
	  omega_diff[2] * r_ij_vec[0] 
	  - omega_diff[0] * r_ij_vec[2];
	dmy_omega[2] = 
	  omega_diff[0] * r_ij_vec[1] 
	  - omega_diff[1] * r_ij_vec[0];

	double force[DIM];
	double normal_force = friction_cst 
	  * Lennard_Jones_rep_f( r_ij, LJ_dia);
	for(int d=0; d < DIM; d++ ){ 
	  dmy_omega[d] *= (RADIUS*ir_ij);
	  double dvt= (v_diff[d] - dv_longitudinal * r_ij_vec[d] + dmy_omega[d]);
	  double sgn = SGN(dvt);
	  //double amp =  -expm1(iv_transition * fabs(dvt));
	  double amp =  1.-exp(iv_transition * fabs(dvt));
	  force[d] = sgn * normal_force;
	}
	double torque[DIM];
	{
	  torque[0] = 
	    r_ij_vec[1] * force[2]
	    - r_ij_vec[2] * force[1];
	  torque[1] = 
	    r_ij_vec[2] * force[0]
	    - r_ij_vec[0] * force[2];
	  torque[2] = 
	    r_ij_vec[0] * force[1]
	    - r_ij_vec[1] * force[0];
	}
	
	for(int d=0; d < DIM; d++ ){ 
	  (*p_n).fr[d] += force[d];
	  p[m].fr[d] -= force[d];
	  
	  //torque[d] *= (RADIUS * ir_ij); 
	  torque[d] *= 0.0;
	  (*p_n).torquer[d] += torque[d];
	  p[m].torquer[d] += torque[d];

	}
	//fprintf(stderr, "%g %g %g\n", r_ij, gap, amplitude);
      }
    }
  }
}

void Calc_f_hydro_correct(Particle *p, const Value *fp, const CTime &jikan){
  static const double dmy0 = -DX3*RHO;
  double dmy = dmy0 /jikan.dt_fluid; 

  int *nlattice;
  if(SW_EQ == Shear_Navier_Stokes){
    nlattice = Ns_shear;
  }else {
    nlattice = Ns;
  } 

  for(int n = 0; n < Particle_Number ; n++){
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell 
      = Particle_cell(p[n].x, DX, x_int, residue);// {1,0} が返ってくる
    sw_in_cell = 1;
    double force[DIM];
    double torque[DIM];
    for(int d=0; d < DIM; d++){
      force[d] = 0.0;
      torque[d] = 0.0;
    }
    int r_mesh[DIM];
    double r[DIM];
    for(int mesh=0; mesh < NP_domain; mesh++){
      Relative_coord(Sekibun_cell[mesh], x_int, residue, sw_in_cell, nlattice, DX, r_mesh, r);
      double dmy_fp[DIM];
      for(int d=0; d < DIM; d++ ){ 
	dmy_fp[d] = fp[d][r_mesh[0]][r_mesh[1]][r_mesh[2]];
	force[d] += dmy_fp[d];
      }
      {// torque
	torque[0] += (r[1] * dmy_fp[2] - r[2] * dmy_fp[1]);
	torque[1] += (r[2] * dmy_fp[0] - r[0] * dmy_fp[2]);
	torque[2] += (r[0] * dmy_fp[1] - r[1] * dmy_fp[0]);
      }
    }
    for(int d=0; d < DIM; d++ ){ 
      p[n].f_hydro[d] = (dmy * force[d] * p[n].eff_mass_ratio);
    }
    if(ROTATION){
      for(int d=0; d < DIM; d++ ){ 
	p[n].torque_hydro[d] = ( dmy * torque[d] * p[n].eff_mass_ratio);
      }
    }
  }
}

////////////////////
void Calc_particle_velocity(const Value &phi, const Value *u, Particle *p){
  for(int n = 0; n < Particle_Number ; n++){
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell 
      = Particle_cell(p[n].x, DX, x_int, residue);// {1,0} が返ってくる
    sw_in_cell = 1;
    double volume_phi=0.0;
    double velocity[DIM];
    for(int d=0; d < DIM; d++){
      velocity[d] = 0.0;
    }
    int r_mesh[DIM];
    double r[DIM];
    for(int mesh=0; mesh < NP_domain; mesh++){
      Relative_coord(Sekibun_cell[mesh], x_int, residue, sw_in_cell, Ns, DX, r_mesh, r);
      double dmy_phi = phi[r_mesh[0]][r_mesh[1]][r_mesh[2]];
      volume_phi +=  dmy_phi;
      for(int d=0; d < DIM; d++ ){ 
	velocity[d] += ( u[d][r_mesh[0]][r_mesh[1]][r_mesh[2]] * dmy_phi );
      }
    }
    double ivolume_phi = 1./volume_phi;
    for(int d=0; d < DIM; d++ ){ 
      p[n].v[d] = (velocity[d] * ivolume_phi);
    }
  }
}
void Calc_f_Lennard_Jones_wall_cap(Particle *p
				   ,const int &dwall
				   ,const double &radius_wall
				   ,const double cap
				   ){
  Min_rij_wall = DBL_MAX;
  Max_force_wall = 0.;
  int SW_minrij = 1;
  // Particle 変数の f に 
  // !! += 
  //で足す. f の初期値 が正しいと仮定している!!
  //const double LJ_cutoff = LJ_dia;
  const double LJ_cutoff = pow(2.0,1./6.) * LJ_dia;
  const double wall_depth = LJ_dia*.5 - radius_wall;

  for(int n=0;n<Particle_Number ; n++){
    Particle *p_n = &p[n];
    double wall_x[DIM];
    {// lower wall
      {
	for(int d=0;d<DIM;d++){
	  wall_x[d] = (*p_n).x[d];
	}
	wall_x[dwall] = -wall_depth;
      }
      double r_ij = (*p_n).x[dwall] - wall_x[dwall];
      if(SW_minrij){
	Min_rij_wall = MIN(r_ij, Min_rij_wall);
      }
      if(r_ij < LJ_cutoff){
	//fprintf(stderr, "aaa:%g %g\n",r_ij, LJ_cutoff);
	double r_ij_vec[DIM]={0.,0.,0.};
	{
	  r_ij_vec[dwall] = -r_ij;
	}

	double dmy = MIN(cap/r_ij,Lennard_Jones_f( r_ij , LJ_dia));
	if(SW_minrij){
	  Max_force_wall = MAX(dmy,Max_force_wall);
	}
	{
	  double dmy1 = dmy * -r_ij_vec[dwall];
	  (*p_n).fr[dwall] += dmy1;
	}
      }
    }
    {// upper wall
      {
	for(int d=0;d<DIM;d++){
	  wall_x[d] = (*p_n).x[d];
	}
	wall_x[dwall] = L_particle[dwall] + wall_depth;
      }
      double r_ij = wall_x[dwall] - (*p_n).x[dwall];
      if(SW_minrij){
	Min_rij_wall = MIN(r_ij, Min_rij_wall);
      }
      if(r_ij < LJ_cutoff){
	double r_ij_vec[DIM]={0.,0.,0.};
	{
	  r_ij_vec[dwall] = r_ij;
	}
	double dmy = MIN(cap/r_ij,Lennard_Jones_f( r_ij , LJ_dia));
	if(SW_minrij){
	  Max_force_wall = MAX(dmy,Max_force_wall);
	}
	{
	  double dmy1 = dmy * -r_ij_vec[dwall];
	  (*p_n).fr[dwall] += dmy1;
	}
      }
    }
  }
}
