/*!
  \file md_force.cxx
  \brief Routines to compute MD forces on particles
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
 */

#include "md_force.h"

double Min_rij;
double Max_force;
double *Hydro_force;

double *Hydro_force_new;
double *Hydro_force_new_u;
double *Hydro_force_new_p;
double *Hydro_force_new_v;
double *Hydro_force_new_w;

#define Cell_length 16

double Calc_f_Lennard_Jones_shear_cap_primitive_lnk(Particle *p
				       ,void (*distance0_func)(const double *x1,const double *x2,double &r12,double *x12)
				      ,const double cap
				       ){
  Min_rij = DBL_MAX;
  Max_force = 0.;
  int SW_minrij = 1;
  // Particle ïœêîÇÃ f Ç… 
  // !! += 
  //Ç≈ë´Ç∑. f ÇÃèâä˙íl Ç™ê≥ÇµÇ¢Ç∆âºíËÇµÇƒÇ¢ÇÈ!!
  const double LJ_cutoff = A_R_cutoff * LJ_dia;
  double r_ij_vec[DIM];
  double r_ij;
  int rigid_pair;
  double shear_stress=0.0;

// List Constructor
  int i,j,k;
  int *lscl;
  lscl = alloc_1d_int(Particle_Number);
  int lc[DIM];
  double lc_r[DIM];
  int mc[DIM];
  int lcyz,lcxyz;
  int ic[DIM],cn;
  lc[0] = int (NX / Cell_length);
  lc[1] = int (NY / Cell_length);
  lc[2] = int (NZ / Cell_length);
  for(int d=0; d < DIM; d++ ){ 
      lc_r[d] = Cell_length * DX;
  }
  lcyz = lc[1]*lc[2];
  lcxyz = lc[0]*lcyz;
  
  int *head;
  head = alloc_1d_int(lcxyz);
  
#pragma omp parallel for schedule(dynamic, 1)
  for(int d=0; d < lcxyz; d++) head[d] = -1;
  for(int n=0;n<Particle_Number ; n++){
      Particle *p_n = &p[n];
      for(int d=0; d < DIM; d++ ){ 
	  mc[d] = int( (*p_n).x[d] / lc_r[d] ) ;
      }
      cn = mc[0]*lcyz+mc[1]*lc[2]+mc[2];
      lscl[n] = head[cn];
      head[cn] = n;
  }
  
  for (ic[0]= 0 ; ic[0] < lc[0] ; ic[0]++){ 
      for (ic[1]= 0 ; ic[1] < lc[1] ; ic[1]++){
	  for (ic[2]= 0 ; ic[2] < lc[2] ; ic[2]++){
	      cn = ic[0]*lcyz+ic[1]*lc[2]+ic[2];
// Scan the neighbor cells (including itself) of cell c
	      int icl[DIM],icls[DIM];
	      for (icl[0]= ic[0]-1 ; icl[0] <= ic[0]+1 ; icl[0]++){ 
		  for (icl[1]= ic[1]-1 ; icl[1] <= ic[1]+1 ; icl[1]++){
		      for (icl[2]= ic[2]-1 ; icl[2] <= ic[2]+1 ; icl[2]++){
// Consider periodic condition
			  for(int d=0; d < DIM; d++ ){ 
			      icls[d] = icl[d];
			      if(icls[d] < 0){
				  icls[d] = lc[d] -1;
			      }
			      if(icls[d] > (lc[d]-1)){
				  icls[d] = 0;
			      }
			  }
// Calculate the scalar cell index of the neighbor cell
			  int cl = ((icls[0]+lc[0])%lc[0])*lcyz+((icls[1]+lc[1])%lc[1])*lc[2]+((icls[2]+lc[2])%lc[2]);
// Scan atom i in cell c
			  i = head[cn];
			  while (i != -1){
			      j = head[cl];
			      while (j != -1){
                                rigid_pair = rigid_chain(i, j);
				  if (i > j && !rigid_pair) {
				      distance0_func( p[i].x, p[j].x, r_ij, r_ij_vec);
				      if(SW_minrij){
					  Min_rij = MIN(r_ij, Min_rij);
				      }
				      if (r_ij < LJ_cutoff){
					  double dmy = MIN(cap/r_ij,Lennard_Jones_f( r_ij , LJ_dia));
					  double dmyf = 0.;
					  for(int d=0; d < DIM; d++ ){ 
					      double dmy1 = dmy * -r_ij_vec[d];
					      if(SW_minrij){
						  dmyf += SQ(dmy1);
					      }
					      p[i].fr[d] += dmy1;
					      p[j].fr[d] -= dmy1;
					  }
					  if(SW_minrij){
					      Max_force = MAX(sqrt(dmyf),Max_force);
					  }
					  shear_stress += ((dmy * -r_ij_vec[0])* (r_ij_vec[1]));
				      }
				  }
				  j = lscl[j];
			      }
			      i = lscl[i];
			  }
		      }
		  }
	      }
	  }
      }
  }
  free_1d_int(lscl);
  free_1d_int(head); 
  return shear_stress;
}

double Calc_f_Lennard_Jones_shear_cap_primitive(Particle *p
						,void (*distance0_func)(const double *x1,const double *x2,double &r12,double *x12)
						,const double cap
				       ){
    Min_rij = DBL_MAX;
    Max_force = 0.;
    int SW_minrij = 1;
    // Particle $BJQ?t$N(B f $B$K(B 
    // !! += 
    //$B$GB-$9(B. f $B$N=i4|CM(B $B$,@5$7$$$H2>Dj$7$F$$$k(B!!
    const double LJ_cutoff = A_R_cutoff * LJ_dia;
    double shear_stress=0.0;
    int rigid_pair;
    for(int n=0;n<Particle_Number ; n++){
	Particle *p_n = &p[n];
	for(int m=n+1; m < Particle_Number ; m++){
	    double r_ij_vec[DIM];
	    double r_ij;
	    distance0_func( (*p_n).x, p[m].x, r_ij, r_ij_vec);
	    if(SW_minrij){
		Min_rij = MIN(r_ij, Min_rij);
	    }
            rigid_pair = rigid_chain(n, m);
	    if(r_ij < LJ_cutoff && !rigid_pair){
		
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
		shear_stress += ((dmy * -r_ij_vec[0])* (r_ij_vec[1]));
	    }
	}
    }
    return shear_stress;
}

void Add_f_gravity(Particle *p){
  static const double Gravity_on_fluid 
    = G*RHO * 4./3.*M_PI * SQ(RADIUS)* RADIUS;
   // Particle $BJQ?t$N(B f $B$K(B 
  // !! += 
  //$B$GB-$9(B. f $B$N=i4|CM(B $B$,@5$7$$$H2>Dj$7$F$$$k(B!!
#pragma omp parallel for schedule(dynamic, 1) 
  for(int n = 0; n < Particle_Number ; n++){
    p[n].fr[G_direction] -= Gravity_on_fluid * (MASS_RATIOS[p[n].spec] -1.0); 
  }
}
void Calc_f_slip_correct_precision(Particle *p, double const* const* u, const CTime &jikan){
    static const double dmy0 = DX3*RHO; 
    double dmy = dmy0 /jikan.dt_fluid; 
    int *nlattice;
    nlattice = Ns;

    double xp[DIM];
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell;
    double force[DIM];
    double torque[DIM];
    int r_mesh[DIM];
    double r[DIM];
    double dmy_fp[DIM];
    double x[DIM];
    double dmyR;
    double dmy_phi;
    double dmy_mass;
    int pspec;
#pragma omp parallel for schedule(dynamic, 1) \
  private(xp,x_int,residue,sw_in_cell,force,torque,r_mesh,r,dmy_fp,x,dmyR,dmy_phi,dmy_mass,pspec) 
    for(int n = 0; n < Particle_Number; n++){
      for (int d = 0; d < DIM; d++) {
	xp[d] = p[n].x[d];

        force[d] = torque[d] = 0.0;
      }
      
      sw_in_cell = Particle_cell(xp, DX, x_int, residue);
      sw_in_cell = 1;
      
      for(int mesh=0; mesh < NP_domain; mesh++){
	Relative_coord(Sekibun_cell[mesh], x_int, residue, sw_in_cell, nlattice, DX, r_mesh, r);
	double x[DIM];
	for(int d=0;d<DIM;d++){
	  x[d] = r_mesh[d] * DX;
	}
	dmyR = Distance(x, xp);
	dmy_phi= Phi(dmyR, RADIUS);
	
	int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
	for(int d=0; d < DIM; d++ ){
	  dmy_fp[d] = u[d][im]*dmy_phi;
	  force[d] += dmy_fp[d];
	}
	{// torque
	  torque[0] += (r[1] * dmy_fp[2] - r[2] * dmy_fp[1]);
	  torque[1] += (r[2] * dmy_fp[0] - r[0] * dmy_fp[2]);
	  torque[2] += (r[0] * dmy_fp[1] - r[1] * dmy_fp[0]);
	}
      }// mesh
      

      pspec = p[n].spec;
      double dmy_mass2;
      if(DBG_MASS_GRID){//correct mass grid resolution
        dmy_mass = MASS[pspec] / (p[n].mass * DX3 * RHO_particle[pspec]);      
        dmy_mass2 = (p[n].inertia[0][0] + p[n].inertia[1][1] + p[n].inertia[2][2])/3.0;
        dmy_mass2 = MOI[pspec] / (dmy_mass2 * DX3 * RHO_particle[pspec]);
      }else{
        dmy_mass = 1.0;
        dmy_mass2 = 1.0;
      }

      for(int d = 0; d < DIM; d++){
	p[n].f_slip[d] = (dmy * force[d]) * dmy_mass;
        p[n].torque_slip[d] = (dmy * torque[d]) * dmy_mass2;
      }
    }//Particle_number
}
void Calc_f_hydro_correct_precision(Particle *p, double const* const* u, const CTime &jikan){
    static const double dmy0 = -DX3*RHO;
    double dmy = dmy0 /jikan.dt_fluid; 
    
    int *nlattice;
    //if (SW_EQ == Shear_Navier_Stokes){
    //nlattice = Ns_shear;
    //}else {
    nlattice = Ns;
    //}
    double xp[DIM],vp[DIM],omega_p[DIM];
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell;
    double force[DIM];
    double torque[DIM];
    int r_mesh[DIM];
    double r[DIM];
    double dmy_fp[DIM];
    double x[DIM];
    double dmyR;
    double dmy_phi;
    double v_rot[DIM];
    double dmy_mass;
    double dmy_inertia[DIM];
    int pspec;

    double forceg[DIM];
    double torqueg[DIM];
    int rigidID;
    // initialize forceGs and torqueGs
    for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
      for(int d=0; d<DIM; d++){
        forceGs[rigidID][d] = 0.0;
        torqueGs[rigidID][d] = 0.0;
      }
    }

#pragma omp parallel for schedule(dynamic, 1) \
  private(xp,vp,omega_p,x_int,residue,sw_in_cell,force,torque,r_mesh,r,dmy_fp,x,dmyR,dmy_phi,v_rot,dmy_mass,dmy_inertia,pspec, \
          rigidID, forceg, torqueg) 
    for(int n = 0; n < Particle_Number ; n++){
	//double xp[DIM],vp[DIM],omega_p[DIM];
	//int x_int[DIM];
	//double residue[DIM];
	if(SW_PT == rigid) rigidID = Particle_RigidID[n];
	for (int d = 0; d < DIM; d++) {
	    xp[d] = p[n].x[d];
	    vp[d] = p[n].v[d];
	    omega_p[d] = p[n].omega[d];

	    force[d] = torque[d] = dmy_inertia[d] = 0.0;
            forceg[d] = torqueg[d] = 0.0;
	}
	dmy_mass = 0.0;

	sw_in_cell 
	    = Particle_cell(xp, DX, x_int, residue);// {1,0} $B$,JV$C$F$/$k(B
	sw_in_cell = 1;
	for(int mesh=0; mesh < NP_domain; mesh++){
	    Relative_coord(Sekibun_cell[mesh], x_int, residue, sw_in_cell, nlattice, DX, r_mesh, r);
	    for(int d=0;d<DIM;d++){
		x[d] = r_mesh[d] * DX;
		//dmyR += SQ(r[d]);
	    }
	    //dmyR = sqrt(dmyR); // vesion2.10 needs this value
	    dmyR = Distance(x, xp); // vesion2.00 needs this value
	    dmy_phi= Phi(dmyR, RADIUS);
	    Angular2v(omega_p, r, v_rot);

	    //Mass grid factors
	    dmy_mass += dmy_phi;
	    for(int d = 0; d < DIM; d++){
	      dmy_inertia[d] += dmy_phi * (dmyR * dmyR - r[d]*r[d]);
	    }
            
	    int im = (r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2];
	    for(int d=0; d < DIM; d++ ){ 
              dmy_fp[d] = ((vp[d]+v_rot[d]) - u[d][im])*dmy_phi;
              force[d] += dmy_fp[d];
	    }
	    {// torque
              torque[0] += (r[1] * dmy_fp[2] - r[2] * dmy_fp[1]);
              torque[1] += (r[2] * dmy_fp[0] - r[0] * dmy_fp[2]);
              torque[2] += (r[0] * dmy_fp[1] - r[1] * dmy_fp[0]);
	    }
            
            if(SW_PT == rigid){
              for(int d = 0; d < DIM; d++) forceg[d] += dmy_fp[d];
              torqueg[0] += ( (GRvecs[n][1] + r[1]) * dmy_fp[2] - (GRvecs[n][2] + r[2]) * dmy_fp[1] );
              torqueg[1] += ( (GRvecs[n][2] + r[2]) * dmy_fp[0] - (GRvecs[n][0] + r[0]) * dmy_fp[2] );
              torqueg[2] += ( (GRvecs[n][0] + r[0]) * dmy_fp[1] - (GRvecs[n][1] + r[1]) * dmy_fp[0] );
            }
        }// mesh

        pspec = p[n].spec;
        double dmy_mass2;
        if(DBG_MASS_GRID){ //correct mass grid resolution
          dmy_mass = MASS[pspec] / (dmy_mass * DX3 * RHO_particle[pspec]);
          dmy_mass2 = (dmy_inertia[0] + dmy_inertia[1] + dmy_inertia[2])/3.0;
          dmy_mass2 = MOI[pspec] / (dmy_mass2 * DX3 * RHO_particle[pspec]);
        }else{
          dmy_mass = 1.0;
          dmy_mass2 =1.0;
        }

        for(int d = 0; d < DIM; d++){
          p[n].f_hydro[d] = (dmy * force[d]) * dmy_mass;
          p[n].f_slip[d] = 0.0;
	}
	if(ROTATION){
	  for(int d = 0; d < DIM; d++){
	    p[n].torque_hydro[d] = (dmy * torque[d]) * dmy_mass2;
	    p[n].torque_slip[d] = 0.0;
	  }
	}
	if(SW_PT == rigid){
          for(int d=0; d<DIM; d++){
            #pragma omp atomic
            forceGs[rigidID][d] += dmy * forceg[d];

            #pragma omp atomic
            torqueGs[rigidID][d] += dmy * torqueg[d];
          }
	}
    }//Particle_Number
}

void Calc_f_hydro_correct_precision_OBL(Particle *p, double const* const* u, const CTime &jikan){
    static const double dmy0 = -DX3*RHO;
    double dmy = dmy0 /jikan.dt_fluid; 
    
    int *nlattice;
    nlattice = Ns;
    
    double xp[DIM],vp[DIM],omega_p[DIM];
    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell;
    double force[DIM];
    double torque[DIM];
    int r_mesh[DIM];
    double r[DIM];
    double dmy_fp[DIM];
    double x[DIM];
    double dmyR;
    double dmy_phi;
    double dmy_rhop;
    double dmy_ry;
    double v_rot[DIM];
    double volume[Particle_Number];
    double Itrace[Particle_Number];
    int sign;
    int im;
    double sum_force = 0.0;
    double sum_volume = 0.0;
    
    double forceg[DIM];
    double torqueg[DIM];
    double dVg[Rigid_Number][DIM];
    double dWg[Rigid_Number][DIM];
    // initialize forceGs and torqueGs
    for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
      for(int d=0; d<DIM; d++){
        forceGs[rigidID][d] = 0.0;
        torqueGs[rigidID][d] = 0.0;
      }
    }
    
    Reset_phi(Hydro_force);
    Reset_phi(Hydro_force_new);
    if(DBG_LE_SHEAR){
      Reset_phi(Hydro_force_new_u);
      Reset_phi(Hydro_force_new_p);
      Reset_phi(Hydro_force_new_v);
      Reset_phi(Hydro_force_new_w);
    }
#pragma omp parallel for schedule(dynamic, 1)\
  private(xp,vp,omega_p,x_int,residue,sw_in_cell,force,torque,r_mesh,r,dmy_fp,x,\
          dmyR,dmy_phi,dmy_rhop,dmy_ry,v_rot,sign,im,forceg,torqueg) 
    for(int n = 0; n < Particle_Number ; n++){
      int rigidID;
      dmy_rhop = RHO_particle[p[n].spec];
      if(SW_PT == rigid) rigidID = Particle_RigidID[n];
      
      for (int d = 0; d < DIM; d++){
        xp[d] = p[n].x[d];
        vp[d] = p[n].v[d];
        omega_p[d] = p[n].omega[d];
        
        force[d] = torque[d] = 0.0;
        forceg[d] = torqueg[d] = 0.0;
      }
      
      volume[n] = 0.0;
      Itrace[n] = 0.0;
      sw_in_cell 
        = Particle_cell(xp, DX, x_int, residue);// {1,0} $B$,JV$C$F$/$k(B
      sw_in_cell = 1;
      
      for(int mesh=0; mesh < NP_domain; mesh++){
        sign = Relative_coord_check_stepover_Y(Sekibun_cell[mesh], x_int, residue, 
                                               sw_in_cell, nlattice, DX, r_mesh, r);
        dmyR = 0.;
        for(int d=0;d<DIM;d++){
          x[d] = r_mesh[d] * DX;
          dmyR += SQ(r[d]);
        }
        dmyR = sqrt(dmyR);
        dmy_phi= Phi(dmyR, RADIUS);
        Angular2v(omega_p, r, v_rot);
	
        im = (r_mesh[0]*NY*NZ_ + r_mesh[1]*NZ_ + r_mesh[2]);
        for(int d=0; d < DIM; d++ ){ 
          if (!(d==0)) {
            dmy_fp[d] =	((vp[d]+v_rot[d]) - u[d][im])*dmy_phi;
          } else {
            dmy_fp[d] =	((vp[d] - sign*Shear_rate_eff*L_particle[1] + v_rot[d]) 
                         - u[d][im])*dmy_phi;
          }
          force[d] += dmy_fp[d];
        }
        {// torque
          torque[0] += (r[1] * dmy_fp[2] - r[2] * dmy_fp[1]);
          torque[1] += (r[2] * dmy_fp[0] - r[0] * dmy_fp[2]);
          torque[2] += (r[0] * dmy_fp[1] - r[1] * dmy_fp[0]);
        }
        if(SW_PT == rigid){
          for(int d=0; d<DIM; d++) forceg[d] += dmy_fp[d];
          torqueg[0] += ( (GRvecs[n][1] + r[1]) * dmy_fp[2] - (GRvecs[n][2] + r[2]) * dmy_fp[1] );
          torqueg[1] += ( (GRvecs[n][2] + r[2]) * dmy_fp[0] - (GRvecs[n][0] + r[0]) * dmy_fp[2] );
          torqueg[2] += ( (GRvecs[n][0] + r[0]) * dmy_fp[1] - (GRvecs[n][1] + r[1]) * dmy_fp[0] );
        }
	
        dmy_ry = (r_mesh[1] + sign*L_particle[1]);
#pragma omp atomic
        Hydro_force[im] += dmy_fp[0]*dmy_ry*dmy_rhop;//viscocity
#pragma omp atomic
        Hydro_force_new[im] += dmy_fp[0]*dmy_ry*dmy_rhop;
        
        if(DBG_LE_SHEAR){
#pragma omp atomic
          Hydro_force_new_u[im] += dmy_fp[0]*dmy_ry*dmy_rhop;
        }

        volume[n] += dmy_phi;
        Itrace[n] += dmy_phi*SQ(dmyR);
      }//mesh
      Itrace[n] *= 2.0/3.0;
      
      for(int d=0; d < DIM; d++ ){ 
        p[n].f_hydro[d] = (dmy * force[d]);
      }
      if(ROTATION){
        for(int d=0; d < DIM; d++ ){ 
          p[n].torque_hydro[d] = ( dmy * torque[d]);
        }
      }
      if(SW_PT == rigid){
        for(int d=0; d<DIM; d++){
#pragma omp atomic
          forceGs[rigidID][d] += dmy * forceg[d];
          
#pragma omp atomic
          torqueGs[rigidID][d] += dmy * torqueg[d];
        }
      }

      for(int mesh=0; mesh < NP_domain; mesh++){
        sign = Relative_coord_check_stepover_Y(Sekibun_cell[mesh], x_int, residue,
                                               sw_in_cell, nlattice, DX, r_mesh, r);
        dmyR = 0.;
        for(int d=0;d<DIM;d++){
          dmyR += SQ(r[d]);
        }
        dmyR = sqrt(dmyR);
        dmy_phi= Phi(dmyR, RADIUS);
        
        im = (r_mesh[0]*NY*NZ_ + r_mesh[1]*NZ_ + r_mesh[2]);
        dmy_ry = (r_mesh[1] + sign*L_particle[1]);
        double dmy_stress_p = p[n].momentum_depend_fr[0] / volume[n];

#pragma omp atomic
        Hydro_force[im] -= dmy_stress_p*dmy_ry*dmy_phi;//viscocity

        if(SW_PT != rigid){
          double dmy_stress_v = -RHO*force[0] / volume[n];
          double dmy_stress_w = -RHO*(torque[1]*r[2] - torque[2]*r[1])/Itrace[n];
#pragma omp atomic
          Hydro_force_new[im] += (dmy_stress_v + dmy_stress_w)*dmy_ry*dmy_phi;
          if(DBG_LE_SHEAR){
#pragma omp atomic
            Hydro_force_new_p[im] += dmy_stress_p*dmy_ry*dmy_phi;
#pragma omp atomic
            Hydro_force_new_v[im] += dmy_stress_v*dmy_ry*dmy_phi;
#pragma omp atomic
            Hydro_force_new_w[im] += dmy_stress_w*dmy_ry*dmy_phi;
          }
        }
      }
      
# pragma omp critical
      {
        sum_force += force[0];
        sum_volume += volume[n];
      }
    }// Particle_Number
    sum_volume /= RHO;

    if(SW_PT == rigid){
#pragma omp parallel for schedule(dynamic, 1)
      for(int rigidID = 0; rigidID < Rigid_Number; rigidID++){
        for(int d = 0; d < DIM; d++){
          dVg[rigidID][d] = jikan.dt_md * Rigid_IMasses[rigidID] * forceGs[rigidID][d];
          dWg[rigidID][d] = jikan.dt_md * (Rigid_IMoments[rigidID][d][0] * torqueGs[rigidID][0] +
                                           Rigid_IMoments[rigidID][d][1] * torqueGs[rigidID][1] +
                                           Rigid_IMoments[rigidID][d][2] * torqueGs[rigidID][2]);
        }
      }
    }
    
#pragma omp parallel for schedule(dynamic, 1) \
  private(xp,x_int,residue,sw_in_cell,r_mesh,r,x,dmyR,dmy_phi,dmy_ry,sign,im,dmy_rhop) 
    for(int n = 0; n < Particle_Number ; n++){
      int rigidID;
      dmy_rhop = RHO_particle[p[n].spec];
      for (int d = 0; d < DIM; d++) {
        xp[d] = p[n].x[d];
      }
      sw_in_cell = Particle_cell(xp, DX, x_int, residue);// {1,0} Ç™ï‘Ç¡ÇƒÇ≠ÇÈ
      sw_in_cell = 1;
      
      for(int mesh=0; mesh < NP_domain; mesh++){
        sign = Relative_coord_check_stepover_Y(Sekibun_cell[mesh], x_int, residue, \
                                               sw_in_cell, nlattice, DX, r_mesh, r);

        dmyR = 0;
        for(int d=0;d<DIM;d++){
          x[d] = r_mesh[d] * DX;
          dmyR += SQ(r[d]);
        }
        dmyR= sqrt(dmyR);
        dmy_phi = Phi(dmyR, RADIUS);

        im = (r_mesh[0]*NY*NZ_ + r_mesh[1]*NZ_ + r_mesh[2]);
        dmy_ry = (r_mesh[1] + sign*L_particle[1]);
#pragma omp atomic	
        Hydro_force[im] -= sum_force*dmy_ry*dmy_phi/sum_volume;//viscocity
        
        if(SW_PT == rigid){
          rigidID = Particle_RigidID[n];
          double dmy_stress_v = dVg[rigidID][0];
          double dmy_stress_w = dWg[rigidID][1]*(GRvecs[n][2] + r[2]) - dWg[rigidID][2]*(GRvecs[n][1] + r[1]);
#pragma omp atomic
          Hydro_force_new[im] += (dmy_stress_v + dmy_stress_w)*dmy_ry*dmy_phi*dmy_rhop;
          if(DBG_LE_SHEAR){
            Hydro_force_new_v[im] += dmy_stress_v*dmy_ry*dmy_phi*dmy_rhop;
            Hydro_force_new_w[im] += dmy_stress_w*dmy_ry*dmy_phi*dmy_rhop;
          }
        }
      }
    }//Particle_Number
}
