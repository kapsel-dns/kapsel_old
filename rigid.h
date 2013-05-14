//
//	rigid.h 12/12/29 T.K
//
#ifndef RIGID_H
#define RIGID_H

#include "input.h"
#include "Matrix_Inverse.h"


inline void init_set_xGs(Particle *p){
	// initialize xGs[][]
#pragma omp parallel for schedule(dynamic, 1)
	for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
		for(int d=0; d<DIM; d++) xGs[rigidID][d] = 0.0;
	}
	
	int rigidID;
	
#pragma omp parallel for schedule(dynamic, 1) private(rigidID)
	for(int n=0; n<Particle_Number; n++){
		rigidID = Particle_RigidID[n];
		for(int d=0; d<DIM; d++) {
#pragma omp atomic
                  xGs[rigidID][d] += p[n].x[d];
                }
	}

#pragma omp parallel for schedule(dynamic, 1)
	for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
		for(int d=0; d<DIM; d++) xGs[rigidID][d] /= Rigid_Particle_Numbers[rigidID];
	}

#pragma omp parallel for schedule(dynamic, 1) private(rigidID)
	for(int n=0; n<Particle_Number; n++){
		rigidID = Particle_RigidID[n];
		for(int d=0; d<DIM; d++) {
                  GRvecs[n][d] = p[n].x[d] - xGs[rigidID][d];
                }
	}
}

inline void set_xGs(Particle *p){
	int rigid_first_n = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+:rigid_first_n)
	for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
		for(int d=0; d<DIM; d++){
			xGs[rigidID][d] = p[rigid_first_n].x[d] - GRvecs[rigid_first_n][d];
			xGs[rigidID][d] = fmod(xGs[rigidID][d] + L_particle[d], L_particle[d]);
		}
		rigid_first_n += Rigid_Particle_Numbers[rigidID];
	}
	
	//check(for debug)
	if(rigid_first_n != Particle_Number) fprintf(stderr, "debug: set_xGs() error");
}

inline void solver_GRvecs(const CTime &jikan, string CASE){
	int rigidID;
	double bufvec[DIM];
	if(CASE == "Euler"){
#pragma omp parallel for schedule(dynamic, 1) private(rigidID, bufvec)
		for(int n=0; n<Particle_Number; n++){
			rigidID = Particle_RigidID[n];
			bufvec[0] = (omegaGs[rigidID][1]*GRvecs[n][2] - omegaGs[rigidID][2]*GRvecs[n][1]) * jikan.dt_md;
			bufvec[1] = (omegaGs[rigidID][2]*GRvecs[n][0] - omegaGs[rigidID][0]*GRvecs[n][2]) * jikan.dt_md;
			bufvec[2] = (omegaGs[rigidID][0]*GRvecs[n][1] - omegaGs[rigidID][1]*GRvecs[n][0]) * jikan.dt_md;
			for(int d=0; d<DIM; d++) GRvecs[n][d] += bufvec[d];
		}
	}
	else if(CASE == "AB2"){
#pragma omp parallel for schedule(dynamic, 1) private(rigidID, bufvec)
		for(int n=0; n<Particle_Number; n++){
			rigidID = Particle_RigidID[n];
			bufvec[0] = ( (3.*omegaGs[rigidID][1]-omegaGs_old[rigidID][1])*GRvecs[n][2] - (3.*omegaGs[rigidID][2]-omegaGs_old[rigidID][2])*GRvecs[n][1] ) * jikan.hdt_md;
			bufvec[1] = ( (3.*omegaGs[rigidID][2]-omegaGs_old[rigidID][2])*GRvecs[n][0] - (3.*omegaGs[rigidID][0]-omegaGs_old[rigidID][0])*GRvecs[n][2] ) * jikan.hdt_md;
			bufvec[2] = ( (3.*omegaGs[rigidID][0]-omegaGs_old[rigidID][0])*GRvecs[n][1] - (3.*omegaGs[rigidID][1]-omegaGs_old[rigidID][1])*GRvecs[n][0] ) * jikan.hdt_md;
			for(int d=0; d<DIM; d++) GRvecs[n][d] += bufvec[d];
		}
	}
	else{
		fprintf(stderr, "error, string CASE in solver_GRvecs()");
		exit_job(EXIT_FAILURE);
	}
}

// set Rigid_Masses and Rigid_IMasses and Rigid_Moments and Rigid_IMoments
// dont use it before (init_set_xGs) or (set_xGs and solver_GRvecs)
inline void set_Rigid_MMs(Particle *p){
	// initialize Rigid_Masses and Rigid_Moments and Rigid_IMoments
#pragma omp parallel for schedule(dynamic, 1)
	for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
		Rigid_Masses[rigidID] = 0.0;
		for(int row=0; row<DIM; row++){
			for(int column=0; column<DIM; column++){
				Rigid_Moments[rigidID][row][column] = 0.0;
				Rigid_IMoments[rigidID][row][column] = 0.0;
			}
		}
	}
	
	int rigidID;
	
#pragma omp parallel for schedule(dynamic, 1) private(rigidID)
	for(int n=0; n<Particle_Number; n++){
		rigidID = Particle_RigidID[n];
		// mass of p[n] adds to Rigid_Masses[rigidID];
#pragma omp atomic
		Rigid_Masses[rigidID] += MASS[ RigidID_Components[rigidID] ];
		// moment of p[n] adds to Rigid_Moments[rigidID]
#pragma omp atomic
		Rigid_Moments[rigidID][0][0] += MOI[ RigidID_Components[rigidID] ]
                  + MASS[ RigidID_Components[rigidID] ] * ( GRvecs[n][1]*GRvecs[n][1] + GRvecs[n][2]*GRvecs[n][2] );
#pragma omp atomic
		Rigid_Moments[rigidID][0][1] += MASS[ RigidID_Components[rigidID] ] * ( - GRvecs[n][0]*GRvecs[n][1] );
#pragma omp atomic
		Rigid_Moments[rigidID][0][2] += MASS[ RigidID_Components[rigidID] ] * ( - GRvecs[n][0]*GRvecs[n][2] );
#pragma omp atomic
		Rigid_Moments[rigidID][1][0] += MASS[ RigidID_Components[rigidID] ] * ( - GRvecs[n][1]*GRvecs[n][0] );
#pragma omp atomic
		Rigid_Moments[rigidID][1][1] += MOI[ RigidID_Components[rigidID] ]
                  + MASS[ RigidID_Components[rigidID] ] * ( GRvecs[n][2]*GRvecs[n][2] + GRvecs[n][0]*GRvecs[n][0] );
#pragma omp atomic
		Rigid_Moments[rigidID][1][2] += MASS[ RigidID_Components[rigidID] ] * ( - GRvecs[n][1]*GRvecs[n][2] );
#pragma omp atomic
		Rigid_Moments[rigidID][2][0] += MASS[ RigidID_Components[rigidID] ] * ( - GRvecs[n][2]*GRvecs[n][0] );
#pragma omp atomic
		Rigid_Moments[rigidID][2][1] += MASS[ RigidID_Components[rigidID] ] * ( - GRvecs[n][2]*GRvecs[n][1] );
#pragma omp atomic
		Rigid_Moments[rigidID][2][2] += MOI[ RigidID_Components[rigidID] ]
                  + MASS[ RigidID_Components[rigidID] ] * ( GRvecs[n][0]*GRvecs[n][0] + GRvecs[n][1]*GRvecs[n][1] );
	}
	
	char str[256];
	// inverse
#pragma omp parallel for schedule(dynamic, 1)
	for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
		Rigid_IMasses[rigidID] = 1. / Rigid_Masses[rigidID];
		Matrix_Inverse(Rigid_Moments[rigidID], Rigid_IMoments[rigidID], DIM);
		
//		sprintf(str, "Rigid_Moments[%d].csv", rigidID);
//		Matrix_output(Rigid_Moments[rigidID], DIM, str); 
//		sprintf(str, "Rigid_IMoments[%d].csv", rigidID);
//		Matrix_output(Rigid_IMoments[rigidID], DIM, str); 
		
		check_Inverse(Rigid_Moments[rigidID], Rigid_IMoments[rigidID], DIM);
	}

}

inline void set_particle_vomegas(Particle *p){

	int rigidID;
#pragma omp parallel for schedule(dynamic, 1) private(rigidID)
	for(int n=0; n<Particle_Number; n++){
		rigidID = Particle_RigidID[n];
		for(int d=0; d<DIM; d++){
			p[n].omega[d] = omegaGs[rigidID][d];
		}
		p[n].v[0] = velocityGs[rigidID][0] + omegaGs[rigidID][1]*GRvecs[n][2] - omegaGs[rigidID][2]*GRvecs[n][1];
		p[n].v[1] = velocityGs[rigidID][1] + omegaGs[rigidID][2]*GRvecs[n][0] - omegaGs[rigidID][0]*GRvecs[n][2];
		p[n].v[2] = velocityGs[rigidID][2] + omegaGs[rigidID][0]*GRvecs[n][1] - omegaGs[rigidID][1]*GRvecs[n][0];
	}
}

// set VelocityGs and OmegaGs
// ### set_Rigid_VOGs() after calculating xGs, Rigid_IMoments, forceGs and torqueGs!! ###
inline void calc_Rigid_VOGs(Particle *p, const CTime &jikan, string CASE){
	set_Rigid_MMs(p);
	
	int rigidID;
	//calc forceGrs, forceGvs
#pragma omp parallel for schedule(dynamic, 1) private(rigidID)
	for(int n=0; n<Particle_Number; n++){
		rigidID = Particle_RigidID[n];
		for(int d=0; d<DIM; d++){
#pragma omp atomic
			forceGrs[rigidID][d] += p[n].fr[d];
			//forceGvs[rigidID][d] += p[n].fv[d];
			//torqueGrs[rigidID][d] += p[n].torquer[d];
			//torqueGvs[rigidID][d] += p[n].torquev[d];	//fv, torquer and torquev are constant zero.
		}
	}
	//set olds
#pragma omp parallel for schedule(dynamic, 1)
	for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
		for(int d=0; d<DIM; d++){
			velocityGs_old[rigidID][d] = velocityGs[rigidID][d];
			omegaGs_old[rigidID][d] = omegaGs[rigidID][d];
		}
	}
	
	//calc velocityGs and omegaGs
	if(CASE == "Euler"){
#pragma omp parallel for schedule(dynamic, 1)
		for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
			if(Rigid_Motions[ RigidID_Components[rigidID] ] == 0) continue;	// if "fix"
			for(int d1=0; d1<DIM; d1++){
				velocityGs[rigidID][d1] += jikan.dt_md * Rigid_IMasses[rigidID]
                                  * ( forceGs[rigidID][d1] + forceGrs[rigidID][d1] );
				for(int d2=0; d2<DIM; d2++) omegaGs[rigidID][d1] += jikan.dt_md * Rigid_IMoments[rigidID][d1][d2]
                                                              * ( torqueGs[rigidID][d2] );
			}
		}
	
	}else if(CASE == "AB2_hydro"){
#pragma omp parallel for schedule(dynamic, 1)
		for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
			if(Rigid_Motions[ RigidID_Components[rigidID] ] == 0) continue;	// if "fix"
			for(int d1=0; d1<DIM; d1++){
				velocityGs[rigidID][d1] += jikan.hdt_md * Rigid_IMasses[rigidID]
                                  * ( 2. * forceGs[rigidID][d1]
                                      + 3. * forceGvs[rigidID][d1] - forceGvs_previous[rigidID][d1]
                                      + forceGrs[rigidID][d1] + forceGrs_previous[rigidID][d1]
                                      );
				for(int d2=0; d2<DIM; d2++) omegaGs[rigidID][d1] += jikan.hdt_md * Rigid_IMoments[rigidID][d1][d2]
                                                              * ( 2. * torqueGs[rigidID][d2] );
			}
		}
	
	}else{
		fprintf(stderr, "error, string CASE in calc_Rigid_VOGs()");
		exit_job(EXIT_FAILURE);
	}


		
	// renew old previous and initialize
#pragma omp parallel for schedule(dynamic, 1)
	for(int rigidID=0; rigidID<Rigid_Number; rigidID++){
		for(int d=0; d<DIM; d++){
			 forceGrs_previous[rigidID][d] = forceGrs[rigidID][d];
			 forceGrs[rigidID][d] = 0.0;
			 forceGs[rigidID][d] = 0.0;
			 torqueGs[rigidID][d] = 0.0;
		}
	}
}
#endif
