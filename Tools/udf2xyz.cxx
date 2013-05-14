#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include "../alloc.h"
#include "../macro.h"
#include "../quaternion.h"
#include "../rigid_body.h"
#include "../lad3.h"
#include "udfmanager.h"
#define MAXBUFFER 100
#define NDIM 3
using namespace std;
enum OPTION {in_format=1, in_fin=2, in_fout=3};

double ex[NDIM] = {1.0, 0.0, 0.0};
double ey[NDIM] = {0.0, 1.0, 0.0};
double ez[NDIM] = {0.0, 0.0, 1.0};
double e_none[NDIM] = {0.0, 0.0, 0.0};

int xyz_format;
const int max_spec = 5;
const char* alphabet[]={"H", "O", "C", "N", "K"};
const char* in_option[]={"-dump", "-xyz"};

void newframe(ofstream &outfile,
	      double &t,
	      int &ntot);

void write_xyz(ofstream &outfile, 
               const int &pid, 
               const int &sid, 
               double r[NDIM]);

void write_xyz_dump(ofstream &outfile, 
                    double &t, 
                    const int &pid, 
                    const int &sid,
                    double r[NDIM], 
                    double r_raw[NDIM], 
                    const quaternion &q,
                    double QR[NDIM][NDIM], 
                    double v[NDIM], 
                    double vs[NDIM],
                    double w[NDIM], 
                    double frc[NDIM], 
                    double frc_slip[NDIM],
                    double tau[NDIM], 
                    double tau_slip[NDIM], 
                    double ni[NDIM],
                    double n0[NDIM],
                    int &ntot);

void init_xyz(ofstream &outfile);

void close_xyz(ofstream &outfile);


inline void wrong_invocation(){
  cout << "Converts udf file format to simple xyz format." << endl;
  cout << "Usage: get_xyz format UDFfile XYZfile" << endl;
  cout << " format:   -dump      dump all particle data"  << endl;
  cout << "           -xyz       xyz file with positions only" << endl;
  exit(1);
}


// convert udf trajectory file to xyz file
int main(int argc, char* argv[])
{
  
  //check invocation

  if(argc != 4){
    cout << "Error: " << endl;
    wrong_invocation();
  }

  //make sure udf file exists
  UDFManager *ufin;
  if(file_check(argv[in_fin])){
    ufin = new UDFManager(argv[in_fin]);
  }
  if(strcmp(argv[in_format],in_option[0]) == 0){
    cout << "Dumping trj data..." << endl;
    xyz_format = 0;
  }else if(strcmp(argv[in_format],in_option[1]) == 0){
    cout << "Generating xyz file..." << endl;
    xyz_format = 1;
  }else{
    cout << "Unknown format option" << endl;
    wrong_invocation();
  }

  ofstream outfile;
  outfile.open(argv[in_fout]);
  init_xyz(outfile);  

  int nx, ny, nz;
  int dt, frames, records;
  int ntotal;
  int nspec;
  int *pnum;
  int *spec_id;
  double **janus_axis;
  double *uslip;
  double uscale;
  {
    //mesh size
    Location target("mesh");
    int npx, npy, npz;
    ufin -> get(target.sub("NPX"), npx);
    ufin -> get(target.sub("NPY"), npy);
    ufin -> get(target.sub("NPZ"), npz);
    nx = 1 << npx;
    ny = 1 << npy;
    nz = 1 << npz;
  }
  {
    //particle species data
    Location target("object_type");
    string str;
    ufin -> get(target.sub("type"), str);
    if(str != "spherical_particle"){
      cout << "Error unknown particle type" << endl;
      exit(1);
    }
    nspec = ufin -> size("object_type.spherical_particle.Particle_spec[]");
    if(nspec > max_spec && xyz_format){
      cout << "Increase max_spec" << endl;
      exit(1);
    }

    pnum = alloc_1d_int(nspec);
    uslip = alloc_1d_double(nspec);
    janus_axis = (double **) malloc(sizeof(double*) * nspec);
    int nmax = 0;
    ntotal = 0;
    uscale = 0.0;
    for(int i = 0; i < nspec; i++){
      char str[256];
      string dmy_str;
      sprintf(str, "object_type.spherical_particle.Particle_spec[%d]",i);
      Location target(str);
      ufin -> get(target.sub("Particle_number"), pnum[i]);
      if(!xyz_format){
        ufin -> get(target.sub("janus_slip_vel"), uslip[i]);
        ufin -> get(target.sub("janus_axis"), dmy_str);
        if(dmy_str == "X"){
          janus_axis[i] = ex;
        }else if(dmy_str == "Y"){
          janus_axis[i] = ey;
        }else if(dmy_str == "Z"){
          janus_axis[i] = ez;
        }else if(dmy_str == "NONE"){
          janus_axis[i] = e_none;
        }else{
          fprintf(stderr, "Error: Unknown janus axis\n");
          exit_job(EXIT_FAILURE);
        }
        fprintf(stderr, "JANUS AXIS: %3.2f %3.2f %3.2f\n", 
                janus_axis[i][0], janus_axis[i][1], janus_axis[i][2]);
        fprintf(stderr, "V_SLIP: %10.8g ", uslip[i]);
        uslip[i] = (ABS(uslip[i]) > 0.0 ? 1.0 / (2./3 * uslip[i]) : 0.0);
        uscale = MAX(ABS(uslip[i]), uscale);
        fprintf(stderr, " (%10.8g) \n", uslip[i]);
      }
      nmax = MAX(nmax, pnum[i]);
      ntotal += pnum[i];
    }
    if(uscale == 0.0) uscale = 1.0;
    fprintf(stderr, "V_SCALE: %10.8g \n", uscale);

    spec_id = (int*) malloc(sizeof(int) * ntotal);
    int sum = 0;
    for(int i = 0; i < nspec; i++){
      for(int j = 0; j < pnum[i]; j++){
	spec_id[sum] = i;
	sum++;
      }
    }
    assert(sum == ntotal);
  }
  {
    //frame data
    Location target("output");
    ufin -> get(target.sub("GTS"), dt);
    ufin -> get(target.sub("Num_snap"), frames);
    records = ufin -> totalRecord();
  }
  {
    //read particle data
    double r[NDIM];
    double r_raw[NDIM];
    double v[NDIM];
    double vs[NDIM];
    double w[NDIM];
    double QR[NDIM][NDIM];
    double frc[NDIM];
    double frc_slip[NDIM];
    double tau[NDIM];
    double tau_slip[NDIM];
    double n0[NDIM];
    double ni[NDIM];
    quaternion q;
    double t;
    double q0,q1,q2,q3;
    for(int i = 0; i < records; i++){
      char str[256];
      sprintf(str, "#%d", i);
      ufin -> jump(str);
      ufin -> get("t", t);

      newframe(outfile, t, ntotal);
      for(int j = 0; j < ntotal; j++){
	char str[256];
	sprintf(str, "Particles[%d]", j);
	Location target(str);

        if(!xyz_format){
          // position
          ufin -> get(target.sub("R.x"), r[0]);
          ufin -> get(target.sub("R.y"), r[1]);
          ufin -> get(target.sub("R.z"), r[2]);
          
          ufin -> get(target.sub("R_raw.x"), r_raw[0]);
          ufin -> get(target.sub("R_raw.y"), r_raw[1]);
          ufin -> get(target.sub("R_raw.z"), r_raw[2]);
          
          // orientation
          ufin -> get(target.sub("q.q0"), q0);
          ufin -> get(target.sub("q.q1"), q1);
          ufin -> get(target.sub("q.q2"), q2);
          ufin -> get(target.sub("q.q3"), q3);
          qtn_init(q, q0 ,q1, q2, q3);
          qtn_normalize(q);
          
          // velocity
          ufin -> get(target.sub("v.x"), v[0]);
          ufin -> get(target.sub("v.y"), v[1]);
          ufin -> get(target.sub("v.z"), v[2]);
          for(int d = 0; d < NDIM; d++){
            vs[d] = v[d] * uscale;
          }
          
          // angular velocity
          ufin -> get(target.sub("omega.x"), w[0]);
          ufin -> get(target.sub("omega.y"), w[1]);
          ufin -> get(target.sub("omega.z"), w[2]);
          
          // forces
          ufin -> get(target.sub("f_hydro.x"), frc[0]);
          ufin -> get(target.sub("f_hydro.y"), frc[1]);
          ufin -> get(target.sub("f_hydro.z"), frc[2]);
          
          ufin -> get(target.sub("f_slip.x"), frc_slip[0]);
          ufin -> get(target.sub("f_slip.y"), frc_slip[1]);
          ufin -> get(target.sub("f_slip.z"), frc_slip[2]);
          
          // torques
          ufin -> get(target.sub("torque_hydro.x"), tau[0]);
          ufin -> get(target.sub("torque_hydro.y"), tau[1]);
          ufin -> get(target.sub("torque_hydro.z"), tau[2]);
          
          ufin -> get(target.sub("torque_slip.x"), tau_slip[0]);
          ufin -> get(target.sub("torque_slip.y"), tau_slip[1]);
          ufin -> get(target.sub("torque_slip.z"), tau_slip[2]);
          
          rqtn_rm(QR, q);	
          rigid_body_rotation(ni, janus_axis[spec_id[j]], q, BODY2SPACE);
          if(i == 0){
            v_copy(n0, ni);
          }
          write_xyz_dump(outfile, t, j+1, spec_id[j], r, r_raw, q, QR, v, 
                         vs, w, frc, frc_slip, tau, tau_slip, ni, n0, ntotal);
        }else{//xyz_format
          ufin -> get(target.sub("R.x"), r[0]);
          ufin -> get(target.sub("R.y"), r[1]);
          ufin -> get(target.sub("R.z"), r[2]);
          ufin -> get(target.sub("v.x"), v[0]);
          ufin -> get(target.sub("v.y"), v[1]);
          ufin -> get(target.sub("v.z"), v[2]);
            write_xyz(outfile, j+1, spec_id[j], r);
        }
      }//j
    }//i
  }

  free_1d_int(pnum);
  close_xyz(outfile);
  return 0;
}

void newframe(ofstream &outfile,
              double &t,
              int &ntot){
  if(xyz_format){
    char dmy_str[256];
    sprintf(dmy_str, "%d\n", ntot);
    outfile << dmy_str << endl;
  }else{
    //  outfile << endl;
  }
}
void init_xyz(ofstream &outfile){

  if(xyz_format){
  }else{
    char str[512];    
    
    //id
    sprintf(str, "## 1:p_id 2:spec_id ");
    
    //position
    strcat(str, "3:r_x 4:r_y 5:r_z ");
    
    //position no pbc
    strcat(str, "6:rr_x 7:rr_y 8:rr_z ");
    
    //velocities
    strcat(str, "9:v_x 10:v_y 11:v_z 12:v 13:vs ");
    
    //orientation matrix
    strcat(str, "14:e1_x 15:e1_y 16:e1_z 17:e2_x 18:e2_y 19:e2_z 20:e3_x 21:e3_y 22:e3_z ");
    
    //angular velocity
    strcat(str, "23:w_x 24:w_y 25:w_z ");
    
    // forces (hydro/slip)
    strcat(str, "26:Fh_x 27:Fh_y 28:Fh_z 29:Fs_x 30:Fs_y 31:Fs_z 32:F_x 33:F_y 34:F_z");
    
    // torques (hydro/slip)
    strcat(str, "35:Nh_x 36:Nh_y 37:Nh_z 38:Ns_x 39:Ns_y 40:Ns_z 41:N_x 42:N_y 43:N_z");
    
    // janus axis
    strcat(str, "44:e_x 45:e_y 46:e_z 47:c_ee 48:theta_0 49:alpha 50:beta 51:gamma ");
    
    outfile << str << endl;
  }
}
void write_xyz(ofstream &outfile, const int &pid, const int &sid, double r[NDIM]){
  char str[256];
  sprintf(str, "%5s %.6g %.6g %.6g", alphabet[sid], r[0], r[1], r[2]);
  outfile << str << endl;
}
void write_xyz_dump(ofstream &outfile, 
                    double &t,
                    const int &pid,
                    const int &sid,
                    double r[NDIM],
                    double r_raw[NDIM],
                    const quaternion &q,
                    double QR[NDIM][NDIM],
                    double v[NDIM],
                    double vs[NDIM],
                    double w[NDIM],
                    double frc[NDIM],
                    double frc_slip[NDIM],
                    double tau[NDIM],
                    double tau_slip[NDIM],
                    double ni[NDIM],
                    double n0[NDIM],
                    int &ntot
                    ){
  char str[4096];
  char dmy_str[256];
  double dmy_ndot = ni[0]*n0[0] + ni[1]*n0[1] + ni[2]*n0[2];
  double dmy_theta = (dmy_ndot < 1.0 ? acos(dmy_ndot) : 0.0);
  double dmy_theta_z = (ni[2] < 1.0 ? acos(ni[2]) : 0.0);

  // particle info
  sprintf(str, "%d %d ", pid, sid);

  // positions
  sprintf(dmy_str, "%.6g  %.6g  %.6g  %.6g  %.6g %.6g ", 
          r[0], r[1], r[2], r_raw[0], r_raw[1], r_raw[2]);
  strcat(str, dmy_str);

  // velocities
  sprintf(dmy_str, "%.6g  %.6g  %.6g  %.6g  %.6g ", v[0], v[1], v[2],
	  sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]),
	  sqrt(vs[0]*vs[0] + vs[1]*vs[1] + vs[2]*vs[2])
	  );
  strcat(str, dmy_str);

  // orientation matrix
  sprintf(dmy_str, "%.6g  %.6g  %.6g  %.6g  %.6g  %.6g  %.6g  %.6g  %.6g  ", 
	  QR[0][0], QR[1][0], QR[2][0],
	  QR[0][1], QR[1][1], QR[2][1],
	  QR[0][2], QR[1][2], QR[2][2]);
  strcat(str, dmy_str);

  // angular velocity
  sprintf(dmy_str, "%.6g  %.6g  %.6g  ", w[0], w[1], w[2]);
  strcat(str, dmy_str);

  // forces and torques
  sprintf(dmy_str, "%.6g  %.6g  %.6g  %.6g  %.6g  %.6g  %.6g  %.6g  %.6g  ", 
          frc[0], frc[1], frc[2], frc_slip[0], frc_slip[1], frc_slip[2], 
          frc[0] + frc_slip[0], frc[1] + frc_slip[1], frc[2] + frc_slip[2]);
  strcat(str, dmy_str);
  sprintf(dmy_str, "%.6g  %.6g  %.6g  %.6g  %.6g  %.6g  %.6g  %.6g  %.6g  ", 
          tau[0], tau[1], tau[2], tau_slip[0], tau_slip[1], tau_slip[2], 
          tau[0] + tau_slip[0], tau[1] + tau_slip[1], tau[2] + tau_slip[2]);
  strcat(str, dmy_str);

  // janus vector
  double alpha, beta, gamma;
  double dmy_angle = 180.0 / M_PI;
  rqtn_euler(alpha, beta, gamma, q);
  if(!equal_tol(dmy_theta_z, beta, HUGE_TOL_MP)){
    fprintf(stderr, "Check euler conversion: %g %g\n", dmy_theta_z, beta);
  }
  sprintf(dmy_str, "%.6g  %.6g  %.6g  %.15g  %.6g  %.6g  %.6g  %.6g  ", 
          ni[0], ni[1], ni[2], dmy_ndot, 
          dmy_theta*dmy_angle, alpha*dmy_angle, beta*dmy_angle, gamma*dmy_angle);
  strcat(str, dmy_str);

  outfile << str << endl;
}


void close_xyz(ofstream &outfile){
  outfile.close();
}
