/*!
  \file avs_output_p.cxx
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
  \brief Output routines for particle data in AVS/Express format
 */
#include "avs_output_p.h"

const int Veclen_p = 10;
const char *Label_p="Radius x1 x2 x3 y1 y2 y3 z1 z2 z3";

void Init_avs_p(const AVS_parameters &Avs_parameters)
{
  FILE *fout;
  fout = filecheckopen(Avs_parameters.pfld_file, "w");
  fprintf(fout, "# AVS field file\n");
  fprintf(fout, "ndim=%d\n", 1);
  fprintf(fout, "dim1=%d\n", Particle_Number);
  fprintf(fout, "nspace=%d\n", DIM);
  fprintf(fout, "veclen=%d\n", Veclen_p);
  fprintf(fout, "data=float\n");
  fprintf(fout, "field=irregular\n");
  fprintf(fout, "nstep=%d\n", Avs_parameters.nstep);
  fprintf(fout, "label = %s\n", Label_p);
  fclose(fout);
}

inline void Add_pfield_description(AVS_parameters &Avs_parameters
                                   , const CTime &time
                                   , const int &veclen
				   )
{
  FILE *fout;
  char line[512];
  fout = filecheckopen(Avs_parameters.pfld_file, "a");
  fprintf(fout, "time value = \"step%dtime%g\"\n"
	  , time.ts, time.time);
  if(BINARY) {
    static const int data_size = sizeof(float) *
      Particle_Number;// * DIM;
    for(int n = 0; n < DIM; n++) {
      fprintf(fout,
	      "coord %d file = %s%d.cod filetype = binary skip = %d\n",
	      DIM - n,
	      Avs_parameters.out_ppfx, time.ts, n * data_size);
    }
    for(int n = 0; n < veclen; n++) {
      fprintf(fout,
	      "variable %d file = %s%d.dat filetype = binary skip = %d\n",
	      n + 1,
	      Avs_parameters.out_ppfx, time.ts, n * data_size);
    }
  }
  else {
    for(int n = 0; n < DIM; n++) {
      fprintf(fout,
	      "coord %d file = %s%d.cod filetype = ascii skip = 0 offset = %d stride =%d\n",
	      DIM - n,
	      Avs_parameters.out_ppfx, time.ts, n, DIM);
    }
    for(int n = 0; n < veclen; n++) {
      fprintf(fout,
	      "variable %d file = %s%d.dat filetype = ascii skip = 1 offset = %d stride = %d\n",
	      n + 1, Avs_parameters.out_ppfx, time.ts, n, veclen);
    }
  }
  fprintf(fout, "EOT\n");
  fclose(fout);
}

void Output_avs_p(AVS_parameters &Avs_parameters
                  , Particle *p
                  , const CTime &time
		  )
{
  Add_pfield_description(Avs_parameters, time, Veclen_p);
  
  FILE *fout;
  char line[512];
  sprintf(Avs_parameters.data_file, "%s/%s%d.cod",
	  Out_dir, Avs_parameters.out_ppfx, time.ts);
  fout = filecheckopen(Avs_parameters.data_file, "wb");

  if(BINARY) {
    float dmy;
    for(int d = 0; d < DIM; d++) {
      for(int n = 0; n < Particle_Number; n++) {
	dmy = (float)p[n].x[DIM - 1 - d];
	fwrite(&dmy, sizeof(float), 1, fout);
      }
    }
  }
  else {
    for(int n = 0; n < Particle_Number; n++) {
      fprintf(fout, "%.3g %.3g %.3g\n"
	      , p[n].x[2]
	      , p[n].x[1]
	      , p[n].x[0]
	      );
    }
  }
  fclose(fout);

  sprintf(line, "timesteps=%d time=%f\n%s", time.ts, time.time, Label_p);
  sprintf(Avs_parameters.data_file, "%s/%s%d.dat",
	  Out_dir, Avs_parameters.out_ppfx, time.ts);
  fout = filecheckopen(Avs_parameters.data_file, "wb");

  //Precompute orientation matrix
  for(int n = 0; n < Particle_Number; n++){
    qtn_isnormal(p[n].q);
    rqtn_rm(p[n].QR, p[n].q);
  }

  if(BINARY) {
    float dmy;
    for(int n = 0; n < Particle_Number; n++) {   //radius
      dmy = (float)RADIUS;
      fwrite(&dmy, sizeof(float), 1, fout);
    }
    for(int d1 = 0; d1 < DIM; d1++){              //body coordinate system (ex, ey, ez)
      for(int d2 = 0; d2 < DIM; d2++){
	for(int n = 0; n < Particle_Number; n++){
	  dmy = (float)p[n].QR[d2][d1];
	  fwrite(&dmy, sizeof(float), 1, fout);
	}
      }
    }

  }
  else {
    fprintf(fout, "%s\n", line);
    for(int n = 0; n < Particle_Number; n++) {
      fprintf(fout, "%.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g\n",
	      RADIUS,
	      p[n].QR[0][0], p[n].QR[1][0], p[n].QR[2][0],
	      p[n].QR[0][1], p[n].QR[1][1], p[n].QR[2][1],
	      p[n].QR[0][2], p[n].QR[1][2], p[n].QR[2][2]);
    }
  }
  fclose(fout);
}


