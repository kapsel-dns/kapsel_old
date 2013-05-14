/*!
  \file macro.h
  \brief Global parameters and macro function definitions
  \details Defines various constants and low-level functions
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
 */

#ifndef MACRO_H_
#define MACRO_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>


const double Euler_cst = 0.57721566490153286060651209008240243104215933593992359880576723488485;

const double PI4 = 4. * M_PI;
const double PI2 = 2. * M_PI;
const double PI_half = .5 * M_PI;
const double One_third= 1./3.;
const double Root_six= sqrt(6.);
const double Root_two= sqrt(2.);
const double Root_three= sqrt(3.);
const double EPSILON_MP= DBL_EPSILON;// E-16
const double TOL_MP=10.0*EPSILON_MP;// E-15
const double LARGE_TOL_MP=1.0e3*TOL_MP;// E-12
const double HUGE_TOL_MP=1.0e6*LARGE_TOL_MP;//E-6
/////////////////////// inline function prototype if necessary
void exit_job(int status);
// 
/////////////////////// 
inline int Nint(const double &a){
    int ans;
    if(a>= 0.0){
	ans = int(a+.5);
    }else {
	ans = int(a-.5);
    }
    return ans;
}
/////////////////////// random number generator:fast linear congruential method
const unsigned long ARA = 1664525UL;
const unsigned long BRA = 1073904223UL;
static unsigned long  IDUM = time(NULL);// seed for semi_ra(),ra()
inline double semi_ra(){ //[0,1)
    IDUM=IDUM * ARA + BRA;
    return (double)(IDUM)/((double)ULONG_MAX+1.);
}
inline double ra(){ //[0,1]
    IDUM=IDUM * ARA + BRA;
    return (double)(IDUM)/ULONG_MAX;
}

/////////////////////// uniform random in [-1.0,1.0] with rand() 
const int GIVEN_SEED = 0;
const int RANDOM_SEED = 1;
inline void SRA(const int &SW_seed, const unsigned int &seed){
    if(SW_seed == RANDOM_SEED){
      srand(time(NULL));
      fprintf(stderr, "# rand_seed= time(NULL)\n");
   }else if(SW_seed == GIVEN_SEED){
      srand(seed);
      fprintf(stderr, "# rand_seed= %u\n",seed);
    }else {
      fprintf(stderr, "SRA(): invalid SW_seed.\n");
      exit_job(EXIT_FAILURE);
    }
}
inline double RA(){// uniform in [-1, 1]
    return (2.0*((double)rand() - RAND_MAX*0.5)/(RAND_MAX));
}
inline double RAx(const double &x){ // uniform in [0, x)
	return (double)rand()/RAND_MAX * x;
}
//random point inside unit circle
inline void RA_circle(double &a, double &b){
  int inside = 0;
  do{
    a = RA();
    b = RA();
    if(a*a + b*b <= 1){
      inside = 1;
    }
  }while(!inside);
}
//random point on unit circle
inline void RA_on_circle(double &a, double &b){
  RA_circle(a, b);
  double norm = 1.0/sqrt(a*a + b*b);
  a *= norm;
  b *= norm;
}
//random point on 3D/4D unit sphere
// Numerical Recipes, Edition 3, pg. 1130
// G. Marsaglia, Annals of Mathematical Statistics, 43(2), 645-646 (1972)
inline void RA_on_sphere3D(double &a, double &b, double &c){
  double u0, u1, sqnorm, dmy;

  RA_circle(u0,u1);
  sqnorm = u0*u0 + u1*u1;
  dmy = sqrt(1.0 - sqnorm);

  a = 2.0 * u0 * dmy;
  b = 2.0 * u1 * dmy;
  c = 1.0 - 2.0 * sqnorm;
}
inline void RA_on_sphere4D(double &a, double &b, double &c, double &d){
  double u0, u1, u2, u3, scale;

  RA_circle(u0,u1);
  RA_circle(u2,u3);
  scale = sqrt( (1.0 - (u0*u0 + u1*u1)) / (u2*u2 + u3*u3)  );
  a = u0;
  b = u1;
  c = u2 * scale;
  d = u3 * scale;
}
/////////////////////// macro for simple arithmetics
inline double POW6(const double x){
  double dmy = x*x*x;
  return dmy*dmy;
}
inline double POW3(const double x){
    return x*x*x;
}
inline double SQ(const double x){
    return x*x;
}
inline int SQ(const int x){
    return x*x;
}
inline double ABS(const double x){
    return x>=0?x:-x;
} 
inline double MIN(const double &a, const double &b){
    return a<b?a:b;
}
inline int MIN(const int &a, const int &b){
    return a<b?a:b;
}
inline int MAX(const int &a, const int &b){
    return a>b?a:b;
}
inline double MAX(const double &a, const double &b){
    return a>b?a:b;
}
inline double SGN(const double &a){
    return a>=0?1.:-1.;
}

inline bool equal_tol(const double &a, const double &b, const double &rtol, 
		      const double &atol = TOL_MP){
  if(a == b) return true;

  double diff = ABS(a-b);
  if(diff <= atol) return true;
  
  double eps = MAX(ABS(a), ABS(b)) * rtol;
  return (diff <= eps ? true : false);
}
inline bool equal_mp(const double &a, const double &b){
  if(a == b) return true;
  double eps = (MAX(ABS(a), ABS(b)) + 10.0)*EPSILON_MP;
  return (ABS(a - b) <= eps ? true: false);
}
inline bool greater_than_mp(const double &a, const double &b){
  double eps = (MAX(ABS(a), ABS(b)) + 10.0)*EPSILON_MP;
  return (a-b) > eps;
}
inline bool less_than_mp(const double &a, const double &b){
  double eps = (MAX(ABS(a), ABS(b)) + 10.0)*EPSILON_MP;
  return (b-a) > eps;
}
inline bool non_zero_mp(const double &a){
  return (a > TOL_MP) || (-TOL_MP > a);
}
inline bool zero_mp(const double &a){
  return !non_zero_mp(a);
}
inline bool positive_mp(const double &a){
  return a > TOL_MP;
}
inline bool negative_mp(const double &a){
  return -TOL_MP > a;
}

///
inline void exit_job(int status){
  exit(status);
}

inline FILE *filecheckopen(const char *fname, const char st[]){
  FILE *fp;
  
  if ((fp=fopen(fname,st))==NULL) {
    fprintf(stderr,"file open %s not succeeded\n",fname);
    exit_job(EXIT_FAILURE);
  }
  return(fp);
}

inline int file_check(const char *filename){
    FILE *ftest;
    if(NULL==(ftest=fopen(filename,"r"))){
	printf("File <%s> does not exist.\n",filename);
	printf("Program terminated\n");
	exit_job(EXIT_FAILURE);
	return false;
    }
    fclose(ftest);
    return true;
}

#endif 
