/*!
  \file parameter_define.h
  \brief Define global system parameters for FFT routines
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
 */
#ifndef PARAMETER_DEFINE_H_
#define PARAMETER_DEFINE_H_

#define DIM 3
#define QDIM ((DIM*(DIM+1))/2-1)
/////// for SW_FFT
#define Ooura 0
#define RFFTW 1
#define MPI_RFFTW 2
#define IMKL_FFT 3

#endif
