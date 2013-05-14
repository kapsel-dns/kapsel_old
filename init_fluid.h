/*!
  \file init_fluid.h
  \brief Initialize fluid velocity fields (header file)
  \author Y. Nakayama
  \date 2006/06/27
  \version 1.1
 */
#ifndef INIT_FLUID_H_
#define INIT_FLUID_H_

#include "variable.h"
#include "macro.h"
#include "operate_omega.h"
#include "input.h"

/*!
  \brief Initializes velocity field for fluid at rest
  \param[out] zeta reduced vorticity field
  \param[out] uk_dc zero-wavenumber Fourier transform of the fluid velocity
 */
void Init_zeta_k(double **zeta, double *uk_dc);

#endif
