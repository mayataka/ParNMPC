//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: Simu_Matlab.h
//
// MATLAB Coder version            : 5.0
// C/C++ source code generated on  : 06-Nov-2020 17:34:51
//
#ifndef SIMU_MATLAB_H
#define SIMU_MATLAB_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "Simu_Matlab_types.h"

// Type Definitions
#include <stdio.h>

// Variable Declarations
extern omp_nest_lock_t emlrtNestLockGlobal;

#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
extern void Simu_Matlab();
extern void Simu_Matlab_initialize();
extern void Simu_Matlab_terminate();

#endif

//
// File trailer for Simu_Matlab.h
//
// [EOF]
//
