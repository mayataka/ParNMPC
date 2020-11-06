/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_Simu_Matlab_api.h
 *
 * MATLAB Coder version            : 5.0
 * C/C++ source code generated on  : 06-Nov-2020 17:34:51
 */

#ifndef _CODER_SIMU_MATLAB_API_H
#define _CODER_SIMU_MATLAB_API_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

#define MAX_THREADS                    omp_get_max_threads()

/* Function Declarations */
extern void MEXGlobalSyncInFunction(const emlrtStack *sp);
extern void MEXGlobalSyncOutFunction(boolean_T skipDirtyCheck);
extern void Simu_Matlab(void);
extern void Simu_Matlab_api(int32_T nlhs);
extern void Simu_Matlab_atexit(void);
extern void Simu_Matlab_initialize(void);
extern void Simu_Matlab_terminate(void);
extern void Simu_Matlab_xil_shutdown(void);
extern void Simu_Matlab_xil_terminate(void);

#endif

/*
 * File trailer for _coder_Simu_Matlab_api.h
 *
 * [EOF]
 */
