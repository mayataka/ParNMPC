/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_Simu_Matlab_api.c
 *
 * MATLAB Coder version            : 5.0
 * C/C++ source code generated on  : 06-Nov-2020 17:34:51
 */

/* Include Files */
#include "_coder_Simu_Matlab_api.h"
#include "_coder_Simu_Matlab_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131594U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "Simu_Matlab",                       /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Definitions */

/*
 * Arguments    : const emlrtStack *sp
 * Return Type  : void
 */
void MEXGlobalSyncInFunction(const emlrtStack *sp)
{
  const mxArray *ParNMPCGlobalVariableMx;
  int32_T iv[1];
  static const char * sv[1] = { "ParNMPCGlobalVariable" };

  static const uint32_T uv[4] = { 1510074372U, 2288097353U, 3044386473U,
    4264122680U };

  /* Marshall in global variables */
  ParNMPCGlobalVariableMx = emlrtMexGetVariablePtr("global",
    "ParNMPCGlobalVariable");
  if (ParNMPCGlobalVariableMx != NULL) {
    iv[0] = 4;
    emlrtCheckArrayChecksumR2018b(sp, ParNMPCGlobalVariableMx, true, iv, sv, uv);
  }
}

/*
 * Arguments    : boolean_T skipDirtyCheck
 * Return Type  : void
 */
void MEXGlobalSyncOutFunction(boolean_T skipDirtyCheck)
{
  (void)skipDirtyCheck;

  /* Marshall out global variables */
}

/*
 * Arguments    : int32_T nlhs
 * Return Type  : void
 */
void Simu_Matlab_api(int32_T nlhs)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  (void)nlhs;
  st.tls = emlrtRootTLSGlobal;

  /* Marshall in global variables */
  MEXGlobalSyncInFunction(&st);

  /* Invoke the target function */
  Simu_Matlab();

  /* Marshall out global variables */
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void Simu_Matlab_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  Simu_Matlab_xil_terminate();
  Simu_Matlab_xil_shutdown();
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void Simu_Matlab_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void Simu_Matlab_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * File trailer for _coder_Simu_Matlab_api.c
 *
 * [EOF]
 */
