/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_Simu_Matlab_mex.c
 *
 * MATLAB Coder version            : 5.0
 * C/C++ source code generated on  : 06-Nov-2020 17:34:51
 */

/* Include Files */
#include "_coder_Simu_Matlab_mex.h"
#include "_coder_Simu_Matlab_api.h"

/* Function Declarations */
MEXFUNCTION_LINKAGE void Simu_Matlab_mexFunction(int32_T nlhs, int32_T nrhs);

/* Function Definitions */

/*
 * Arguments    : int32_T nlhs
 *                int32_T nrhs
 * Return Type  : void
 */
void Simu_Matlab_mexFunction(int32_T nlhs, int32_T nrhs)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 0) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 0, 4,
                        11, "Simu_Matlab");
  }

  if (nlhs > 0) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 11,
                        "Simu_Matlab");
  }

  /* Call the function. */
  Simu_Matlab_api(nlhs);
}

/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[]
 *                int32_T nrhs
 *                const mxArray *prhs[]
 * Return Type  : void
 */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  (void)plhs;
  (void)prhs;
  mexAtExit(&Simu_Matlab_atexit);

  /* Module initialization. */
  Simu_Matlab_initialize();

  /* Dispatch the entry-point. */
  Simu_Matlab_mexFunction(nlhs, nrhs);

  /* Module termination. */
  Simu_Matlab_terminate();
}

/*
 * Arguments    : void
 * Return Type  : emlrtCTX
 */
emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/*
 * File trailer for _coder_Simu_Matlab_mex.c
 *
 * [EOF]
 */
