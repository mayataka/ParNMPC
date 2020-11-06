//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: Simu_Matlab.cpp
//
// MATLAB Coder version            : 5.0
// C/C++ source code generated on  : 06-Nov-2020 17:34:51
//

// Include Files
#include "Simu_Matlab.h"
#include "iiwa14.h"
#include "omp.h"
#include "stdio.h"
#include <cmath>
#include <cstring>
#include <stdio.h>

// Type Definitions
#include <stdio.h>

// Type Definitions
struct cell_wrap_7
{
  char f1[3];
};

struct emxArray_real_T_14x6x4
{
  double data[336];
  int size[3];
};

struct emxArray_real_T_14x14x6x4
{
  double data[4704];
  int size[4];
};

struct emxArray_real_T_29x6x4
{
  double data[696];
  int size[3];
};

struct emxArray_real_T_8x6x4
{
  double data[192];
  int size[3];
};

// Variable Definitions
static emxArray_real_T_14x6x4 lambdaInit;
static emxArray_real_T_8x6x4 uInit;
static emxArray_real_T_14x6x4 xInit;
static boolean_T xInit_not_empty;
static emxArray_real_T_29x6x4 zInit;
static emxArray_real_T_14x14x6x4 LAMBDAInit;
static FILE * eml_openfiles[20];
static boolean_T eml_autoflush[20];
omp_nest_lock_t emlrtNestLockGlobal;
static boolean_T isInitialized_Simu_Matlab = false;

// Function Declarations
static void NMPC_Solve(const double x0[14], const double p[168], double
  solution_lambda[336], double solution_u[192], double solution_x[336], double
  solution_z[696], double solution_LAMBDA[4704], double *output_iterInit, double
  *output_cost, double *output_KKTError, double *c_output_timeElapsed_searchDire,
  double *output_timeElapsed_lineSearch, double *c_output_timeElapsed_KKTErrorCh,
  double *output_timeElapsed_total, double *output_rho, double *output_iterTotal,
  double *output_exitflag);
static void NMPC_Solve_SearchDirection(const double x0[14], const double p_data[],
  double rho, double lambda_data[], int lambda_size[3], int mu_size[3], double
  u_data[], int u_size[3], double x_data[], int x_size[3], double z_data[], int
  z_size[3], double LAMBDA_data[], double *KKTError_stateEquation, double
  *KKTError_C, double *KKTError_Hu, double *KKTError_costateEquation, double
  *costL, double *timeElapsed);
static void NMPC_Solve_init();
static void OCP_F_Fu_Fx(const double u[8], const double x[14], double parIdx,
  double F[14], double Fu[112], double Fx[196]);
static void OCP_G(const double u_data[], const double x_data[], double G_data[],
                  int G_size[2]);
static void OCP_GEN_fdt_fudt_fxdt(const double u[8], const double x[14], double
  parIdx, double fdt[14], double fudt[112], double fxdt[196]);
static void SIM_Plant_RK4(const double u[7], const double x[14], double xNext[14]);
static void b_inv(const double x[64], double y[64]);
static void c_fraction_to_boundary_parallel(const double u_i_data[], const
  double x_i_data[], double z_i_data[], const double u_k_i_data[], const double
  x_k_i_data[], const double z_k_i_data[], double rho, double *stepSizeMaxZ_i,
  double *stepSizeMaxG_i);
static int cfclose(double fid);
static signed char cfopen(const char * cfilename, const char * cpermission);
static void coarse_update_func(double lambda_i_data[], double u_i_data[], double
  x_i_data[], const double z_i_data[], const double p_i_data[], double
  xPrev_i_data[], double lambdaNext_i_data[], double LAMBDA_i_data[], double rho,
  double i, double p_muu_F_i_data[], int p_muu_F_i_size[3], double
  p_muu_Lambda_i_data[], int p_muu_Lambda_i_size[3], double
  p_lambda_Lambda_i_data[], int p_lambda_Lambda_i_size[3], double
  p_x_Lambda_i_data[], int p_x_Lambda_i_size[3], double p_x_F_i_data[], int
  p_x_F_i_size[3], double KKTxEquation_i_data[], int KKTxEquation_i_size[2],
  double KKTC_i_data[], int KKTC_i_size[2], double KKTHu_i_data[], int
  KKTHu_i_size[2], double KKTlambdaEquation_i_data[], int
  KKTlambdaEquation_i_size[2], double L_i_data[], int L_i_size[2], double
  LB_i_data[], int LB_i_size[2]);
static void f_fu_fx_Wrapper(const double u[8], const double x[14], double parIdx,
  double f[14], double fu[112], double fx[196]);
static void fileManager(double varargin_1, FILE * *f, boolean_T *a);
static signed char filedata();
static void filedata_init();
static void inv(const double x[196], double y[196]);
static void mldivide(const double A[196], double B[14]);
static double rt_roundd(double u);
static void xzgetrf(double A[196], int ipiv[14], int *info);

// Function Definitions

//
// Arguments    : const double x0[14]
//                const double p[168]
//                double solution_lambda[336]
//                double solution_u[192]
//                double solution_x[336]
//                double solution_z[696]
//                double solution_LAMBDA[4704]
//                double *output_iterInit
//                double *output_cost
//                double *output_KKTError
//                double *c_output_timeElapsed_searchDire
//                double *output_timeElapsed_lineSearch
//                double *c_output_timeElapsed_KKTErrorCh
//                double *output_timeElapsed_total
//                double *output_rho
//                double *output_iterTotal
//                double *output_exitflag
// Return Type  : void
//
static void NMPC_Solve(const double x0[14], const double p[168], double
  solution_lambda[336], double solution_u[192], double solution_x[336], double
  solution_z[696], double solution_LAMBDA[4704], double *output_iterInit, double
  *output_cost, double *output_KKTError, double *c_output_timeElapsed_searchDire,
  double *output_timeElapsed_lineSearch, double *c_output_timeElapsed_KKTErrorCh,
  double *output_timeElapsed_total, double *output_rho, double *output_iterTotal,
  double *output_exitflag)
{
  double tStart;
  int mode;
  int lambdaSplit_size[3];
  static const signed char iv[192] = { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1 };

  int muSplit_size[3];
  int uSplit_size[3];
  int xSplit_size[3];
  int zSplit_size[3];
  int iter;
  int b_iter;
  boolean_T exitg1;
  double KKTError_stateEquation;
  double lambda_L1Norm;
  double KKTError_Hu;
  double KKTError_costateEquation;
  double tSearchDirection;
  tStart = omp_get_wtime();

  //  only for the very first problem
  //  Reshape
  if (!xInit_not_empty) {
    lambdaInit.size[0] = 14;
    lambdaInit.size[1] = 6;
    lambdaInit.size[2] = 4;
    std::memset(&lambdaInit.data[0], 0, 336U * sizeof(double));
    uInit.size[0] = 8;
    uInit.size[1] = 6;
    uInit.size[2] = 4;
    for (mode = 0; mode < 192; mode++) {
      uInit.data[mode] = iv[mode];
    }

    xInit.size[0] = 14;
    xInit.size[1] = 6;
    xInit.size[2] = 4;
    std::memset(&xInit.data[0], 0, 336U * sizeof(double));
    xInit_not_empty = true;
    LAMBDAInit.size[0] = 14;
    LAMBDAInit.size[1] = 14;
    LAMBDAInit.size[2] = 6;
    LAMBDAInit.size[3] = 4;
    std::memset(&LAMBDAInit.data[0], 0, 4704U * sizeof(double));

    //  assert
  }

  zInit.size[0] = 29;
  zInit.size[1] = 6;
  zInit.size[2] = 4;
  for (mode = 0; mode < 696; mode++) {
    zInit.data[mode] = 1.0;
  }

  //  Init
  lambdaSplit_size[0] = 14;
  lambdaSplit_size[1] = 6;
  lambdaSplit_size[2] = 4;
  std::memcpy(&solution_lambda[0], &lambdaInit.data[0], 336U * sizeof(double));
  muSplit_size[0] = 0;
  muSplit_size[1] = 6;
  muSplit_size[2] = 4;
  uSplit_size[0] = 8;
  uSplit_size[1] = 6;
  uSplit_size[2] = 4;
  std::memcpy(&solution_u[0], &uInit.data[0], 192U * sizeof(double));
  xSplit_size[0] = 14;
  xSplit_size[1] = 6;
  xSplit_size[2] = 4;
  std::memcpy(&solution_x[0], &xInit.data[0], 336U * sizeof(double));
  zSplit_size[0] = 29;
  zSplit_size[1] = 6;
  zSplit_size[2] = 4;
  std::memcpy(&solution_z[0], &zInit.data[0], 696U * sizeof(double));
  std::memcpy(&solution_LAMBDA[0], &LAMBDAInit.data[0], 4704U * sizeof(double));
  *output_iterInit = 0.0;
  *c_output_timeElapsed_searchDire = 0.0;
  *output_timeElapsed_lineSearch = 0.0;
  *c_output_timeElapsed_KKTErrorCh = 0.0;
  *output_exitflag = 0.0;
  *output_cost = 0.0;
  *output_KKTError = 0.0;
  mode = 1;

  //  initialize a filter [LAll;xEq+C;flag]
  //  Iteration
  iter = 1;
  b_iter = 0;
  exitg1 = false;
  while ((!exitg1) && (b_iter < 10)) {
    iter = b_iter + 1;

    //  backup
    //         %% Search direction
    NMPC_Solve_SearchDirection(x0, p, 0.001, solution_lambda, lambdaSplit_size,
      muSplit_size, solution_u, uSplit_size, solution_x, xSplit_size, solution_z,
      zSplit_size, solution_LAMBDA, &KKTError_stateEquation, &lambda_L1Norm,
      &KKTError_Hu, &KKTError_costateEquation, output_cost, &tSearchDirection);
    *c_output_timeElapsed_searchDire += tSearchDirection;

    //         %% Line search
    //         %% KKT error
    //  avoid chattering of the num of iter
    lambda_L1Norm = 0.0;
    for (int k = 0; k < 336; k++) {
      lambda_L1Norm += std::abs(solution_lambda[k]);
    }

    lambda_L1Norm = lambda_L1Norm / 14.0 / 24.0;
    if (100.0 > lambda_L1Norm) {
      lambda_L1Norm = 100.0;
    }

    lambda_L1Norm /= 100.0;
    tSearchDirection = KKTError_Hu / lambda_L1Norm / 5.0;
    lambda_L1Norm = KKTError_costateEquation / lambda_L1Norm / 5.0;
    *output_KKTError = KKTError_stateEquation * 5.0;
    if (KKTError_stateEquation * 5.0 < 0.0) {
      *output_KKTError = 0.0;
    }

    if (*output_KKTError < tSearchDirection) {
      *output_KKTError = tSearchDirection;
    }

    if (*output_KKTError < lambda_L1Norm) {
      *output_KKTError = lambda_L1Norm;
    }

    //         %% print
    //         %% barrier parameter update
    if (mode == 1) {
      //  init
      if ((*output_KKTError < 0.001) || (b_iter + 1 >= 10)) {
        mode = 2;
        lambdaInit.size[0] = 14;
        lambdaInit.size[1] = 6;
        lambdaInit.size[2] = 4;
        std::memcpy(&lambdaInit.data[0], &solution_lambda[0], 336U * sizeof
                    (double));
        uInit.size[0] = 8;
        uInit.size[1] = 6;
        uInit.size[2] = 4;
        std::memcpy(&uInit.data[0], &solution_u[0], 192U * sizeof(double));
        xInit.size[0] = 14;
        xInit.size[1] = 6;
        xInit.size[2] = 4;
        std::memcpy(&xInit.data[0], &solution_x[0], 336U * sizeof(double));
        xInit_not_empty = true;
        zInit.size[0] = 29;
        zInit.size[1] = 6;
        zInit.size[2] = 4;
        std::memcpy(&zInit.data[0], &solution_z[0], 696U * sizeof(double));
        LAMBDAInit.size[0] = 14;
        LAMBDAInit.size[1] = 14;
        LAMBDAInit.size[2] = 6;
        LAMBDAInit.size[3] = 4;
        std::memcpy(&LAMBDAInit.data[0], &solution_LAMBDA[0], 4704U * sizeof
                    (double));
        *output_iterInit = static_cast<double>(b_iter) + 1.0;
        if (*output_KKTError < 0.001) {
          *output_exitflag = 1.0;
          exitg1 = true;
        } else {
          b_iter++;
        }
      } else {
        b_iter++;
      }
    } else {
      //  decay
      if (*output_KKTError < 0.001) {
        *output_exitflag = 1.0;
        exitg1 = true;
      } else {
        b_iter++;
      }
    }
  }

  //     %% Reshape
  *output_rho = 0.001;
  *output_iterTotal = iter;
  lambda_L1Norm = omp_get_wtime();
  *output_timeElapsed_total = lambda_L1Norm - tStart;
}

//
// Arguments    : const double x0[14]
//                const double p_data[]
//                double rho
//                double lambda_data[]
//                int lambda_size[3]
//                int mu_size[3]
//                double u_data[]
//                int u_size[3]
//                double x_data[]
//                int x_size[3]
//                double z_data[]
//                int z_size[3]
//                double LAMBDA_data[]
//                double *KKTError_stateEquation
//                double *KKTError_C
//                double *KKTError_Hu
//                double *KKTError_costateEquation
//                double *costL
//                double *timeElapsed
// Return Type  : void
//
static void NMPC_Solve_SearchDirection(const double x0[14], const double p_data[],
  double rho, double lambda_data[], int lambda_size[3], int mu_size[3], double
  u_data[], int u_size[3], double x_data[], int x_size[3], double z_data[], int
  z_size[3], double LAMBDA_data[], double *KKTError_stateEquation, double
  *KKTError_C, double *KKTError_Hu, double *KKTError_costateEquation, double
  *costL, double *timeElapsed)
{
  double tStart;
  double lambda_k_data[336];
  double u_k_data[192];
  double x_k_data[336];
  double z_k_data[696];
  double lambdaNext_data[336];
  double xPrev_data[336];
  double dlambda_data[336];
  int i;
  int b_i;
  int c_i;
  int xPrev_data_tmp;
  int b_xPrev_data_tmp;
  double L_data[24];
  double KKTxEquation_data[24];
  double KKTlambdaEquation_data[24];
  double KKTHu_data[24];
  int i1;
  double p_x_F_data[4704];
  double p_x_Lambda_data[4704];
  double p_lambda_Lambda_data[4704];
  double p_muu_Lambda_data[2688];
  double tEnd;
  double p_muu_F_data[2688];
  double LAMBDA_i_data[1176];
  double lambdaNext_i_data[84];
  double xPrev_i_data[84];
  double x_k_i_data[84];
  double dmu_u_j_i[8];
  double u_k_i_data[48];
  double dlambda_j_i[14];
  double lambda_next[14];
  double lambda_i_data[84];
  double LB_i_data[6];
  int LB_i_size[2];
  double L_i_data[6];
  double stepSizeG_data[4];
  int L_i_size[2];
  double stepSizeZ_data[4];
  double KKTlambdaEquation_i_data[6];
  double z_i_data[174];
  int KKTlambdaEquation_i_size[2];
  double stepSizeG_i;
  double KKTHu_i_data[6];
  double stepSizeZ_i;
  int KKTHu_i_size[2];
  double KKTC_i_data[6];
  int KKTC_i_size[2];
  int i2;
  double KKTxEquation_i_data[6];
  int KKTxEquation_i_size[2];
  int j;
  double p_x_F_i_data[1176];
  int p_x_F_i_size[3];
  double p_x_Lambda_i_data[1176];
  int p_x_Lambda_i_size[3];
  int i3;
  double p_lambda_Lambda_i_data[1176];
  int p_lambda_Lambda_i_size[3];
  double p_muu_Lambda_i_data[672];
  int p_muu_Lambda_i_size[3];
  double p_muu_F_i_data[672];
  int p_muu_F_i_size[3];
  double p_i_data[42];
  double b_u_k_data[48];
  double b_z_k_data[174];
  int lambda_i_data_tmp;
  tStart = omp_get_wtime();

  //  global variables
  //  parallel seg
  //  Code generation
  //  backup
  std::memcpy(&lambda_k_data[0], &lambda_data[0], 336U * sizeof(double));
  std::memcpy(&u_k_data[0], &u_data[0], 192U * sizeof(double));
  std::memcpy(&x_k_data[0], &x_data[0], 336U * sizeof(double));
  std::memcpy(&z_k_data[0], &z_data[0], 696U * sizeof(double));

  //  line search parameters
  //  local variables
  std::memset(&lambdaNext_data[0], 0, 336U * sizeof(double));
  std::memset(&xPrev_data[0], 0, 336U * sizeof(double));
  std::memset(&dlambda_data[0], 0, 336U * sizeof(double));

  //  coupling variable for each segment
  std::memcpy(&xPrev_data[0], &x0[0], 14U * sizeof(double));
  for (i = 0; i < 3; i++) {
    for (b_i = 0; b_i < 14; b_i++) {
      xPrev_data_tmp = b_i + 84 * (i + 1);
      b_xPrev_data_tmp = (b_i + 84 * i) + 70;
      xPrev_data[xPrev_data_tmp] = x_data[b_xPrev_data_tmp];
      lambdaNext_data[b_xPrev_data_tmp] = lambda_data[xPrev_data_tmp];
      for (i1 = 0; i1 < 14; i1++) {
        b_xPrev_data_tmp = i1 + 14 * b_i;
        LAMBDA_data[(b_xPrev_data_tmp + 1176 * i) + 980] =
          LAMBDA_data[b_xPrev_data_tmp + 1176 * (i + 1)];
      }
    }
  }

  for (b_i = 0; b_i < 14; b_i++) {
    std::memset(&LAMBDA_data[b_i * 14 + 4508], 0, 14U * sizeof(double));
  }

#pragma omp parallel for \
 num_threads(4 > omp_get_max_threads() ? omp_get_max_threads() : 4) \
 private(LAMBDA_i_data,lambdaNext_i_data,xPrev_i_data,x_k_i_data,u_k_i_data,lambda_i_data,LB_i_data,LB_i_size,L_i_data,L_i_size,KKTlambdaEquation_i_data,KKTlambdaEquation_i_size,KKTHu_i_data,KKTHu_i_size,KKTC_i_data,KKTC_i_size,KKTxEquation_i_data,KKTxEquation_i_size,p_x_F_i_data,p_x_F_i_size,p_x_Lambda_i_data,p_x_Lambda_i_size,p_lambda_Lambda_i_data,p_lambda_Lambda_i_size,p_muu_Lambda_i_data,p_muu_Lambda_i_size,p_muu_F_i_data,p_muu_F_i_size,p_i_data,z_i_data,i2,i3,lambda_i_data_tmp,j)

  for (c_i = 0; c_i < 4; c_i++) {
    //     %% V(:,index_inside_sig_j,which_sigment_i)
    //      for i=1:1:DoP
    for (i2 = 0; i2 < 6; i2++) {
      std::memcpy(&lambda_i_data[i2 * 14], &lambda_data[c_i * 84 + i2 * 14], 14U
                  * sizeof(double));
      std::memcpy(&u_k_i_data[i2 * 8], &u_data[c_i * 48 + i2 * 8], 8U * sizeof
                  (double));
      for (i3 = 0; i3 < 14; i3++) {
        lambda_i_data_tmp = i3 + 14 * i2;
        j = lambda_i_data_tmp + 84 * c_i;
        x_k_i_data[lambda_i_data_tmp] = x_data[j];
        xPrev_i_data[lambda_i_data_tmp] = xPrev_data[j];
        lambdaNext_i_data[lambda_i_data_tmp] = lambdaNext_data[j];
        std::memcpy(&LAMBDA_i_data[i2 * 196 + i3 * 14], &LAMBDA_data[(c_i * 1176
          + i2 * 196) + i3 * 14], 14U * sizeof(double));
      }

      std::memcpy(&z_i_data[i2 * 29], &z_data[c_i * 174 + i2 * 29], 29U * sizeof
                  (double));
      for (i3 = 0; i3 < 7; i3++) {
        lambda_i_data_tmp = i3 + 7 * i2;
        p_i_data[lambda_i_data_tmp] = p_data[lambda_i_data_tmp + 42 * c_i];
      }
    }

    coarse_update_func(lambda_i_data, u_k_i_data, x_k_i_data, z_i_data, p_i_data,
                       xPrev_i_data, lambdaNext_i_data, LAMBDA_i_data, rho,
                       static_cast<double>(c_i) + 1.0, p_muu_F_i_data,
                       p_muu_F_i_size, p_muu_Lambda_i_data, p_muu_Lambda_i_size,
                       p_lambda_Lambda_i_data, p_lambda_Lambda_i_size,
                       p_x_Lambda_i_data, p_x_Lambda_i_size, p_x_F_i_data,
                       p_x_F_i_size, KKTxEquation_i_data, KKTxEquation_i_size,
                       KKTC_i_data, KKTC_i_size, KKTHu_i_data, KKTHu_i_size,
                       KKTlambdaEquation_i_data, KKTlambdaEquation_i_size,
                       L_i_data, L_i_size, LB_i_data, LB_i_size);

    //  Recover
    //
    for (i2 = 0; i2 < 6; i2++) {
      std::memcpy(&lambda_data[c_i * 84 + i2 * 14], &lambda_i_data[i2 * 14], 14U
                  * sizeof(double));
      std::memcpy(&u_data[c_i * 48 + i2 * 8], &u_k_i_data[i2 * 8], 8U * sizeof
                  (double));
      for (i3 = 0; i3 < 14; i3++) {
        lambda_i_data_tmp = i3 + 14 * i2;
        j = lambda_i_data_tmp + 84 * c_i;
        x_data[j] = x_k_i_data[lambda_i_data_tmp];
        xPrev_data[j] = xPrev_i_data[lambda_i_data_tmp];
        std::memcpy(&LAMBDA_data[(c_i * 1176 + i2 * 196) + i3 * 14],
                    &LAMBDA_i_data[i2 * 196 + i3 * 14], 14U * sizeof(double));
        lambdaNext_data[j] = lambdaNext_i_data[lambda_i_data_tmp];
        std::memcpy(&p_muu_F_data[(c_i * 672 + i2 * 112) + i3 * 8],
                    &p_muu_F_i_data[i2 * 112 + i3 * 8], 8U * sizeof(double));
        std::memcpy(&p_muu_Lambda_data[(c_i * 672 + i2 * 112) + i3 * 8],
                    &p_muu_Lambda_i_data[i2 * 112 + i3 * 8], 8U * sizeof(double));
        std::memcpy(&p_lambda_Lambda_data[(c_i * 1176 + i2 * 196) + i3 * 14],
                    &p_lambda_Lambda_i_data[i2 * 196 + i3 * 14], 14U * sizeof
                    (double));
        std::memcpy(&p_x_Lambda_data[(c_i * 1176 + i2 * 196) + i3 * 14],
                    &p_x_Lambda_i_data[i2 * 196 + i3 * 14], 14U * sizeof(double));
        std::memcpy(&p_x_F_data[(c_i * 1176 + i2 * 196) + i3 * 14],
                    &p_x_F_i_data[i2 * 196 + i3 * 14], 14U * sizeof(double));
      }

      lambda_i_data_tmp = i2 + 6 * c_i;
      KKTxEquation_data[lambda_i_data_tmp] = KKTxEquation_i_data[i2];
      KKTHu_data[lambda_i_data_tmp] = KKTHu_i_data[i2];
      KKTlambdaEquation_data[lambda_i_data_tmp] = KKTlambdaEquation_i_data[i2];
      L_data[lambda_i_data_tmp] = L_i_data[i2];
    }
  }

  //     %%
  *KKTError_stateEquation = KKTxEquation_data[0];
  *KKTError_C = 0.0;
  *KKTError_Hu = KKTHu_data[0];
  *KKTError_costateEquation = KKTlambdaEquation_data[0];
  *costL = L_data[0];
  for (xPrev_data_tmp = 0; xPrev_data_tmp < 23; xPrev_data_tmp++) {
    tEnd = KKTxEquation_data[xPrev_data_tmp + 1];
    if (*KKTError_stateEquation < tEnd) {
      *KKTError_stateEquation = tEnd;
    }

    tEnd = KKTHu_data[xPrev_data_tmp + 1];
    if (*KKTError_Hu < tEnd) {
      *KKTError_Hu = tEnd;
    }

    tEnd = KKTlambdaEquation_data[xPrev_data_tmp + 1];
    if (*KKTError_costateEquation < tEnd) {
      *KKTError_costateEquation = tEnd;
    }

    *costL += L_data[xPrev_data_tmp + 1];
  }

  //     %% Step 2: Backward correction due to the approximation of lambda
  for (i = 0; i < 3; i++) {
    for (b_xPrev_data_tmp = 0; b_xPrev_data_tmp < 6; b_xPrev_data_tmp++) {
      if (6 - b_xPrev_data_tmp == 6) {
        std::memcpy(&lambda_next[0], &lambda_data[i * -84 + 252], 14U * sizeof
                    (double));
      } else {
        std::memcpy(&lambda_next[0], &lambda_data[(i * -84 + b_xPrev_data_tmp *
          -14) + 252], 14U * sizeof(double));
      }

      for (b_i = 0; b_i < 14; b_i++) {
        xPrev_data_tmp = (b_i + 14 * (5 - b_xPrev_data_tmp)) + 84 * (2 - i);
        dlambda_data[xPrev_data_tmp] = lambda_next[b_i] -
          lambdaNext_data[xPrev_data_tmp];
      }

      for (b_i = 0; b_i < 14; b_i++) {
        tEnd = 0.0;
        for (i1 = 0; i1 < 14; i1++) {
          tEnd += p_lambda_Lambda_data[((b_i + 14 * i1) + 196 * (5 -
            b_xPrev_data_tmp)) + 1176 * (2 - i)] * dlambda_data[(i1 + 14 * (5 -
            b_xPrev_data_tmp)) + 84 * (2 - i)];
        }

        i1 = (b_i + 14 * (5 - b_xPrev_data_tmp)) + 84 * (2 - i);
        lambda_data[i1] -= tEnd;
      }
    }
  }

#pragma omp parallel for \
 num_threads(4 > omp_get_max_threads() ? omp_get_max_threads() : 4) \
 private(dmu_u_j_i,dlambda_j_i,x_k_i_data,u_k_i_data,i2,j,stepSizeG_i,i3)

  for (c_i = 0; c_i < 4; c_i++) {
    //      for i=1:1:DoP
    for (i2 = 0; i2 < 6; i2++) {
      std::memcpy(&u_k_i_data[i2 * 8], &u_data[c_i * 48 + i2 * 8], 8U * sizeof
                  (double));
      std::memcpy(&x_k_i_data[i2 * 14], &x_data[c_i * 84 + i2 * 14], 14U *
                  sizeof(double));
    }

    for (j = 0; j < 6; j++) {
      for (i2 = 0; i2 < 8; i2++) {
        stepSizeG_i = 0.0;
        for (i3 = 0; i3 < 14; i3++) {
          stepSizeG_i += p_muu_Lambda_data[((i2 + 8 * i3) + 112 * (5 - j)) + 672
            * c_i] * dlambda_data[(i3 + 14 * (5 - j)) + 84 * c_i];
        }

        dmu_u_j_i[i2] = stepSizeG_i;
      }

      for (i2 = 0; i2 < 14; i2++) {
        stepSizeG_i = 0.0;
        for (i3 = 0; i3 < 14; i3++) {
          stepSizeG_i += p_x_Lambda_data[((i2 + 14 * i3) + 196 * (5 - j)) + 1176
            * c_i] * dlambda_data[(i3 + 14 * (5 - j)) + 84 * c_i];
        }

        dlambda_j_i[i2] = stepSizeG_i;
      }

      for (i2 = 0; i2 < 8; i2++) {
        i3 = i2 + 8 * (5 - j);
        stepSizeG_i = u_k_i_data[i3] - dmu_u_j_i[i2];
        dmu_u_j_i[i2] = stepSizeG_i;
        u_k_i_data[i3] = stepSizeG_i;
      }

      for (i2 = 0; i2 < 14; i2++) {
        i3 = i2 + 14 * (5 - j);
        stepSizeG_i = x_k_i_data[i3] - dlambda_j_i[i2];
        dlambda_j_i[i2] = stepSizeG_i;
        x_k_i_data[i3] = stepSizeG_i;
      }
    }

    for (i2 = 0; i2 < 6; i2++) {
      std::memcpy(&u_data[c_i * 48 + i2 * 8], &u_k_i_data[i2 * 8], 8U * sizeof
                  (double));
      std::memcpy(&x_data[c_i * 84 + i2 * 14], &x_k_i_data[i2 * 14], 14U *
                  sizeof(double));
    }
  }

  //     %% Step 3: Forward correction due to the approximation of x
  for (i = 0; i < 4; i++) {
    for (b_xPrev_data_tmp = 0; b_xPrev_data_tmp < 6; b_xPrev_data_tmp++) {
      if (b_xPrev_data_tmp + 1 == 1) {
        if (i + 1 == 1) {
          std::memcpy(&lambda_next[0], &x0[0], 14U * sizeof(double));
        } else {
          std::memcpy(&lambda_next[0], &x_data[i * 84 + -14], 14U * sizeof
                      (double));
        }
      } else {
        std::memcpy(&lambda_next[0], &x_data[(i * 84 + b_xPrev_data_tmp * 14) +
                    -14], 14U * sizeof(double));
      }

      for (b_i = 0; b_i < 14; b_i++) {
        xPrev_data_tmp = (b_i + 14 * b_xPrev_data_tmp) + 84 * i;
        lambdaNext_data[xPrev_data_tmp] = lambda_next[b_i] -
          xPrev_data[xPrev_data_tmp];
      }

      for (b_i = 0; b_i < 14; b_i++) {
        tEnd = 0.0;
        for (i1 = 0; i1 < 14; i1++) {
          tEnd += p_x_F_data[((b_i + 14 * i1) + 196 * b_xPrev_data_tmp) + 1176 *
            i] * lambdaNext_data[(i1 + 14 * b_xPrev_data_tmp) + 84 * i];
        }

        i1 = (b_i + 14 * b_xPrev_data_tmp) + 84 * i;
        x_data[i1] -= tEnd;
      }
    }
  }

#pragma omp parallel for \
 num_threads(4 > omp_get_max_threads() ? omp_get_max_threads() : 4) \
 private(z_i_data,stepSizeG_i,stepSizeZ_i,dmu_u_j_i,x_k_i_data,u_k_i_data,lambda_i_data,i2,j,i3,lambda_i_data_tmp,b_u_k_data,b_z_k_data) \
 firstprivate(lambdaNext_i_data)

  for (c_i = 0; c_i < 4; c_i++) {
    //      for i=1:1:DoP
    for (i2 = 0; i2 < 6; i2++) {
      std::memcpy(&lambda_i_data[i2 * 14], &lambda_data[c_i * 84 + i2 * 14], 14U
                  * sizeof(double));
      std::memcpy(&u_k_i_data[i2 * 8], &u_data[c_i * 48 + i2 * 8], 8U * sizeof
                  (double));
    }

    //  line search variables
    for (j = 0; j < 6; j++) {
      for (i2 = 0; i2 < 8; i2++) {
        stepSizeG_i = 0.0;
        for (i3 = 0; i3 < 14; i3++) {
          stepSizeG_i += p_muu_F_data[((i2 + 8 * i3) + 112 * (5 - j)) + 672 *
            c_i] * lambdaNext_data[(i3 + 14 * (5 - j)) + 84 * c_i];
        }

        dmu_u_j_i[i2] = stepSizeG_i;
      }

      for (i2 = 0; i2 < 14; i2++) {
        stepSizeG_i = 0.0;
        for (i3 = 0; i3 < 14; i3++) {
          stepSizeG_i += LAMBDA_data[((i2 + 14 * i3) + 196 * (5 - j)) + 1176 *
            c_i] * lambdaNext_data[(i3 + 14 * (5 - j)) + 84 * c_i];
        }

        lambda_i_data_tmp = i2 + 14 * (5 - j);
        lambda_i_data[lambda_i_data_tmp] -= stepSizeG_i;
      }

      for (i2 = 0; i2 < 8; i2++) {
        i3 = i2 + 8 * (5 - j);
        stepSizeG_i = u_k_i_data[i3] - dmu_u_j_i[i2];
        dmu_u_j_i[i2] = stepSizeG_i;
        u_k_i_data[i3] = stepSizeG_i;
      }
    }

    for (i2 = 0; i2 < 6; i2++) {
      std::memcpy(&lambda_data[c_i * 84 + i2 * 14], &lambda_i_data[i2 * 14], 14U
                  * sizeof(double));
      std::memcpy(&u_data[c_i * 48 + i2 * 8], &u_k_i_data[i2 * 8], 8U * sizeof
                  (double));
      std::memcpy(&z_i_data[i2 * 29], &z_data[c_i * 174 + i2 * 29], 29U * sizeof
                  (double));
      std::memcpy(&x_k_i_data[i2 * 14], &x_data[c_i * 84 + i2 * 14], 14U *
                  sizeof(double));
      std::memcpy(&b_u_k_data[i2 * 8], &u_k_data[c_i * 48 + i2 * 8], 8U * sizeof
                  (double));
      std::memcpy(&lambdaNext_i_data[i2 * 14], &x_k_data[c_i * 84 + i2 * 14],
                  14U * sizeof(double));
      std::memcpy(&b_z_k_data[i2 * 29], &z_k_data[c_i * 174 + i2 * 29], 29U *
                  sizeof(double));
    }

    c_fraction_to_boundary_parallel(u_k_i_data, x_k_i_data, z_i_data, b_u_k_data,
      lambdaNext_i_data, b_z_k_data, rho, &stepSizeZ_i, &stepSizeG_i);

    //  Recover
    for (i2 = 0; i2 < 6; i2++) {
      std::memcpy(&lambda_data[c_i * 84 + i2 * 14], &lambda_i_data[i2 * 14], 14U
                  * sizeof(double));
      std::memcpy(&u_data[c_i * 48 + i2 * 8], &u_k_i_data[i2 * 8], 8U * sizeof
                  (double));
      std::memcpy(&z_data[c_i * 174 + i2 * 29], &z_i_data[i2 * 29], 29U * sizeof
                  (double));
    }

    stepSizeZ_data[c_i] = stepSizeZ_i;
    stepSizeG_data[c_i] = stepSizeG_i;
  }

  //     %% Line Search to Guarantee Primal Stability
  tEnd = stepSizeG_data[0];
  if (stepSizeG_data[0] > stepSizeG_data[1]) {
    tEnd = stepSizeG_data[1];
  }

  if (tEnd > stepSizeG_data[2]) {
    tEnd = stepSizeG_data[2];
  }

  if (tEnd > stepSizeG_data[3]) {
    tEnd = stepSizeG_data[3];
  }

  if (tEnd != 1.0) {
    lambda_size[0] = 14;
    lambda_size[1] = 6;
    lambda_size[2] = 4;
    mu_size[0] = 0;
    mu_size[1] = 6;
    mu_size[2] = 4;
    u_size[0] = 8;
    u_size[1] = 6;
    u_size[2] = 4;
    x_size[0] = 14;
    x_size[1] = 6;
    x_size[2] = 4;
    for (b_i = 0; b_i < 4; b_i++) {
      for (i1 = 0; i1 < 6; i1++) {
        for (b_xPrev_data_tmp = 0; b_xPrev_data_tmp < 14; b_xPrev_data_tmp++) {
          xPrev_data_tmp = (b_xPrev_data_tmp + 14 * i1) + 84 * b_i;
          lambda_data[xPrev_data_tmp] = (1.0 - tEnd) *
            lambda_k_data[xPrev_data_tmp] + tEnd * lambda_data[xPrev_data_tmp];
        }

        for (b_xPrev_data_tmp = 0; b_xPrev_data_tmp < 8; b_xPrev_data_tmp++) {
          xPrev_data_tmp = (b_xPrev_data_tmp + 8 * i1) + 48 * b_i;
          u_data[xPrev_data_tmp] = (1.0 - tEnd) * u_k_data[xPrev_data_tmp] +
            tEnd * u_data[xPrev_data_tmp];
        }

        for (b_xPrev_data_tmp = 0; b_xPrev_data_tmp < 14; b_xPrev_data_tmp++) {
          xPrev_data_tmp = (b_xPrev_data_tmp + 14 * i1) + 84 * b_i;
          x_data[xPrev_data_tmp] = (1.0 - tEnd) * x_k_data[xPrev_data_tmp] +
            tEnd * x_data[xPrev_data_tmp];
        }
      }
    }
  }

  tEnd = stepSizeZ_data[0];
  if (stepSizeZ_data[0] > stepSizeZ_data[1]) {
    tEnd = stepSizeZ_data[1];
  }

  if (tEnd > stepSizeZ_data[2]) {
    tEnd = stepSizeZ_data[2];
  }

  if (tEnd > stepSizeZ_data[3]) {
    tEnd = stepSizeZ_data[3];
  }

  if (tEnd != 1.0) {
    z_size[0] = 29;
    z_size[1] = 6;
    z_size[2] = 4;
    for (b_i = 0; b_i < 4; b_i++) {
      for (i1 = 0; i1 < 6; i1++) {
        for (b_xPrev_data_tmp = 0; b_xPrev_data_tmp < 29; b_xPrev_data_tmp++) {
          xPrev_data_tmp = (b_xPrev_data_tmp + 29 * i1) + 174 * b_i;
          z_data[xPrev_data_tmp] = (1.0 - tEnd) * z_k_data[xPrev_data_tmp] +
            tEnd * z_data[xPrev_data_tmp];
        }
      }
    }
  }

  //     %%
  tEnd = omp_get_wtime();
  *timeElapsed = tEnd - tStart;
}

//
// Arguments    : void
// Return Type  : void
//
static void NMPC_Solve_init()
{
  xInit_not_empty = false;
}

//
// Arguments    : const double u[8]
//                const double x[14]
//                double parIdx
//                double F[14]
//                double Fu[112]
//                double Fx[196]
// Return Type  : void
//
static void OCP_F_Fu_Fx(const double u[8], const double x[14], double parIdx,
  double F[14], double Fu[112], double Fx[196])
{
  signed char Ix[196];
  int k;
  std::memset(&Ix[0], 0, 196U * sizeof(signed char));
  for (k = 0; k < 14; k++) {
    Ix[k + 14 * k] = 1;
  }

  //  M disabled
  //  'Euler'
  OCP_GEN_fdt_fudt_fxdt(u, x, parIdx, F, Fu, Fx);
  for (k = 0; k < 14; k++) {
    F[k] -= x[k];
  }

  for (k = 0; k < 196; k++) {
    Fx[k] -= static_cast<double>(Ix[k]);
  }
}

//
// Arguments    : const double u_data[]
//                const double x_data[]
//                double G_data[]
//                int G_size[2]
// Return Type  : void
//
static void OCP_G(const double u_data[], const double x_data[], double G_data[],
                  int G_size[2])
{
  // OCP_GEN_G
  //     G = OCP_GEN_G(IN1,IN2,IN3)
  //     This function was generated by the Symbolic Math Toolbox version 8.5.
  //     06-Nov-2020 17:34:34
  G_size[0] = 29;
  G_size[1] = 6;
  for (int i = 0; i < 6; i++) {
    double G_data_tmp;
    int b_G_data_tmp;
    int c_G_data_tmp;
    int d_G_data_tmp;
    int e_G_data_tmp;
    int f_G_data_tmp;
    int g_G_data_tmp;
    double d;
    double d1;
    double d2;
    double d3;
    double d4;
    double d5;
    G_data_tmp = u_data[8 * i];
    G_data[29 * i] = G_data_tmp + 10.0;
    b_G_data_tmp = 8 * i + 1;
    G_data[29 * i + 1] = u_data[b_G_data_tmp] + 10.0;
    c_G_data_tmp = 8 * i + 2;
    G_data[29 * i + 2] = u_data[c_G_data_tmp] + 10.0;
    d_G_data_tmp = 8 * i + 3;
    G_data[29 * i + 3] = u_data[d_G_data_tmp] + 10.0;
    e_G_data_tmp = 8 * i + 4;
    G_data[29 * i + 4] = u_data[e_G_data_tmp] + 10.0;
    f_G_data_tmp = 8 * i + 5;
    G_data[29 * i + 5] = u_data[f_G_data_tmp] + 10.0;
    g_G_data_tmp = 8 * i + 6;
    G_data[29 * i + 6] = u_data[g_G_data_tmp] + 10.0;
    G_data[29 * i + 7] = -G_data_tmp + 10.0;
    G_data[29 * i + 8] = -u_data[b_G_data_tmp] + 10.0;
    G_data[29 * i + 9] = -u_data[c_G_data_tmp] + 10.0;
    G_data[29 * i + 10] = -u_data[d_G_data_tmp] + 10.0;
    G_data[29 * i + 11] = -u_data[e_G_data_tmp] + 10.0;
    G_data[29 * i + 12] = -u_data[f_G_data_tmp] + 10.0;
    G_data[29 * i + 13] = -u_data[g_G_data_tmp] + 10.0;
    b_G_data_tmp = 8 * i + 7;
    G_data[29 * i + 14] = u_data[b_G_data_tmp];
    G_data_tmp = x_data[14 * i + 7];
    G_data[29 * i + 15] = (u_data[b_G_data_tmp] + 1.5707963267948966) +
      G_data_tmp;
    d = x_data[14 * i + 8];
    G_data[29 * i + 16] = (u_data[b_G_data_tmp] + 1.5707963267948966) + d;
    d1 = x_data[14 * i + 9];
    G_data[29 * i + 17] = (u_data[b_G_data_tmp] + 1.5707963267948966) + d1;
    d2 = x_data[14 * i + 10];
    G_data[29 * i + 18] = (u_data[b_G_data_tmp] + 1.5707963267948966) + d2;
    d3 = x_data[14 * i + 11];
    G_data[29 * i + 19] = (u_data[b_G_data_tmp] + 1.5707963267948966) + d3;
    d4 = x_data[14 * i + 12];
    G_data[29 * i + 20] = (u_data[b_G_data_tmp] + 1.5707963267948966) + d4;
    d5 = x_data[14 * i + 13];
    G_data[29 * i + 21] = (u_data[b_G_data_tmp] + 1.5707963267948966) + d5;
    G_data[29 * i + 22] = (u_data[b_G_data_tmp] + 1.5707963267948966) -
      G_data_tmp;
    G_data[29 * i + 23] = (u_data[b_G_data_tmp] + 1.5707963267948966) - d;
    G_data[29 * i + 24] = (u_data[b_G_data_tmp] + 1.5707963267948966) - d1;
    G_data[29 * i + 25] = (u_data[b_G_data_tmp] + 1.5707963267948966) - d2;
    G_data[29 * i + 26] = (u_data[b_G_data_tmp] + 1.5707963267948966) - d3;
    G_data[29 * i + 27] = (u_data[b_G_data_tmp] + 1.5707963267948966) - d4;
    G_data[29 * i + 28] = (u_data[b_G_data_tmp] + 1.5707963267948966) - d5;
  }
}

//
// Arguments    : const double u[8]
//                const double x[14]
//                double parIdx
//                double fdt[14]
//                double fudt[112]
//                double fxdt[196]
// Return Type  : void
//
static void OCP_GEN_fdt_fudt_fxdt(const double u[8], const double x[14], double
  parIdx, double fdt[14], double fudt[112], double fxdt[196])
{
  int i;
  f_fu_fx_Wrapper(u, x, parIdx, fdt, fudt, fxdt);
  for (i = 0; i < 14; i++) {
    fdt[i] *= 0.041667;
  }

  for (i = 0; i < 112; i++) {
    fudt[i] *= 0.041667;
  }

  for (i = 0; i < 196; i++) {
    fxdt[i] *= 0.041667;
  }
}

//
// Arguments    : const double u[7]
//                const double x[14]
//                double xNext[14]
// Return Type  : void
//
static void SIM_Plant_RK4(const double u[7], const double x[14], double xNext[14])
{
  int i;
  double q[7];
  double qd[7];
  double qdd[7];
  double tau[7];
  static const double dv[196] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };

  double d;
  double b_x[14];
  double k2[14];
  double k3[14];

  //  dx = f(u,x,p)
  //  Specify your own f(u,x,p) function for code generation
  for (i = 0; i < 7; i++) {
    q[i] = x[i];
    qd[i] = x[i + 7];
    qdd[i] = 0.0;
    tau[i] = u[i];
  }

  sim_qdd_cal(q, qd, qdd, tau);
  for (i = 0; i < 7; i++) {
    xNext[i] = qd[i];
    xNext[i + 7] = qdd[i];
  }

  mldivide(dv, xNext);
  for (i = 0; i < 14; i++) {
    d = 0.001 * xNext[i];
    xNext[i] = d;
    b_x[i] = x[i] + d / 2.0;
  }

  //  dx = f(u,x,p)
  //  Specify your own f(u,x,p) function for code generation
  for (i = 0; i < 7; i++) {
    q[i] = b_x[i];
    qd[i] = b_x[i + 7];
    qdd[i] = 0.0;
    tau[i] = u[i];
  }

  sim_qdd_cal(q, qd, qdd, tau);
  for (i = 0; i < 7; i++) {
    k2[i] = qd[i];
    k2[i + 7] = qdd[i];
  }

  mldivide(dv, k2);
  for (i = 0; i < 14; i++) {
    d = 0.001 * k2[i];
    k2[i] = d;
    b_x[i] = x[i] + d / 2.0;
  }

  //  dx = f(u,x,p)
  //  Specify your own f(u,x,p) function for code generation
  for (i = 0; i < 7; i++) {
    q[i] = b_x[i];
    qd[i] = b_x[i + 7];
    qdd[i] = 0.0;
    tau[i] = u[i];
  }

  sim_qdd_cal(q, qd, qdd, tau);
  for (i = 0; i < 7; i++) {
    k3[i] = qd[i];
    k3[i + 7] = qdd[i];
  }

  mldivide(dv, k3);
  for (i = 0; i < 14; i++) {
    d = 0.001 * k3[i];
    k3[i] = d;
    b_x[i] = x[i] + d;
  }

  //  dx = f(u,x,p)
  //  Specify your own f(u,x,p) function for code generation
  for (i = 0; i < 7; i++) {
    q[i] = b_x[i];
    qd[i] = b_x[i + 7];
    qdd[i] = 0.0;
    tau[i] = u[i];
  }

  sim_qdd_cal(q, qd, qdd, tau);
  for (i = 0; i < 7; i++) {
    b_x[i] = qd[i];
    b_x[i + 7] = qdd[i];
  }

  mldivide(dv, b_x);
  for (i = 0; i < 14; i++) {
    xNext[i] = x[i] + (((xNext[i] + 2.0 * k2[i]) + 2.0 * k3[i]) + 0.001 * b_x[i])
      / 6.0;
  }
}

//
// Arguments    : const double x[64]
//                double y[64]
// Return Type  : void
//
static void b_inv(const double x[64], double y[64])
{
  int i;
  double b_x[64];
  int j;
  signed char ipiv[8];
  int b;
  int k;
  signed char p[8];
  int jy;
  int iy;
  int ix;
  int i1;
  for (i = 0; i < 64; i++) {
    y[i] = 0.0;
    b_x[i] = x[i];
  }

  for (i = 0; i < 8; i++) {
    ipiv[i] = static_cast<signed char>(i + 1);
  }

  for (j = 0; j < 7; j++) {
    int mmj_tmp;
    int jj;
    int jp1j;
    double smax;
    mmj_tmp = 6 - j;
    b = j * 9;
    jj = j * 9;
    jp1j = b + 2;
    jy = 8 - j;
    iy = 0;
    ix = b;
    smax = std::abs(b_x[jj]);
    for (k = 2; k <= jy; k++) {
      double s;
      ix++;
      s = std::abs(b_x[ix]);
      if (s > smax) {
        iy = k - 1;
        smax = s;
      }
    }

    if (b_x[jj + iy] != 0.0) {
      if (iy != 0) {
        iy += j;
        ipiv[j] = static_cast<signed char>(iy + 1);
        ix = j;
        for (k = 0; k < 8; k++) {
          smax = b_x[ix];
          b_x[ix] = b_x[iy];
          b_x[iy] = smax;
          ix += 8;
          iy += 8;
        }
      }

      i = (jj - j) + 8;
      for (ix = jp1j; ix <= i; ix++) {
        b_x[ix - 1] /= b_x[jj];
      }
    }

    jy = b + 8;
    iy = jj;
    for (b = 0; b <= mmj_tmp; b++) {
      smax = b_x[jy];
      if (b_x[jy] != 0.0) {
        ix = jj + 1;
        i = iy + 10;
        i1 = (iy - j) + 16;
        for (jp1j = i; jp1j <= i1; jp1j++) {
          b_x[jp1j - 1] += b_x[ix] * -smax;
          ix++;
        }
      }

      jy += 8;
      iy += 8;
    }
  }

  for (i = 0; i < 8; i++) {
    p[i] = static_cast<signed char>(i + 1);
  }

  for (k = 0; k < 7; k++) {
    if (ipiv[k] > k + 1) {
      jy = ipiv[k] - 1;
      iy = p[jy];
      p[jy] = p[k];
      p[k] = static_cast<signed char>(iy);
    }
  }

  for (k = 0; k < 8; k++) {
    b = (p[k] - 1) << 3;
    y[k + b] = 1.0;
    for (j = k + 1; j < 9; j++) {
      i = (j + b) - 1;
      if (y[i] != 0.0) {
        i1 = j + 1;
        for (ix = i1; ix < 9; ix++) {
          jy = (ix + b) - 1;
          y[jy] -= y[i] * b_x[(ix + ((j - 1) << 3)) - 1];
        }
      }
    }
  }

  for (j = 0; j < 8; j++) {
    jy = j << 3;
    for (k = 7; k >= 0; k--) {
      iy = k << 3;
      i = k + jy;
      if (y[i] != 0.0) {
        y[i] /= b_x[k + iy];
        for (ix = 0; ix < k; ix++) {
          b = ix + jy;
          y[b] -= y[i] * b_x[ix + iy];
        }
      }
    }
  }
}

//
// Arguments    : const double u_i_data[]
//                const double x_i_data[]
//                double z_i_data[]
//                const double u_k_i_data[]
//                const double x_k_i_data[]
//                const double z_k_i_data[]
//                double rho
//                double *stepSizeMaxZ_i
//                double *stepSizeMaxG_i
// Return Type  : void
//
static void c_fraction_to_boundary_parallel(const double u_i_data[], const
  double x_i_data[], double z_i_data[], const double u_k_i_data[], const double
  x_k_i_data[], const double z_k_i_data[], double rho, double *stepSizeMaxZ_i,
  double *stepSizeMaxG_i)
{
  int z_i_tmp;
  double stepSizeZ_i_data[174];
  double b_z_i_tmp;
  int stepSizeZ_i_size[2];
  double z_i[29];
  double tmp_data[174];
  double u_k_i[29];

  //  line search variables
  //  update z
  for (int j = 0; j < 6; j++) {
    int c_z_i_tmp;
    int d_z_i_tmp;
    int e_z_i_tmp;
    int f_z_i_tmp;
    int g_z_i_tmp;
    int h_z_i_tmp;
    int i_z_i_tmp;
    int j_z_i_tmp;
    int k_z_i_tmp;
    int l_z_i_tmp;
    int m_z_i_tmp;
    int n_z_i_tmp;
    int o_z_i_tmp;

    // OCP_GEN_G
    //     G = OCP_GEN_G(IN1,IN2,IN3)
    //     This function was generated by the Symbolic Math Toolbox version 8.5. 
    //     06-Nov-2020 17:34:34
    // OCP_GEN_G
    //     G = OCP_GEN_G(IN1,IN2,IN3)
    //     This function was generated by the Symbolic Math Toolbox version 8.5. 
    //     06-Nov-2020 17:34:34
    b_z_i_tmp = u_i_data[8 * j];
    z_i[0] = z_i_data[29 * j] * (b_z_i_tmp + 10.0) - rho;
    z_i_tmp = 8 * j + 1;
    z_i[1] = z_i_data[29 * j + 1] * (u_i_data[z_i_tmp] + 10.0) - rho;
    c_z_i_tmp = 8 * j + 2;
    z_i[2] = z_i_data[29 * j + 2] * (u_i_data[c_z_i_tmp] + 10.0) - rho;
    d_z_i_tmp = 8 * j + 3;
    z_i[3] = z_i_data[29 * j + 3] * (u_i_data[d_z_i_tmp] + 10.0) - rho;
    e_z_i_tmp = 8 * j + 4;
    z_i[4] = z_i_data[29 * j + 4] * (u_i_data[e_z_i_tmp] + 10.0) - rho;
    f_z_i_tmp = 8 * j + 5;
    z_i[5] = z_i_data[29 * j + 5] * (u_i_data[f_z_i_tmp] + 10.0) - rho;
    g_z_i_tmp = 8 * j + 6;
    z_i[6] = z_i_data[29 * j + 6] * (u_i_data[g_z_i_tmp] + 10.0) - rho;
    z_i[7] = z_i_data[29 * j + 7] * (-b_z_i_tmp + 10.0) - rho;
    z_i[8] = z_i_data[29 * j + 8] * (-u_i_data[z_i_tmp] + 10.0) - rho;
    z_i[9] = z_i_data[29 * j + 9] * (-u_i_data[c_z_i_tmp] + 10.0) - rho;
    z_i[10] = z_i_data[29 * j + 10] * (-u_i_data[d_z_i_tmp] + 10.0) - rho;
    z_i[11] = z_i_data[29 * j + 11] * (-u_i_data[e_z_i_tmp] + 10.0) - rho;
    z_i[12] = z_i_data[29 * j + 12] * (-u_i_data[f_z_i_tmp] + 10.0) - rho;
    z_i[13] = z_i_data[29 * j + 13] * (-u_i_data[g_z_i_tmp] + 10.0) - rho;
    h_z_i_tmp = 8 * j + 7;
    z_i[14] = z_i_data[29 * j + 14] * u_i_data[h_z_i_tmp] - rho;
    i_z_i_tmp = 14 * j + 7;
    z_i[15] = z_i_data[29 * j + 15] * ((u_i_data[h_z_i_tmp] + 1.5707963267948966)
      + x_i_data[i_z_i_tmp]) - rho;
    j_z_i_tmp = 14 * j + 8;
    z_i[16] = z_i_data[29 * j + 16] * ((u_i_data[h_z_i_tmp] + 1.5707963267948966)
      + x_i_data[j_z_i_tmp]) - rho;
    k_z_i_tmp = 14 * j + 9;
    z_i[17] = z_i_data[29 * j + 17] * ((u_i_data[h_z_i_tmp] + 1.5707963267948966)
      + x_i_data[k_z_i_tmp]) - rho;
    l_z_i_tmp = 14 * j + 10;
    z_i[18] = z_i_data[29 * j + 18] * ((u_i_data[h_z_i_tmp] + 1.5707963267948966)
      + x_i_data[l_z_i_tmp]) - rho;
    m_z_i_tmp = 14 * j + 11;
    z_i[19] = z_i_data[29 * j + 19] * ((u_i_data[h_z_i_tmp] + 1.5707963267948966)
      + x_i_data[m_z_i_tmp]) - rho;
    n_z_i_tmp = 14 * j + 12;
    z_i[20] = z_i_data[29 * j + 20] * ((u_i_data[h_z_i_tmp] + 1.5707963267948966)
      + x_i_data[n_z_i_tmp]) - rho;
    o_z_i_tmp = 14 * j + 13;
    z_i[21] = z_i_data[29 * j + 21] * ((u_i_data[h_z_i_tmp] + 1.5707963267948966)
      + x_i_data[o_z_i_tmp]) - rho;
    z_i[22] = z_i_data[29 * j + 22] * ((u_i_data[h_z_i_tmp] + 1.5707963267948966)
      - x_i_data[i_z_i_tmp]) - rho;
    z_i[23] = z_i_data[29 * j + 23] * ((u_i_data[h_z_i_tmp] + 1.5707963267948966)
      - x_i_data[j_z_i_tmp]) - rho;
    z_i[24] = z_i_data[29 * j + 24] * ((u_i_data[h_z_i_tmp] + 1.5707963267948966)
      - x_i_data[k_z_i_tmp]) - rho;
    z_i[25] = z_i_data[29 * j + 25] * ((u_i_data[h_z_i_tmp] + 1.5707963267948966)
      - x_i_data[l_z_i_tmp]) - rho;
    z_i[26] = z_i_data[29 * j + 26] * ((u_i_data[h_z_i_tmp] + 1.5707963267948966)
      - x_i_data[m_z_i_tmp]) - rho;
    z_i[27] = z_i_data[29 * j + 27] * ((u_i_data[h_z_i_tmp] + 1.5707963267948966)
      - x_i_data[n_z_i_tmp]) - rho;
    z_i[28] = z_i_data[29 * j + 28] * ((u_i_data[h_z_i_tmp] + 1.5707963267948966)
      - x_i_data[o_z_i_tmp]) - rho;
    b_z_i_tmp = u_k_i_data[8 * j];
    u_k_i[0] = b_z_i_tmp + 10.0;
    u_k_i[1] = u_k_i_data[z_i_tmp] + 10.0;
    u_k_i[2] = u_k_i_data[c_z_i_tmp] + 10.0;
    u_k_i[3] = u_k_i_data[d_z_i_tmp] + 10.0;
    u_k_i[4] = u_k_i_data[e_z_i_tmp] + 10.0;
    u_k_i[5] = u_k_i_data[f_z_i_tmp] + 10.0;
    u_k_i[6] = u_k_i_data[g_z_i_tmp] + 10.0;
    u_k_i[7] = -b_z_i_tmp + 10.0;
    u_k_i[8] = -u_k_i_data[z_i_tmp] + 10.0;
    u_k_i[9] = -u_k_i_data[c_z_i_tmp] + 10.0;
    u_k_i[10] = -u_k_i_data[d_z_i_tmp] + 10.0;
    u_k_i[11] = -u_k_i_data[e_z_i_tmp] + 10.0;
    u_k_i[12] = -u_k_i_data[f_z_i_tmp] + 10.0;
    u_k_i[13] = -u_k_i_data[g_z_i_tmp] + 10.0;
    u_k_i[14] = u_k_i_data[h_z_i_tmp];
    u_k_i[15] = (u_k_i_data[h_z_i_tmp] + 1.5707963267948966) +
      x_k_i_data[i_z_i_tmp];
    u_k_i[16] = (u_k_i_data[h_z_i_tmp] + 1.5707963267948966) +
      x_k_i_data[j_z_i_tmp];
    u_k_i[17] = (u_k_i_data[h_z_i_tmp] + 1.5707963267948966) +
      x_k_i_data[k_z_i_tmp];
    u_k_i[18] = (u_k_i_data[h_z_i_tmp] + 1.5707963267948966) +
      x_k_i_data[l_z_i_tmp];
    u_k_i[19] = (u_k_i_data[h_z_i_tmp] + 1.5707963267948966) +
      x_k_i_data[m_z_i_tmp];
    u_k_i[20] = (u_k_i_data[h_z_i_tmp] + 1.5707963267948966) +
      x_k_i_data[n_z_i_tmp];
    u_k_i[21] = (u_k_i_data[h_z_i_tmp] + 1.5707963267948966) +
      x_k_i_data[o_z_i_tmp];
    u_k_i[22] = (u_k_i_data[h_z_i_tmp] + 1.5707963267948966) -
      x_k_i_data[i_z_i_tmp];
    u_k_i[23] = (u_k_i_data[h_z_i_tmp] + 1.5707963267948966) -
      x_k_i_data[j_z_i_tmp];
    u_k_i[24] = (u_k_i_data[h_z_i_tmp] + 1.5707963267948966) -
      x_k_i_data[k_z_i_tmp];
    u_k_i[25] = (u_k_i_data[h_z_i_tmp] + 1.5707963267948966) -
      x_k_i_data[l_z_i_tmp];
    u_k_i[26] = (u_k_i_data[h_z_i_tmp] + 1.5707963267948966) -
      x_k_i_data[m_z_i_tmp];
    u_k_i[27] = (u_k_i_data[h_z_i_tmp] + 1.5707963267948966) -
      x_k_i_data[n_z_i_tmp];
    u_k_i[28] = (u_k_i_data[h_z_i_tmp] + 1.5707963267948966) -
      x_k_i_data[o_z_i_tmp];
    for (c_z_i_tmp = 0; c_z_i_tmp < 29; c_z_i_tmp++) {
      z_i_tmp = c_z_i_tmp + 29 * j;
      b_z_i_tmp = z_i_data[z_i_tmp] - z_i[c_z_i_tmp] / u_k_i[c_z_i_tmp];
      z_i[c_z_i_tmp] = b_z_i_tmp;
      z_i_data[z_i_tmp] = b_z_i_tmp;
    }
  }

  //         %% Line Search for feasibility
  //  z
  for (z_i_tmp = 0; z_i_tmp < 174; z_i_tmp++) {
    b_z_i_tmp = -0.95 * (z_k_i_data[z_i_tmp] / (z_i_data[z_i_tmp] -
      z_k_i_data[z_i_tmp]));
    stepSizeZ_i_data[z_i_tmp] = b_z_i_tmp;
    if ((b_z_i_tmp > 1.0) || (b_z_i_tmp < 0.0)) {
      stepSizeZ_i_data[z_i_tmp] = 1.0;
    }
  }

  *stepSizeMaxZ_i = stepSizeZ_i_data[0];
  for (z_i_tmp = 0; z_i_tmp < 173; z_i_tmp++) {
    b_z_i_tmp = stepSizeZ_i_data[z_i_tmp + 1];
    if (*stepSizeMaxZ_i > b_z_i_tmp) {
      *stepSizeMaxZ_i = b_z_i_tmp;
    }
  }

  //  G
  OCP_G(u_k_i_data, x_k_i_data, stepSizeZ_i_data, stepSizeZ_i_size);

  //  GMin
  OCP_G(u_i_data, x_i_data, tmp_data, stepSizeZ_i_size);
  for (z_i_tmp = 0; z_i_tmp < 174; z_i_tmp++) {
    b_z_i_tmp = -0.95 * (stepSizeZ_i_data[z_i_tmp] / (tmp_data[z_i_tmp] -
      stepSizeZ_i_data[z_i_tmp]));
    stepSizeZ_i_data[z_i_tmp] = b_z_i_tmp;
    if ((b_z_i_tmp > 1.0) || (b_z_i_tmp < 0.0)) {
      stepSizeZ_i_data[z_i_tmp] = 1.0;
    }
  }

  *stepSizeMaxG_i = stepSizeZ_i_data[0];
  for (z_i_tmp = 0; z_i_tmp < 173; z_i_tmp++) {
    b_z_i_tmp = stepSizeZ_i_data[z_i_tmp + 1];
    if (*stepSizeMaxG_i > b_z_i_tmp) {
      *stepSizeMaxG_i = b_z_i_tmp;
    }
  }
}

//
// Arguments    : double fid
// Return Type  : int
//
static int cfclose(double fid)
{
  int st;
  signed char fileid;
  signed char b_fileid;
  FILE * filestar;
  st = -1;
  fileid = static_cast<signed char>(rt_roundd(fid));
  if ((fileid < 0) || (fid != fileid)) {
    fileid = -1;
  }

  b_fileid = fileid;
  if (fileid < 0) {
    b_fileid = -1;
  }

  if (b_fileid >= 3) {
    filestar = eml_openfiles[b_fileid - 3];
  } else if (b_fileid == 0) {
    filestar = stdin;
  } else if (b_fileid == 1) {
    filestar = stdout;
  } else if (b_fileid == 2) {
    filestar = stderr;
  } else {
    filestar = NULL;
  }

  if ((filestar != NULL) && (fileid >= 3)) {
    int cst;
    cst = fclose(filestar);
    if (cst == 0) {
      st = 0;
      cst = fileid - 3;
      eml_openfiles[cst] = NULL;
      eml_autoflush[cst] = true;
    }
  }

  return st;
}

//
// Arguments    : const char * cfilename
//                const char * cpermission
// Return Type  : signed char
//
static signed char cfopen(const char * cfilename, const char * cpermission)
{
  signed char fileid;
  signed char j;
  fileid = -1;
  j = filedata();
  if (j >= 1) {
    FILE * filestar;
    filestar = fopen(cfilename, cpermission);
    if (filestar != NULL) {
      int i;
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      i = j + 2;
      if (i > 127) {
        i = 127;
      }

      fileid = static_cast<signed char>(i);
    }
  }

  return fileid;
}

//
// Arguments    : double lambda_i_data[]
//                double u_i_data[]
//                double x_i_data[]
//                const double z_i_data[]
//                const double p_i_data[]
//                double xPrev_i_data[]
//                double lambdaNext_i_data[]
//                double LAMBDA_i_data[]
//                double rho
//                double i
//                double p_muu_F_i_data[]
//                int p_muu_F_i_size[3]
//                double p_muu_Lambda_i_data[]
//                int p_muu_Lambda_i_size[3]
//                double p_lambda_Lambda_i_data[]
//                int p_lambda_Lambda_i_size[3]
//                double p_x_Lambda_i_data[]
//                int p_x_Lambda_i_size[3]
//                double p_x_F_i_data[]
//                int p_x_F_i_size[3]
//                double KKTxEquation_i_data[]
//                int KKTxEquation_i_size[2]
//                double KKTC_i_data[]
//                int KKTC_i_size[2]
//                double KKTHu_i_data[]
//                int KKTHu_i_size[2]
//                double KKTlambdaEquation_i_data[]
//                int KKTlambdaEquation_i_size[2]
//                double L_i_data[]
//                int L_i_size[2]
//                double LB_i_data[]
//                int LB_i_size[2]
// Return Type  : void
//
static void coarse_update_func(double lambda_i_data[], double u_i_data[], double
  x_i_data[], const double z_i_data[], const double p_i_data[], double
  xPrev_i_data[], double lambdaNext_i_data[], double LAMBDA_i_data[], double rho,
  double i, double p_muu_F_i_data[], int p_muu_F_i_size[3], double
  p_muu_Lambda_i_data[], int p_muu_Lambda_i_size[3], double
  p_lambda_Lambda_i_data[], int p_lambda_Lambda_i_size[3], double
  p_x_Lambda_i_data[], int p_x_Lambda_i_size[3], double p_x_F_i_data[], int
  p_x_F_i_size[3], double KKTxEquation_i_data[], int KKTxEquation_i_size[2],
  double KKTC_i_data[], int KKTC_i_size[2], double KKTHu_i_data[], int
  KKTHu_i_size[2], double KKTlambdaEquation_i_data[], int
  KKTlambdaEquation_i_size[2], double L_i_data[], int L_i_size[2], double
  LB_i_data[], int LB_i_size[2])
{
  double dv[8];
  double dv1[8];
  double dv2[8];
  double dv3[8];
  double dv4[8];
  double t17;
  double t18;
  double t19;
  double t21;
  double xEq_j_i[14];
  double p_x_muu_j_i[112];
  double Fx_j_i[196];
  double lambdaEq_j_i[14];
  double LAMBDAUncrt_j_i[14];
  double HAlluT_j_i[8];
  double p_muu_F_j_i[112];
  double dmu_u_j_i[8];
  double invFx_j_i[196];
  double b_LAMBDAUncrt_j_i[196];
  double FuT_invFxT_j_i[112];
  double b_Fx_j_i[196];
  double b_p_muu_F_j_i[8];
  double Aux_invFx_j_i[112];
  double FuT_LAMBDAUncrt_m_Aux_invFx_j_i[112];
  double MT_j_i[64];
  double b_HAlluT_j_i[8];
  double z_i[64];
  double reshapes_f2[64];
  double b_reshapes_f2[112];

  //  regularization
  p_muu_F_i_size[0] = 8;
  p_muu_F_i_size[1] = 14;
  p_muu_F_i_size[2] = 6;
  p_muu_Lambda_i_size[0] = 8;
  p_muu_Lambda_i_size[1] = 14;
  p_muu_Lambda_i_size[2] = 6;
  p_lambda_Lambda_i_size[0] = 14;
  p_lambda_Lambda_i_size[1] = 14;
  p_lambda_Lambda_i_size[2] = 6;
  p_x_Lambda_i_size[0] = 14;
  p_x_Lambda_i_size[1] = 14;
  p_x_Lambda_i_size[2] = 6;
  p_x_F_i_size[0] = 14;
  p_x_F_i_size[1] = 14;
  p_x_F_i_size[2] = 6;
  L_i_size[0] = 1;
  L_i_size[1] = 6;
  LB_i_size[0] = 1;
  LB_i_size[1] = 6;
  KKTxEquation_i_size[0] = 1;
  KKTxEquation_i_size[1] = 6;
  KKTC_i_size[0] = 1;
  KKTC_i_size[1] = 6;
  KKTHu_i_size[0] = 1;
  KKTHu_i_size[1] = 6;
  KKTlambdaEquation_i_size[0] = 1;
  KKTlambdaEquation_i_size[1] = 6;
  dv[0] = 0.0;
  dv[1] = 0.0;
  dv[2] = 0.0;
  dv[4] = 0.0;
  dv[5] = 0.0;
  dv[6] = 0.0;
  dv[7] = 0.0;
  dv1[0] = 0.0;
  dv1[1] = 0.0;
  dv1[2] = 0.0;
  dv1[3] = 0.0;
  dv1[5] = 0.0;
  dv1[6] = 0.0;
  dv1[7] = 0.0;
  dv2[0] = 0.0;
  dv2[1] = 0.0;
  dv2[2] = 0.0;
  dv2[3] = 0.0;
  dv2[4] = 0.0;
  dv2[6] = 0.0;
  dv2[7] = 0.0;
  dv3[0] = 0.0;
  dv3[1] = 0.0;
  dv3[2] = 0.0;
  dv3[3] = 0.0;
  dv3[4] = 0.0;
  dv3[5] = 0.0;
  dv3[7] = 0.0;
  dv4[0] = 0.0;
  dv4[1] = 0.0;
  dv4[2] = 0.0;
  dv4[3] = 0.0;
  dv4[4] = 0.0;
  dv4[5] = 0.0;
  dv4[6] = 0.0;
  for (int j = 0; j < 6; j++) {
    int coffset;
    int t2_tmp_tmp;
    double t2;
    double t20;
    double t4;
    double t22;
    double t6;
    double t23;
    double t24;
    double t8;
    double t25;
    double t26;
    double t10;
    double t27;
    double t28;
    double t12;
    double t29;
    double t30;
    double t14;
    int L_i_tmp_tmp;
    double L_i_tmp;
    int b_L_i_tmp_tmp;
    double b_L_i_tmp;
    int c_L_i_tmp_tmp;
    double c_L_i_tmp;
    int d_L_i_tmp_tmp;
    double d_L_i_tmp;
    int e_L_i_tmp_tmp;
    double e_L_i_tmp;
    int f_L_i_tmp_tmp;
    double f_L_i_tmp;
    int g_L_i_tmp_tmp;
    double g_L_i_tmp;
    int aoffset;
    double h_L_i_tmp;
    int h_L_i_tmp_tmp;
    double i_L_i_tmp;
    int i_L_i_tmp_tmp;
    double j_L_i_tmp;
    int j_L_i_tmp_tmp;
    double k_L_i_tmp;
    int k_L_i_tmp_tmp;
    double l_L_i_tmp;
    int l_L_i_tmp_tmp;
    double m_L_i_tmp;
    int m_L_i_tmp_tmp;
    double n_L_i_tmp;
    int n_L_i_tmp_tmp;
    double o_L_i_tmp;
    double t31;
    double t32;
    double t33;
    double t34;
    double t35;
    double t36;
    double t37;
    double t39;
    double t41;
    double t43;
    double t45;
    double t47;
    double t49;
    double t51;
    int b_i;
    int i1;
    int t31_tmp;

    //  Function and Jacobian
    // OCP_GEN_L_LU_LX
    //     [L,LU,LX] = OCP_GEN_L_LU_LX(IN1,IN2,IN3)
    //     This function was generated by the Symbolic Math Toolbox version 8.5. 
    //     06-Nov-2020 17:34:34
    coffset = 7 * (5 - j);
    t17 = p_i_data[coffset];
    t2_tmp_tmp = 14 * (5 - j);
    t18 = x_i_data[t2_tmp_tmp];
    t2 = t17 - t18;
    t19 = p_i_data[coffset + 1];
    t20 = x_i_data[t2_tmp_tmp + 1];
    t4 = t19 - t20;
    t21 = p_i_data[coffset + 2];
    t22 = x_i_data[t2_tmp_tmp + 2];
    t6 = t21 - t22;
    t23 = p_i_data[coffset + 3];
    t24 = x_i_data[t2_tmp_tmp + 3];
    t8 = t23 - t24;
    t25 = p_i_data[coffset + 4];
    t26 = x_i_data[t2_tmp_tmp + 4];
    t10 = t25 - t26;
    t27 = p_i_data[coffset + 5];
    t28 = x_i_data[t2_tmp_tmp + 5];
    t12 = t27 - t28;
    t29 = p_i_data[coffset + 6];
    t30 = x_i_data[t2_tmp_tmp + 6];
    t14 = t29 - t30;
    L_i_tmp_tmp = 8 * (5 - j);
    L_i_tmp = u_i_data[L_i_tmp_tmp];
    b_L_i_tmp_tmp = L_i_tmp_tmp + 1;
    b_L_i_tmp = u_i_data[b_L_i_tmp_tmp];
    c_L_i_tmp_tmp = L_i_tmp_tmp + 2;
    c_L_i_tmp = u_i_data[c_L_i_tmp_tmp];
    d_L_i_tmp_tmp = L_i_tmp_tmp + 3;
    d_L_i_tmp = u_i_data[d_L_i_tmp_tmp];
    e_L_i_tmp_tmp = L_i_tmp_tmp + 4;
    e_L_i_tmp = u_i_data[e_L_i_tmp_tmp];
    f_L_i_tmp_tmp = L_i_tmp_tmp + 5;
    f_L_i_tmp = u_i_data[f_L_i_tmp_tmp];
    g_L_i_tmp_tmp = L_i_tmp_tmp + 6;
    g_L_i_tmp = u_i_data[g_L_i_tmp_tmp];
    aoffset = L_i_tmp_tmp + 7;
    h_L_i_tmp = u_i_data[aoffset];
    h_L_i_tmp_tmp = t2_tmp_tmp + 7;
    i_L_i_tmp = x_i_data[h_L_i_tmp_tmp];
    i_L_i_tmp_tmp = t2_tmp_tmp + 8;
    j_L_i_tmp = x_i_data[i_L_i_tmp_tmp];
    j_L_i_tmp_tmp = t2_tmp_tmp + 9;
    k_L_i_tmp = x_i_data[j_L_i_tmp_tmp];
    k_L_i_tmp_tmp = t2_tmp_tmp + 10;
    l_L_i_tmp = x_i_data[k_L_i_tmp_tmp];
    l_L_i_tmp_tmp = t2_tmp_tmp + 11;
    m_L_i_tmp = x_i_data[l_L_i_tmp_tmp];
    m_L_i_tmp_tmp = t2_tmp_tmp + 12;
    n_L_i_tmp = x_i_data[m_L_i_tmp_tmp];
    n_L_i_tmp_tmp = t2_tmp_tmp + 13;
    o_L_i_tmp = x_i_data[n_L_i_tmp_tmp];
    L_i_data[5 - j] = ((((((((((((((((((((L_i_tmp * L_i_tmp / 2000.0 + b_L_i_tmp
      * b_L_i_tmp / 2000.0) + c_L_i_tmp * c_L_i_tmp / 2000.0) + d_L_i_tmp *
      d_L_i_tmp / 2000.0) + e_L_i_tmp * e_L_i_tmp / 2000.0) + f_L_i_tmp *
      f_L_i_tmp / 2000.0) + g_L_i_tmp * g_L_i_tmp / 2000.0) + h_L_i_tmp *
      h_L_i_tmp * 1000.0) + i_L_i_tmp * i_L_i_tmp / 20.0) + j_L_i_tmp *
      j_L_i_tmp / 20.0) + k_L_i_tmp * k_L_i_tmp / 20.0) + l_L_i_tmp * l_L_i_tmp /
      20.0) + m_L_i_tmp * m_L_i_tmp / 20.0) + n_L_i_tmp * n_L_i_tmp / 20.0) +
      o_L_i_tmp * o_L_i_tmp / 20.0) + t2 * (t17 / 2.0 - t18 / 2.0)) + t4 * (t19 /
      2.0 - t20 / 2.0)) + t6 * (t21 / 2.0 - t22 / 2.0)) + t8 * (t23 / 2.0 - t24 /
      2.0)) + t10 * (t25 / 2.0 - t26 / 2.0)) + t12 * (t27 / 2.0 - t28 / 2.0)) +
      t14 * (t29 / 2.0 - t30 / 2.0);

    // OCP_GEN_LB_LBU_LBX
    //     [LB,LBU,LBX] = OCP_GEN_LB_LBU_LBX(IN1,IN2,IN3)
    //     This function was generated by the Symbolic Math Toolbox version 8.5. 
    //     06-Nov-2020 17:34:33
    t17 = (h_L_i_tmp + 1.5707963267948966) + i_L_i_tmp;
    t18 = (u_i_data[aoffset] + 1.5707963267948966) + j_L_i_tmp;
    t19 = (u_i_data[aoffset] + 1.5707963267948966) + k_L_i_tmp;
    t20 = (u_i_data[aoffset] + 1.5707963267948966) + l_L_i_tmp;
    t21 = (u_i_data[aoffset] + 1.5707963267948966) + m_L_i_tmp;
    t22 = (u_i_data[aoffset] + 1.5707963267948966) + n_L_i_tmp;
    t23 = (u_i_data[aoffset] + 1.5707963267948966) + o_L_i_tmp;
    t24 = (-i_L_i_tmp + 1.5707963267948966) + h_L_i_tmp;
    t25 = (-j_L_i_tmp + 1.5707963267948966) + h_L_i_tmp;
    t26 = (-k_L_i_tmp + 1.5707963267948966) + h_L_i_tmp;
    t27 = (-l_L_i_tmp + 1.5707963267948966) + h_L_i_tmp;
    t28 = (-m_L_i_tmp + 1.5707963267948966) + h_L_i_tmp;
    t29 = (-n_L_i_tmp + 1.5707963267948966) + h_L_i_tmp;
    t30 = (-o_L_i_tmp + 1.5707963267948966) + h_L_i_tmp;
    t31 = 1.0 / t17;
    t32 = 1.0 / t18;
    t33 = 1.0 / t19;
    t34 = 1.0 / t20;
    t35 = 1.0 / t21;
    t36 = 1.0 / t22;
    t37 = 1.0 / t23;
    t39 = 1.0 / t24;
    t41 = 1.0 / t25;
    t43 = 1.0 / t26;
    t45 = 1.0 / t27;
    t47 = 1.0 / t28;
    t49 = 1.0 / t29;
    t51 = 1.0 / t30;
    LB_i_data[5 - j] = ((((((((((((((((((((((((((((((h_L_i_tmp * 0.0015 +
      0.0021991148575128553) - std::log(-L_i_tmp + 10.0)) - std::log(-b_L_i_tmp
      + 10.0)) - std::log(-c_L_i_tmp + 10.0)) - std::log(-d_L_i_tmp + 10.0)) -
      std::log(-e_L_i_tmp + 10.0)) - std::log(-f_L_i_tmp + 10.0)) - std::log
      (-g_L_i_tmp + 10.0)) - std::log(L_i_tmp + 10.0)) - std::log(b_L_i_tmp +
      10.0)) - std::log(c_L_i_tmp + 10.0)) - std::log(d_L_i_tmp + 10.0)) - std::
      log(e_L_i_tmp + 10.0)) - std::log(f_L_i_tmp + 10.0)) - std::log(g_L_i_tmp
      + 10.0)) - std::log(t17)) - std::log(t18)) - std::log(t19)) - std::log(t20))
      - std::log(t21)) - std::log(t22)) - std::log(t23)) - std::log(t24)) - std::
      log(t25)) - std::log(t26)) - std::log(t27)) - std::log(t28)) - std::log
                          (t29)) - std::log(t30)) - std::log(h_L_i_tmp)) + 0.014;
    OCP_F_Fu_Fx(*(double (*)[8])&u_i_data[8 * (5 - j)], *(double (*)[14])&
                x_i_data[14 * (5 - j)], i, xEq_j_i, p_x_muu_j_i, Fx_j_i);

    //  KKT
    if (6 - j > 1) {
      std::memcpy(&xPrev_i_data[t2_tmp_tmp], &x_i_data[j * -14 + 56], 14U *
                  sizeof(double));
    }

    if (6 - j < 6) {
      std::memcpy(&lambdaNext_i_data[t2_tmp_tmp], &lambda_i_data[j * -14 + 84],
                  14U * sizeof(double));
    }

    lambdaEq_j_i[0] = -t2;
    lambdaEq_j_i[1] = -t4;
    lambdaEq_j_i[2] = -t6;
    lambdaEq_j_i[3] = -t8;
    lambdaEq_j_i[4] = -t10;
    lambdaEq_j_i[5] = -t12;
    lambdaEq_j_i[6] = -t14;
    lambdaEq_j_i[7] = i_L_i_tmp / 10.0;
    lambdaEq_j_i[8] = j_L_i_tmp / 10.0;
    lambdaEq_j_i[9] = k_L_i_tmp / 10.0;
    lambdaEq_j_i[10] = l_L_i_tmp / 10.0;
    lambdaEq_j_i[11] = m_L_i_tmp / 10.0;
    lambdaEq_j_i[12] = n_L_i_tmp / 10.0;
    lambdaEq_j_i[13] = o_L_i_tmp / 10.0;
    LAMBDAUncrt_j_i[0] = 0.0;
    LAMBDAUncrt_j_i[1] = 0.0;
    LAMBDAUncrt_j_i[2] = 0.0;
    LAMBDAUncrt_j_i[3] = 0.0;
    LAMBDAUncrt_j_i[4] = 0.0;
    LAMBDAUncrt_j_i[5] = 0.0;
    LAMBDAUncrt_j_i[6] = 0.0;
    LAMBDAUncrt_j_i[7] = rho * (-t31 + t39);
    LAMBDAUncrt_j_i[8] = rho * (-t32 + t41);
    LAMBDAUncrt_j_i[9] = rho * (-t33 + t43);
    LAMBDAUncrt_j_i[10] = rho * (-t34 + t45);
    LAMBDAUncrt_j_i[11] = rho * (-t35 + t47);
    LAMBDAUncrt_j_i[12] = rho * (-t36 + t49);
    LAMBDAUncrt_j_i[13] = rho * (-t37 + t51);
    for (b_i = 0; b_i < 14; b_i++) {
      coffset = b_i + t2_tmp_tmp;
      xEq_j_i[b_i] += xPrev_i_data[coffset];
      for (i1 = 0; i1 < 8; i1++) {
        p_muu_F_j_i[i1 + (b_i << 3)] = p_x_muu_j_i[b_i + 14 * i1];
      }

      t21 = 0.0;
      for (i1 = 0; i1 < 14; i1++) {
        t21 += Fx_j_i[i1 + 14 * b_i] * lambda_i_data[i1 + t2_tmp_tmp];
      }

      lambdaEq_j_i[b_i] = (lambdaNext_i_data[coffset] + (lambdaEq_j_i[b_i] + t21))
        + LAMBDAUncrt_j_i[b_i];
    }

    HAlluT_j_i[0] = L_i_tmp / 1000.0;
    HAlluT_j_i[1] = b_L_i_tmp / 1000.0;
    HAlluT_j_i[2] = c_L_i_tmp / 1000.0;
    HAlluT_j_i[3] = d_L_i_tmp / 1000.0;
    HAlluT_j_i[4] = e_L_i_tmp / 1000.0;
    HAlluT_j_i[5] = f_L_i_tmp / 1000.0;
    HAlluT_j_i[6] = g_L_i_tmp / 1000.0;
    HAlluT_j_i[7] = h_L_i_tmp * 2000.0;
    dmu_u_j_i[0] = rho * (-1.0 / (L_i_tmp - 10.0) - 1.0 / (u_i_data[L_i_tmp_tmp]
      + 10.0));
    dmu_u_j_i[1] = rho * (-1.0 / (b_L_i_tmp - 10.0) - 1.0 /
                          (u_i_data[b_L_i_tmp_tmp] + 10.0));
    dmu_u_j_i[2] = rho * (-1.0 / (c_L_i_tmp - 10.0) - 1.0 /
                          (u_i_data[c_L_i_tmp_tmp] + 10.0));
    dmu_u_j_i[3] = rho * (-1.0 / (d_L_i_tmp - 10.0) - 1.0 /
                          (u_i_data[d_L_i_tmp_tmp] + 10.0));
    dmu_u_j_i[4] = rho * (-1.0 / (e_L_i_tmp - 10.0) - 1.0 /
                          (u_i_data[e_L_i_tmp_tmp] + 10.0));
    dmu_u_j_i[5] = rho * (-1.0 / (f_L_i_tmp - 10.0) - 1.0 /
                          (u_i_data[f_L_i_tmp_tmp] + 10.0));
    dmu_u_j_i[6] = rho * (-1.0 / (g_L_i_tmp - 10.0) - 1.0 /
                          (u_i_data[g_L_i_tmp_tmp] + 10.0));
    dmu_u_j_i[7] = rho * (((((((((((((((-t31 - t39) + -t32) - t41) + -t33) - t43)
      + -t34) - t45) + -t35) - t47) + -t36) - t49) + -t37) - t51) - 1.0 /
      h_L_i_tmp) + 0.0015);
    for (b_i = 0; b_i < 8; b_i++) {
      t21 = 0.0;
      for (i1 = 0; i1 < 14; i1++) {
        t21 += p_muu_F_j_i[b_i + (i1 << 3)] * lambda_i_data[i1 + t2_tmp_tmp];
      }

      HAlluT_j_i[b_i] = (HAlluT_j_i[b_i] + t21) + dmu_u_j_i[b_i];
    }

    //  Hessian
    //  --- AuuCondensed_j_i = Auu_j_i + Gu_j_i.'*(z_j_i./G_j_i.*Gu_j_i);
    //  --- AuxCondensed_j_i = Aux_j_i + Gu_j_i.'*(z_j_i./G_j_i.*Gx_j_i);
    //  --- AxxCondensed_j_i = Axx_j_i + Gx_j_i.'*(z_j_i./G_j_i.*Gx_j_i);
    // OCP_GEN_AUU_AUX_AXX_CONDENSED
    //     [AUUCONDENSED,AUXCONDENSED,AXXCONDENSED] = OCP_GEN_AUU_AUX_AXX_CONDENSED(IN1,IN2,IN3,IN4,IN5,IN6) 
    //     This function was generated by the Symbolic Math Toolbox version 8.5. 
    //     06-Nov-2020 17:34:36
    t31_tmp = 29 * (5 - j);
    t31 = 1.0 / ((u_i_data[aoffset] + 1.5707963267948966) +
                 x_i_data[h_L_i_tmp_tmp]) * z_i_data[t31_tmp + 15];
    t32 = 1.0 / ((u_i_data[aoffset] + 1.5707963267948966) +
                 x_i_data[i_L_i_tmp_tmp]) * z_i_data[t31_tmp + 16];
    t33 = 1.0 / ((u_i_data[aoffset] + 1.5707963267948966) +
                 x_i_data[j_L_i_tmp_tmp]) * z_i_data[t31_tmp + 17];
    t34 = 1.0 / ((u_i_data[aoffset] + 1.5707963267948966) +
                 x_i_data[k_L_i_tmp_tmp]) * z_i_data[t31_tmp + 18];
    t35 = 1.0 / ((u_i_data[aoffset] + 1.5707963267948966) +
                 x_i_data[l_L_i_tmp_tmp]) * z_i_data[t31_tmp + 19];
    t36 = 1.0 / ((u_i_data[aoffset] + 1.5707963267948966) +
                 x_i_data[m_L_i_tmp_tmp]) * z_i_data[t31_tmp + 20];
    t37 = 1.0 / ((u_i_data[aoffset] + 1.5707963267948966) +
                 x_i_data[n_L_i_tmp_tmp]) * z_i_data[t31_tmp + 21];
    t45 = 1.0 / ((-x_i_data[h_L_i_tmp_tmp] + 1.5707963267948966) +
                 u_i_data[aoffset]) * z_i_data[t31_tmp + 22];
    t18 = 1.0 / ((-x_i_data[i_L_i_tmp_tmp] + 1.5707963267948966) +
                 u_i_data[aoffset]) * z_i_data[t31_tmp + 23];
    t47 = 1.0 / ((-x_i_data[j_L_i_tmp_tmp] + 1.5707963267948966) +
                 u_i_data[aoffset]) * z_i_data[t31_tmp + 24];
    t19 = 1.0 / ((-x_i_data[k_L_i_tmp_tmp] + 1.5707963267948966) +
                 u_i_data[aoffset]) * z_i_data[t31_tmp + 25];
    t49 = 1.0 / ((-x_i_data[l_L_i_tmp_tmp] + 1.5707963267948966) +
                 u_i_data[aoffset]) * z_i_data[t31_tmp + 26];
    t20 = 1.0 / ((-x_i_data[m_L_i_tmp_tmp] + 1.5707963267948966) +
                 u_i_data[aoffset]) * z_i_data[t31_tmp + 27];
    t51 = 1.0 / ((-x_i_data[n_L_i_tmp_tmp] + 1.5707963267948966) +
                 u_i_data[aoffset]) * z_i_data[t31_tmp + 28];

    //  descent regularization
    //  descent regularization
    //  nonsingular regularization
    //  Intermediate Variables
    if (6 - j < 6) {
      for (b_i = 0; b_i < 14; b_i++) {
        for (i1 = 0; i1 < 14; i1++) {
          h_L_i_tmp_tmp = i1 + 14 * b_i;
          LAMBDA_i_data[h_L_i_tmp_tmp + 196 * (5 - j)] =
            LAMBDA_i_data[h_L_i_tmp_tmp + 196 * (6 - j)];
        }
      }
    }

    inv(Fx_j_i, invFx_j_i);
    b_LAMBDAUncrt_j_i[0] = 1.0;
    std::memset(&b_LAMBDAUncrt_j_i[1], 0, 14U * sizeof(double));
    b_LAMBDAUncrt_j_i[15] = 1.0;
    std::memset(&b_LAMBDAUncrt_j_i[16], 0, 14U * sizeof(double));
    b_LAMBDAUncrt_j_i[30] = 1.0;
    std::memset(&b_LAMBDAUncrt_j_i[31], 0, 14U * sizeof(double));
    b_LAMBDAUncrt_j_i[45] = 1.0;
    std::memset(&b_LAMBDAUncrt_j_i[46], 0, 14U * sizeof(double));
    b_LAMBDAUncrt_j_i[60] = 1.0;
    std::memset(&b_LAMBDAUncrt_j_i[61], 0, 14U * sizeof(double));
    b_LAMBDAUncrt_j_i[75] = 1.0;
    std::memset(&b_LAMBDAUncrt_j_i[76], 0, 14U * sizeof(double));
    b_LAMBDAUncrt_j_i[90] = 1.0;
    std::memset(&b_LAMBDAUncrt_j_i[91], 0, 14U * sizeof(double));
    b_LAMBDAUncrt_j_i[105] = (t31 + t45) + 0.1;
    std::memset(&b_LAMBDAUncrt_j_i[106], 0, 14U * sizeof(double));
    b_LAMBDAUncrt_j_i[120] = (t32 + t18) + 0.1;
    std::memset(&b_LAMBDAUncrt_j_i[121], 0, 14U * sizeof(double));
    b_LAMBDAUncrt_j_i[135] = (t33 + t47) + 0.1;
    std::memset(&b_LAMBDAUncrt_j_i[136], 0, 14U * sizeof(double));
    b_LAMBDAUncrt_j_i[150] = (t34 + t19) + 0.1;
    std::memset(&b_LAMBDAUncrt_j_i[151], 0, 14U * sizeof(double));
    b_LAMBDAUncrt_j_i[165] = (t35 + t49) + 0.1;
    std::memset(&b_LAMBDAUncrt_j_i[166], 0, 14U * sizeof(double));
    b_LAMBDAUncrt_j_i[180] = (t36 + t20) + 0.1;
    std::memset(&b_LAMBDAUncrt_j_i[181], 0, 14U * sizeof(double));
    b_LAMBDAUncrt_j_i[195] = (t37 + t51) + 0.1;
    for (b_i = 0; b_i < 14; b_i++) {
      for (i1 = 0; i1 < 14; i1++) {
        coffset = i1 + 14 * b_i;
        Fx_j_i[coffset] = invFx_j_i[b_i + 14 * i1];
        b_LAMBDAUncrt_j_i[coffset] -= LAMBDA_i_data[coffset + 196 * (5 - j)];
      }
    }

    for (b_i = 0; b_i < 14; b_i++) {
      for (i1 = 0; i1 < 14; i1++) {
        t21 = 0.0;
        for (h_L_i_tmp_tmp = 0; h_L_i_tmp_tmp < 14; h_L_i_tmp_tmp++) {
          t21 += Fx_j_i[b_i + 14 * h_L_i_tmp_tmp] *
            b_LAMBDAUncrt_j_i[h_L_i_tmp_tmp + 14 * i1];
        }

        b_Fx_j_i[b_i + 14 * i1] = t21;
      }
    }

    for (b_i = 0; b_i < 14; b_i++) {
      for (i1 = 0; i1 < 14; i1++) {
        t21 = 0.0;
        for (h_L_i_tmp_tmp = 0; h_L_i_tmp_tmp < 14; h_L_i_tmp_tmp++) {
          t21 += b_Fx_j_i[b_i + 14 * h_L_i_tmp_tmp] * invFx_j_i[h_L_i_tmp_tmp +
            14 * i1];
        }

        b_LAMBDAUncrt_j_i[b_i + 14 * i1] = t21;
      }
    }

    std::memset(&FuT_invFxT_j_i[0], 0, 63U * sizeof(double));
    FuT_invFxT_j_i[63] = t31 - t45;
    FuT_invFxT_j_i[64] = 0.0;
    FuT_invFxT_j_i[65] = 0.0;
    FuT_invFxT_j_i[66] = 0.0;
    FuT_invFxT_j_i[67] = 0.0;
    FuT_invFxT_j_i[68] = 0.0;
    FuT_invFxT_j_i[69] = 0.0;
    FuT_invFxT_j_i[70] = 0.0;
    FuT_invFxT_j_i[71] = t32 - t18;
    FuT_invFxT_j_i[72] = 0.0;
    FuT_invFxT_j_i[73] = 0.0;
    FuT_invFxT_j_i[74] = 0.0;
    FuT_invFxT_j_i[75] = 0.0;
    FuT_invFxT_j_i[76] = 0.0;
    FuT_invFxT_j_i[77] = 0.0;
    FuT_invFxT_j_i[78] = 0.0;
    FuT_invFxT_j_i[79] = t33 - t47;
    FuT_invFxT_j_i[80] = 0.0;
    FuT_invFxT_j_i[81] = 0.0;
    FuT_invFxT_j_i[82] = 0.0;
    FuT_invFxT_j_i[83] = 0.0;
    FuT_invFxT_j_i[84] = 0.0;
    FuT_invFxT_j_i[85] = 0.0;
    FuT_invFxT_j_i[86] = 0.0;
    FuT_invFxT_j_i[87] = t34 - t19;
    FuT_invFxT_j_i[88] = 0.0;
    FuT_invFxT_j_i[89] = 0.0;
    FuT_invFxT_j_i[90] = 0.0;
    FuT_invFxT_j_i[91] = 0.0;
    FuT_invFxT_j_i[92] = 0.0;
    FuT_invFxT_j_i[93] = 0.0;
    FuT_invFxT_j_i[94] = 0.0;
    FuT_invFxT_j_i[95] = t35 - t49;
    FuT_invFxT_j_i[96] = 0.0;
    FuT_invFxT_j_i[97] = 0.0;
    FuT_invFxT_j_i[98] = 0.0;
    FuT_invFxT_j_i[99] = 0.0;
    FuT_invFxT_j_i[100] = 0.0;
    FuT_invFxT_j_i[101] = 0.0;
    FuT_invFxT_j_i[102] = 0.0;
    FuT_invFxT_j_i[103] = t36 - t20;
    FuT_invFxT_j_i[104] = 0.0;
    FuT_invFxT_j_i[105] = 0.0;
    FuT_invFxT_j_i[106] = 0.0;
    FuT_invFxT_j_i[107] = 0.0;
    FuT_invFxT_j_i[108] = 0.0;
    FuT_invFxT_j_i[109] = 0.0;
    FuT_invFxT_j_i[110] = 0.0;
    FuT_invFxT_j_i[111] = t37 - t51;
    for (b_i = 0; b_i < 8; b_i++) {
      for (i1 = 0; i1 < 14; i1++) {
        t21 = 0.0;
        for (h_L_i_tmp_tmp = 0; h_L_i_tmp_tmp < 14; h_L_i_tmp_tmp++) {
          t21 += FuT_invFxT_j_i[b_i + (h_L_i_tmp_tmp << 3)] *
            invFx_j_i[h_L_i_tmp_tmp + 14 * i1];
        }

        Aux_invFx_j_i[b_i + (i1 << 3)] = t21;
      }
    }

    for (j_L_i_tmp_tmp = 0; j_L_i_tmp_tmp < 14; j_L_i_tmp_tmp++) {
      coffset = j_L_i_tmp_tmp << 3;
      for (h_L_i_tmp_tmp = 0; h_L_i_tmp_tmp < 8; h_L_i_tmp_tmp++) {
        aoffset = h_L_i_tmp_tmp * 14;
        t17 = 0.0;
        for (i_L_i_tmp_tmp = 0; i_L_i_tmp_tmp < 14; i_L_i_tmp_tmp++) {
          t17 += p_x_muu_j_i[aoffset + i_L_i_tmp_tmp] * invFx_j_i[i_L_i_tmp_tmp *
            14 + j_L_i_tmp_tmp];
        }

        FuT_invFxT_j_i[coffset + h_L_i_tmp_tmp] = t17;
      }
    }

    for (j_L_i_tmp_tmp = 0; j_L_i_tmp_tmp < 8; j_L_i_tmp_tmp++) {
      for (b_i = 0; b_i < 14; b_i++) {
        t21 = 0.0;
        for (i1 = 0; i1 < 14; i1++) {
          t21 += p_muu_F_j_i[j_L_i_tmp_tmp + (i1 << 3)] * b_LAMBDAUncrt_j_i[i1 +
            14 * b_i];
        }

        coffset = j_L_i_tmp_tmp + (b_i << 3);
        FuT_LAMBDAUncrt_m_Aux_invFx_j_i[coffset] = t21 - Aux_invFx_j_i[coffset];
      }

      coffset = j_L_i_tmp_tmp << 3;
      for (h_L_i_tmp_tmp = 0; h_L_i_tmp_tmp < 8; h_L_i_tmp_tmp++) {
        aoffset = h_L_i_tmp_tmp * 14;
        t17 = 0.0;
        for (i_L_i_tmp_tmp = 0; i_L_i_tmp_tmp < 14; i_L_i_tmp_tmp++) {
          t17 += p_x_muu_j_i[aoffset + i_L_i_tmp_tmp] * Aux_invFx_j_i
            [(i_L_i_tmp_tmp << 3) + j_L_i_tmp_tmp];
        }

        MT_j_i[coffset + h_L_i_tmp_tmp] = t17;
      }
    }

    b_p_muu_F_j_i[0] = (z_i_data[t31_tmp] / (u_i_data[L_i_tmp_tmp] + 10.0) -
                        z_i_data[t31_tmp + 7] / (u_i_data[L_i_tmp_tmp] - 10.0))
      + 0.001;
    b_p_muu_F_j_i[1] = 0.0;
    b_p_muu_F_j_i[2] = 0.0;
    b_p_muu_F_j_i[3] = 0.0;
    b_p_muu_F_j_i[4] = 0.0;
    b_p_muu_F_j_i[5] = 0.0;
    b_p_muu_F_j_i[6] = 0.0;
    b_p_muu_F_j_i[7] = 0.0;
    dmu_u_j_i[0] = 0.0;
    dmu_u_j_i[1] = (z_i_data[t31_tmp + 1] / (u_i_data[b_L_i_tmp_tmp] + 10.0) -
                    z_i_data[t31_tmp + 8] / (u_i_data[b_L_i_tmp_tmp] - 10.0)) +
      0.001;
    dmu_u_j_i[2] = 0.0;
    dmu_u_j_i[3] = 0.0;
    dmu_u_j_i[4] = 0.0;
    dmu_u_j_i[5] = 0.0;
    dmu_u_j_i[6] = 0.0;
    dmu_u_j_i[7] = 0.0;
    b_HAlluT_j_i[0] = 0.0;
    b_HAlluT_j_i[1] = 0.0;
    b_HAlluT_j_i[2] = (z_i_data[t31_tmp + 2] / (u_i_data[c_L_i_tmp_tmp] + 10.0)
                       - z_i_data[t31_tmp + 9] / (u_i_data[c_L_i_tmp_tmp] - 10.0))
      + 0.001;
    b_HAlluT_j_i[3] = 0.0;
    b_HAlluT_j_i[4] = 0.0;
    b_HAlluT_j_i[5] = 0.0;
    b_HAlluT_j_i[6] = 0.0;
    b_HAlluT_j_i[7] = 0.0;
    dv[3] = (z_i_data[t31_tmp + 3] / (u_i_data[d_L_i_tmp_tmp] + 10.0) -
             z_i_data[t31_tmp + 10] / (u_i_data[d_L_i_tmp_tmp] - 10.0)) + 0.001;
    dv1[4] = (z_i_data[t31_tmp + 4] / (u_i_data[e_L_i_tmp_tmp] + 10.0) -
              z_i_data[t31_tmp + 11] / (u_i_data[e_L_i_tmp_tmp] - 10.0)) + 0.001;
    dv2[5] = (z_i_data[t31_tmp + 5] / (u_i_data[f_L_i_tmp_tmp] + 10.0) -
              z_i_data[t31_tmp + 12] / (u_i_data[f_L_i_tmp_tmp] - 10.0)) + 0.001;
    dv3[6] = (z_i_data[t31_tmp + 6] / (u_i_data[g_L_i_tmp_tmp] + 10.0) -
              z_i_data[t31_tmp + 13] / (u_i_data[g_L_i_tmp_tmp] - 10.0)) + 0.001;
    dv4[7] = ((((((((((((((t31 + t32) + t33) + t34) + t35) + t36) + t37) + t45)
                    + t18) + t47) + t19) + t49) + t20) + t51) + z_i_data[t31_tmp
              + 14] / h_L_i_tmp) + 2000.0;

    //  Sensitivities
    for (b_i = 0; b_i < 8; b_i++) {
      for (i1 = 0; i1 < 8; i1++) {
        t21 = 0.0;
        for (h_L_i_tmp_tmp = 0; h_L_i_tmp_tmp < 14; h_L_i_tmp_tmp++) {
          t21 += FuT_LAMBDAUncrt_m_Aux_invFx_j_i[b_i + (h_L_i_tmp_tmp << 3)] *
            p_x_muu_j_i[h_L_i_tmp_tmp + 14 * i1];
        }

        coffset = b_i + (i1 << 3);
        reshapes_f2[coffset] = t21 - MT_j_i[coffset];
      }

      z_i[b_i] = b_p_muu_F_j_i[b_i] + reshapes_f2[b_i];
      z_i[b_i + 8] = dmu_u_j_i[b_i] + reshapes_f2[b_i + 8];
      z_i[b_i + 16] = b_HAlluT_j_i[b_i] + reshapes_f2[b_i + 16];
      z_i[b_i + 24] = dv[b_i] + reshapes_f2[b_i + 24];
      z_i[b_i + 32] = dv1[b_i] + reshapes_f2[b_i + 32];
      z_i[b_i + 40] = dv2[b_i] + reshapes_f2[b_i + 40];
      z_i[b_i + 48] = dv3[b_i] + reshapes_f2[b_i + 48];
      z_i[b_i + 56] = dv4[b_i] + reshapes_f2[b_i + 56];
      for (i1 = 0; i1 < 14; i1++) {
        b_reshapes_f2[i1 + 14 * b_i] = FuT_LAMBDAUncrt_m_Aux_invFx_j_i[b_i + (i1
          << 3)];
      }
    }

    b_inv(z_i, MT_j_i);
    for (b_i = 0; b_i < 8; b_i++) {
      for (i1 = 0; i1 < 14; i1++) {
        coffset = b_i + (i1 << 3);
        p_x_muu_j_i[i1 + 14 * b_i] = -FuT_invFxT_j_i[coffset];
        t21 = 0.0;
        for (h_L_i_tmp_tmp = 0; h_L_i_tmp_tmp < 8; h_L_i_tmp_tmp++) {
          t21 += MT_j_i[b_i + (h_L_i_tmp_tmp << 3)] * b_reshapes_f2[i1 + 14 *
            h_L_i_tmp_tmp];
        }

        p_muu_F_j_i[coffset] = t21;
      }
    }

    for (b_i = 0; b_i < 8; b_i++) {
      for (i1 = 0; i1 < 14; i1++) {
        t21 = 0.0;
        for (h_L_i_tmp_tmp = 0; h_L_i_tmp_tmp < 8; h_L_i_tmp_tmp++) {
          t21 += MT_j_i[b_i + (h_L_i_tmp_tmp << 3)] * p_x_muu_j_i[i1 + 14 *
            h_L_i_tmp_tmp];
        }

        Aux_invFx_j_i[b_i + (i1 << 3)] = t21;
      }
    }

    //  --- LAMBDA = p_lambda_F
    for (b_i = 0; b_i < 14; b_i++) {
      for (i1 = 0; i1 < 14; i1++) {
        t21 = 0.0;
        t19 = 0.0;
        t17 = 0.0;
        t18 = 0.0;
        for (h_L_i_tmp_tmp = 0; h_L_i_tmp_tmp < 8; h_L_i_tmp_tmp++) {
          coffset = b_i + 14 * h_L_i_tmp_tmp;
          aoffset = h_L_i_tmp_tmp + (i1 << 3);
          t21 += b_reshapes_f2[coffset] * Aux_invFx_j_i[aoffset];
          t18 += p_x_muu_j_i[coffset] * Aux_invFx_j_i[aoffset];
          t19 += p_x_muu_j_i[coffset] * p_muu_F_j_i[aoffset];
          t17 += b_reshapes_f2[coffset] * p_muu_F_j_i[aoffset];
        }

        coffset = b_i + 14 * i1;
        aoffset = coffset + 196 * (5 - j);
        p_x_Lambda_i_data[aoffset] = t18;
        p_lambda_Lambda_i_data[aoffset] = t21 + Fx_j_i[coffset];
        p_x_F_i_data[aoffset] = t19 + invFx_j_i[coffset];
        LAMBDA_i_data[aoffset] = t17 - b_LAMBDAUncrt_j_i[coffset];
      }
    }

    //  Coarse Iteration
    for (b_i = 0; b_i < 8; b_i++) {
      t21 = 0.0;
      t19 = 0.0;
      for (i1 = 0; i1 < 14; i1++) {
        h_L_i_tmp_tmp = b_i + (i1 << 3);
        t21 += FuT_invFxT_j_i[h_L_i_tmp_tmp] * lambdaEq_j_i[i1];
        t19 += FuT_LAMBDAUncrt_m_Aux_invFx_j_i[h_L_i_tmp_tmp] * xEq_j_i[i1];
      }

      b_HAlluT_j_i[b_i] = (HAlluT_j_i[b_i] - t21) + t19;
    }

    for (b_i = 0; b_i < 8; b_i++) {
      t21 = 0.0;
      for (i1 = 0; i1 < 8; i1++) {
        t21 += MT_j_i[b_i + (i1 << 3)] * b_HAlluT_j_i[i1];
      }

      dmu_u_j_i[b_i] = t21;
    }

    for (b_i = 0; b_i < 196; b_i++) {
      b_LAMBDAUncrt_j_i[b_i] = -b_LAMBDAUncrt_j_i[b_i];
    }

    for (b_i = 0; b_i < 14; b_i++) {
      t21 = 0.0;
      t19 = 0.0;
      for (i1 = 0; i1 < 14; i1++) {
        h_L_i_tmp_tmp = b_i + 14 * i1;
        t21 += b_LAMBDAUncrt_j_i[h_L_i_tmp_tmp] * xEq_j_i[i1];
        t19 += Fx_j_i[h_L_i_tmp_tmp] * lambdaEq_j_i[i1];
      }

      t17 = 0.0;
      for (i1 = 0; i1 < 8; i1++) {
        t17 += b_reshapes_f2[b_i + 14 * i1] * dmu_u_j_i[i1];
      }

      i1 = b_i + t2_tmp_tmp;
      lambda_i_data[i1] -= (t21 + t19) + t17;
    }

    for (b_i = 0; b_i < 8; b_i++) {
      i1 = b_i + L_i_tmp_tmp;
      u_i_data[i1] -= dmu_u_j_i[b_i];
    }

    //  Recover
    //
    t18 = 0.0;
    for (i_L_i_tmp_tmp = 0; i_L_i_tmp_tmp < 14; i_L_i_tmp_tmp++) {
      t21 = 0.0;
      for (b_i = 0; b_i < 8; b_i++) {
        t21 += p_x_muu_j_i[i_L_i_tmp_tmp + 14 * b_i] * dmu_u_j_i[b_i];
      }

      t19 = 0.0;
      for (b_i = 0; b_i < 14; b_i++) {
        t19 += invFx_j_i[i_L_i_tmp_tmp + 14 * b_i] * xEq_j_i[b_i];
      }

      b_i = i_L_i_tmp_tmp + t2_tmp_tmp;
      x_i_data[b_i] -= t21 + t19;
      std::memcpy(&p_muu_F_i_data[(j * -112 + i_L_i_tmp_tmp * 8) + 560],
                  &p_muu_F_j_i[i_L_i_tmp_tmp * 8], 8U * sizeof(double));
      std::memcpy(&p_muu_Lambda_i_data[(j * -112 + i_L_i_tmp_tmp * 8) + 560],
                  &Aux_invFx_j_i[i_L_i_tmp_tmp * 8], 8U * sizeof(double));
      t17 = std::abs(xEq_j_i[i_L_i_tmp_tmp]);
      if (t17 > t18) {
        t18 = t17;
      }
    }

    KKTxEquation_i_data[5 - j] = t18;
    KKTC_i_data[5 - j] = 0.0;
    t18 = 0.0;
    for (i_L_i_tmp_tmp = 0; i_L_i_tmp_tmp < 8; i_L_i_tmp_tmp++) {
      t17 = std::abs(HAlluT_j_i[i_L_i_tmp_tmp]);
      if (t17 > t18) {
        t18 = t17;
      }
    }

    KKTHu_i_data[5 - j] = t18;
    t18 = 0.0;
    for (i_L_i_tmp_tmp = 0; i_L_i_tmp_tmp < 14; i_L_i_tmp_tmp++) {
      t17 = std::abs(lambdaEq_j_i[i_L_i_tmp_tmp]);
      if (t17 > t18) {
        t18 = t17;
      }
    }

    KKTlambdaEquation_i_data[5 - j] = t18;
  }
}

//
// dx = f(u,x,p)
// Arguments    : const double u[8]
//                const double x[14]
//                double parIdx
//                double f[14]
//                double fu[112]
//                double fx[196]
// Return Type  : void
//
static void f_fu_fx_Wrapper(const double u[8], const double x[14], double parIdx,
  double f[14], double fu[112], double fx[196])
{
  int i;
  double q[7];
  double qd[7];
  double qdd[7];
  double tau[7];
  double dq[49];
  double dqd[49];
  double dtau[49];
  int b_i;
  static const signed char iv[98] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 1 };

  //  dx = f(u,x,p)
  //  Specify your own f(u,x,p) function for code generation
  for (i = 0; i < 7; i++) {
    q[i] = x[i];
    qd[i] = x[i + 7];
    qdd[i] = 0.0;
    tau[i] = u[i];
  }

  qdd_cal(q, qd, qdd, tau, parIdx);
  for (i = 0; i < 7; i++) {
    f[i] = qd[i];
    f[i + 7] = qdd[i];
  }

  std::memset(&fu[0], 0, 112U * sizeof(double));
  for (i = 0; i < 7; i++) {
    q[i] = x[i];
    qd[i] = x[i + 7];
    tau[i] = u[i];
  }

  std::memset(&dq[0], 0, 49U * sizeof(double));
  std::memset(&dqd[0], 0, 49U * sizeof(double));
  std::memset(&dtau[0], 0, 49U * sizeof(double));
  derivatives_cal(q, qd, tau, dq, dqd, dtau, parIdx);
  for (i = 0; i < 7; i++) {
    for (b_i = 0; b_i < 7; b_i++) {
      fu[(b_i + 14 * i) + 7] = dtau[b_i + 7 * i];
    }
  }

  for (i = 0; i < 14; i++) {
    for (b_i = 0; b_i < 7; b_i++) {
      fx[b_i + 14 * i] = iv[b_i + 7 * i];
    }
  }

  for (i = 0; i < 7; i++) {
    for (b_i = 0; b_i < 7; b_i++) {
      int fx_tmp;
      fx_tmp = b_i + 7 * i;
      fx[(b_i + 14 * i) + 7] = dq[fx_tmp];
      fx[(b_i + 14 * (i + 7)) + 7] = dqd[fx_tmp];
    }
  }
}

//
// Arguments    : double varargin_1
//                FILE * *f
//                boolean_T *a
// Return Type  : void
//
static void fileManager(double varargin_1, FILE * *f, boolean_T *a)
{
  signed char fileid;
  fileid = static_cast<signed char>(rt_roundd(varargin_1));
  if ((fileid < 0) || (varargin_1 != fileid)) {
    fileid = -1;
  }

  if (fileid >= 3) {
    *f = eml_openfiles[fileid - 3];
    *a = eml_autoflush[fileid - 3];
  } else if (fileid == 0) {
    *f = stdin;
    *a = true;
  } else if (fileid == 1) {
    *f = stdout;
    *a = true;
  } else if (fileid == 2) {
    *f = stderr;
    *a = true;
  } else {
    *f = NULL;
    *a = true;
  }
}

//
// Arguments    : void
// Return Type  : signed char
//
static signed char filedata()
{
  signed char f;
  int k;
  boolean_T exitg1;
  f = 0;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 20)) {
    if (eml_openfiles[k] == NULL) {
      f = static_cast<signed char>(k + 1);
      exitg1 = true;
    } else {
      k++;
    }
  }

  return f;
}

//
// Arguments    : void
// Return Type  : void
//
static void filedata_init()
{
  FILE * a;
  a = NULL;
  for (int i = 0; i < 20; i++) {
    eml_autoflush[i] = false;
    eml_openfiles[i] = a;
  }
}

//
// Arguments    : const double x[196]
//                double y[196]
// Return Type  : void
//
static void inv(const double x[196], double y[196])
{
  int i;
  double b_x[196];
  int ipiv[14];
  int pipk;
  int k;
  signed char p[14];
  int j;
  int y_tmp;
  int kAcol;
  int b_i;
  for (i = 0; i < 196; i++) {
    y[i] = 0.0;
    b_x[i] = x[i];
  }

  xzgetrf(b_x, ipiv, &pipk);
  for (i = 0; i < 14; i++) {
    p[i] = static_cast<signed char>(i + 1);
  }

  for (k = 0; k < 13; k++) {
    if (ipiv[k] > k + 1) {
      pipk = p[ipiv[k] - 1];
      p[ipiv[k] - 1] = p[k];
      p[k] = static_cast<signed char>(pipk);
    }
  }

  for (k = 0; k < 14; k++) {
    y_tmp = 14 * (p[k] - 1);
    y[k + y_tmp] = 1.0;
    for (j = k + 1; j < 15; j++) {
      i = (j + y_tmp) - 1;
      if (y[i] != 0.0) {
        pipk = j + 1;
        for (b_i = pipk; b_i < 15; b_i++) {
          kAcol = (b_i + y_tmp) - 1;
          y[kAcol] -= y[i] * b_x[(b_i + 14 * (j - 1)) - 1];
        }
      }
    }
  }

  for (j = 0; j < 14; j++) {
    pipk = 14 * j;
    for (k = 13; k >= 0; k--) {
      kAcol = 14 * k;
      i = k + pipk;
      if (y[i] != 0.0) {
        y[i] /= b_x[k + kAcol];
        for (b_i = 0; b_i < k; b_i++) {
          y_tmp = b_i + pipk;
          y[y_tmp] -= y[i] * b_x[b_i + kAcol];
        }
      }
    }
  }
}

//
// Arguments    : const double A[196]
//                double B[14]
// Return Type  : void
//
static void mldivide(const double A[196], double B[14])
{
  double b_A[196];
  int ipiv[14];
  int kAcol;
  int i;
  int k;
  std::memcpy(&b_A[0], &A[0], 196U * sizeof(double));
  xzgetrf(b_A, ipiv, &kAcol);
  for (i = 0; i < 13; i++) {
    if (ipiv[i] != i + 1) {
      double temp;
      temp = B[i];
      B[i] = B[ipiv[i] - 1];
      B[ipiv[i] - 1] = temp;
    }
  }

  for (k = 0; k < 14; k++) {
    kAcol = 14 * k;
    if (B[k] != 0.0) {
      int b_i;
      b_i = k + 2;
      for (i = b_i; i < 15; i++) {
        B[i - 1] -= B[k] * b_A[(i + kAcol) - 1];
      }
    }
  }

  for (k = 13; k >= 0; k--) {
    kAcol = 14 * k;
    if (B[k] != 0.0) {
      B[k] /= b_A[k + kAcol];
      for (i = 0; i < k; i++) {
        B[i] -= B[k] * b_A[i + kAcol];
      }
    }
  }
}

//
// Arguments    : double u
// Return Type  : double
//
static double rt_roundd(double u)
{
  double y;
  if (std::abs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = std::floor(u + 0.5);
    } else if (u > -0.5) {
      y = 0.0;
    } else {
      y = std::ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

//
// Arguments    : double A[196]
//                int ipiv[14]
//                int *info
// Return Type  : void
//
static void xzgetrf(double A[196], int ipiv[14], int *info)
{
  int i;
  int iy;
  int jA;
  int ix;
  for (i = 0; i < 14; i++) {
    ipiv[i] = i + 1;
  }

  *info = 0;
  for (int j = 0; j < 13; j++) {
    int mmj_tmp;
    int b;
    int jj;
    int jp1j;
    double smax;
    int k;
    mmj_tmp = 12 - j;
    b = j * 15;
    jj = j * 15;
    jp1j = b + 2;
    iy = 14 - j;
    jA = 0;
    ix = b;
    smax = std::abs(A[jj]);
    for (k = 2; k <= iy; k++) {
      double s;
      ix++;
      s = std::abs(A[ix]);
      if (s > smax) {
        jA = k - 1;
        smax = s;
      }
    }

    if (A[jj + jA] != 0.0) {
      if (jA != 0) {
        iy = j + jA;
        ipiv[j] = iy + 1;
        ix = j;
        for (k = 0; k < 14; k++) {
          smax = A[ix];
          A[ix] = A[iy];
          A[iy] = smax;
          ix += 14;
          iy += 14;
        }
      }

      i = (jj - j) + 14;
      for (iy = jp1j; iy <= i; iy++) {
        A[iy - 1] /= A[jj];
      }
    } else {
      *info = j + 1;
    }

    iy = b + 14;
    jA = jj;
    for (b = 0; b <= mmj_tmp; b++) {
      smax = A[iy];
      if (A[iy] != 0.0) {
        ix = jj + 1;
        i = jA + 16;
        k = (jA - j) + 28;
        for (jp1j = i; jp1j <= k; jp1j++) {
          A[jp1j - 1] += A[ix] * -smax;
          ix++;
        }
      }

      iy += 14;
      jA += 14;
    }
  }

  if ((*info == 0) && (A[195] == 0.0)) {
    *info = 14;
  }
}

//
// Arguments    : void
// Return Type  : void
//
void Simu_Matlab()
{
  double x0[14];
  int j;
  static double rec_x[112014];
  double tTotal;
  int step;
  double k;
  signed char fileid;
  cell_wrap_7 validatedHoleFilling[1];
  double ref[7];
  static const double b[7] = { 1.5707963267948966, 0.0, 1.5707963267948966, 0.0,
    1.5707963267948966, 0.0, 1.5707963267948966 };

  static const double b_b[7] = { 0.0, 1.5707963267948966, 0.0,
    1.5707963267948966, 0.0, 1.5707963267948966, 0.0 };

  FILE * b_NULL;
  double x0Measured[14];
  double p[168];
  double expl_temp[336];
  double solution_u[192];
  double b_expl_temp[336];
  double c_expl_temp[696];
  static double d_expl_temp[4704];
  double e_expl_temp;
  double output_KKTError;
  double c_output_timeElapsed_searchDire;
  double f_expl_temp;
  double g_expl_temp;
  double output_timeElapsed_total;
  double h_expl_temp;
  double output_iterTotal;
  double i_expl_temp;
  FILE * filestar;
  boolean_T autoflush;
  double b_x0[14];
  char b_validatedHoleFilling[3];
  static double rec_error[8000];
  static double rec_u[64000];
  static double rec_cpuTime[8000];
  static double rec_cpuTimeSearchDirection[8000];
  static double rec_numIter[8000];
  FILE * c_NULL;
  FILE * d_NULL;
  FILE * e_NULL;
  FILE * f_NULL;
  FILE * g_NULL;
  FILE * h_NULL;
  FILE * i_NULL;
  if (!isInitialized_Simu_Matlab) {
    Simu_Matlab_initialize();
  }

  //  For code generation
  //  sampling interval
  //  Load data
  std::memset(&x0[0], 0, 14U * sizeof(double));

  //  define record variables
  for (j = 0; j < 112014; j++) {
    rec_x[j] = 1.0;
  }

  for (j = 0; j < 14; j++) {
    rec_x[8001 * j] = 0.0;
  }

  //  Simulation
  //  Create options for NMPC_Solve
  //     %% options for the initial rho
  //  initial rho
  //  max number of iterations for the initial rho problem
  //  KKT tolerence for the initial rho problem
  //     %% barrier parameter decaying rate
  //     %% options for the target rho
  //  target/end rho
  //  max number of iterations for the initial rho problem
  //  KKT tolerence for the end rho problem
  //     %% line search parameters
  //  enable or disable line search
  //  line search method
  //     %% degree of parallism
  //  1: in serial, otherwise in parallel
  //     %% display
  //     %% check KKT error after iteration
  //  whether to check the KKT error after iteration
  //  degree of parallism: 1 = in serial, otherwise in parallel
  //  do not check the KKT error after iteration
  //  init
  tTotal = 0.0;
  for (step = 0; step < 8000; step++) {
    // simulation steps
    //  Set the reference angle q_ref with a rate limit
    //  0<t<4: q_ref = [0,   pi/2,0,    pi/2, 0,    pi/2, 0]
    //  4<t<8: q_ref = [pi/2,0,   pi/2, 0,    pi/2, 0,    pi/2]
    if (step + 1 < 4000) {
      k = (static_cast<double>(step) + 1.0) * 0.008;
      if (k > 1.0) {
        k = 1.0;
      }

      for (int iRef = 0; iRef < 7; iRef++) {
        ref[iRef] = k * b_b[iRef];
        for (j = 0; j < 24; j++) {
          p[iRef + 7 * j] = ref[iRef];
        }
      }
    } else {
      k = ((static_cast<double>(step) + 1.0) - 4000.0) * 0.008;
      if (k > 1.0) {
        k = 1.0;
      }

      for (int iRef = 0; iRef < 7; iRef++) {
        ref[iRef] = k * b[iRef] + (1.0 - k) * b_b[iRef];
        for (j = 0; j < 24; j++) {
          p[iRef + 7 * j] = ref[iRef];
        }
      }
    }

    //  Solve the optimal control problem
    std::memcpy(&x0Measured[0], &x0[0], 14U * sizeof(double));
    NMPC_Solve(x0, p, expl_temp, solution_u, b_expl_temp, c_expl_temp,
               d_expl_temp, &k, &e_expl_temp, &output_KKTError,
               &c_output_timeElapsed_searchDire, &f_expl_temp, &g_expl_temp,
               &output_timeElapsed_total, &h_expl_temp, &output_iterTotal,
               &i_expl_temp);
    tTotal += output_timeElapsed_total;

    //  Obtain the first optimal control input
    //  System simulation by the 4th-order Explicit Runge-Kutta Method
    std::memcpy(&b_x0[0], &x0[0], 14U * sizeof(double));
    SIM_Plant_RK4(*(double (*)[7])&(*(double (*)[8])&solution_u[0])[0], b_x0, x0);

    //  Record data
    for (j = 0; j < 14; j++) {
      rec_x[(step + 8001 * j) + 1] = x0Measured[j];
    }

    for (j = 0; j < 8; j++) {
      rec_u[step + 8000 * j] = solution_u[j];
    }

    rec_error[step] = output_KKTError;
    rec_cpuTime[step] = output_timeElapsed_total * 1.0E+6;
    rec_cpuTimeSearchDirection[step] = c_output_timeElapsed_searchDire * 1.0E+6;
    rec_numIter[step] = output_iterTotal;
  }

  const char * fmt1;

  //  Log to file
  //  Code generation
  //  show Time Elapsed for RTI
  fmt1 = "Time Elapsed for NMPC_Solve: %f s\r\n";
  printf(fmt1, tTotal);

  //  Log to file
  fileid = cfopen("GEN_log_rec.txt", "wb");

  //  printf header
  validatedHoleFilling[0].f1[0] = 'x';
  validatedHoleFilling[0].f1[2] = '\x00';
  b_NULL = NULL;
  for (j = 0; j < 14; j++) {
    validatedHoleFilling[0].f1[1] = static_cast<signed char>(j + 49);
    fileManager(static_cast<double>(fileid), &filestar, &autoflush);
    if (!(filestar == b_NULL)) {
      b_validatedHoleFilling[0] = validatedHoleFilling[0].f1[0];
      b_validatedHoleFilling[1] = validatedHoleFilling[0].f1[1];
      b_validatedHoleFilling[2] = validatedHoleFilling[0].f1[2];
      fprintf(filestar, "%s\t", b_validatedHoleFilling);
      if (autoflush) {
        fflush(filestar);
      }
    }
  }

  validatedHoleFilling[0].f1[0] = 'u';
  validatedHoleFilling[0].f1[2] = '\x00';
  b_NULL = NULL;
  for (j = 0; j < 8; j++) {
    validatedHoleFilling[0].f1[1] = static_cast<signed char>(j + 49);
    fileManager(static_cast<double>(fileid), &filestar, &autoflush);
    if (!(filestar == b_NULL)) {
      b_validatedHoleFilling[0] = validatedHoleFilling[0].f1[0];
      b_validatedHoleFilling[1] = validatedHoleFilling[0].f1[1];
      b_validatedHoleFilling[2] = validatedHoleFilling[0].f1[2];
      fprintf(filestar, "%s\t", b_validatedHoleFilling);
      if (autoflush) {
        fflush(filestar);
      }
    }
  }

  b_NULL = NULL;
  fileManager(static_cast<double>(fileid), &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "%s\t", "error");
    if (autoflush) {
      fflush(filestar);
    }
  }

  b_NULL = NULL;
  fileManager(static_cast<double>(fileid), &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "%s\t", "numIter");
    if (autoflush) {
      fflush(filestar);
    }
  }

  b_NULL = NULL;
  fileManager(static_cast<double>(fileid), &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "%s\t", "cpuTime");
    if (autoflush) {
      fflush(filestar);
    }
  }

  b_NULL = NULL;
  fileManager(static_cast<double>(fileid), &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "%s\t", "cpuTimeSearchDirection");
    if (autoflush) {
      fflush(filestar);
    }
  }

  b_NULL = NULL;
  fileManager(static_cast<double>(fileid), &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "%s\t", "cpuTimeLineSearch");
    if (autoflush) {
      fflush(filestar);
    }
  }

  b_NULL = NULL;
  fileManager(static_cast<double>(fileid), &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "%s\n", "cpuTimeKKTError");
    if (autoflush) {
      fflush(filestar);
    }
  }

  //  printf data
  b_NULL = NULL;
  c_NULL = NULL;
  d_NULL = NULL;
  e_NULL = NULL;
  f_NULL = NULL;
  g_NULL = NULL;
  h_NULL = NULL;
  i_NULL = NULL;
  for (step = 0; step < 8000; step++) {
    for (j = 0; j < 14; j++) {
      fileManager(static_cast<double>(fileid), &filestar, &autoflush);
      if (!(filestar == b_NULL)) {
        fprintf(filestar, "%f\t", rec_x[step + 8001 * j]);
        if (autoflush) {
          fflush(filestar);
        }
      }
    }

    for (j = 0; j < 8; j++) {
      fileManager(static_cast<double>(fileid), &filestar, &autoflush);
      if (!(filestar == c_NULL)) {
        fprintf(filestar, "%f\t", rec_u[step + 8000 * j]);
        if (autoflush) {
          fflush(filestar);
        }
      }
    }

    fileManager(static_cast<double>(fileid), &filestar, &autoflush);
    if (!(filestar == d_NULL)) {
      fprintf(filestar, "%f\t", rec_error[step]);
      if (autoflush) {
        fflush(filestar);
      }
    }

    fileManager(static_cast<double>(fileid), &filestar, &autoflush);
    if (!(filestar == e_NULL)) {
      fprintf(filestar, "%f\t", rec_numIter[step]);
      if (autoflush) {
        fflush(filestar);
      }
    }

    fileManager(static_cast<double>(fileid), &filestar, &autoflush);
    if (!(filestar == f_NULL)) {
      fprintf(filestar, "%f\t", rec_cpuTime[step]);
      if (autoflush) {
        fflush(filestar);
      }
    }

    fileManager(static_cast<double>(fileid), &filestar, &autoflush);
    if (!(filestar == g_NULL)) {
      fprintf(filestar, "%f\t", rec_cpuTimeSearchDirection[step]);
      if (autoflush) {
        fflush(filestar);
      }
    }

    fileManager(static_cast<double>(fileid), &filestar, &autoflush);
    if (!(filestar == h_NULL)) {
      fprintf(filestar, "%f\t", 0.0);
      if (autoflush) {
        fflush(filestar);
      }
    }

    fileManager(static_cast<double>(fileid), &filestar, &autoflush);
    if (!(filestar == i_NULL)) {
      fprintf(filestar, "%f\n", 0.0);
      if (autoflush) {
        fflush(filestar);
      }
    }
  }

  cfclose(static_cast<double>(fileid));
}

//
// Arguments    : void
// Return Type  : void
//
void Simu_Matlab_initialize()
{
  omp_init_nest_lock(&emlrtNestLockGlobal);
  NMPC_Solve_init();
  filedata_init();

  // user code (Initialize function Body)
  {
    iiwa14_init();
  }

  isInitialized_Simu_Matlab = true;
}

//
// Arguments    : void
// Return Type  : void
//
void Simu_Matlab_terminate()
{
  omp_destroy_nest_lock(&emlrtNestLockGlobal);
  isInitialized_Simu_Matlab = false;
}

//
// File trailer for Simu_Matlab.cpp
//
// [EOF]
//
