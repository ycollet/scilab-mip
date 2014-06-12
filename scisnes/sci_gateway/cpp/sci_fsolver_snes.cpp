////////////////////////////////////////
// Interface to the PETSC/SNES solver //
////////////////////////////////////////

#include <math.h>
#include <string.h>
#include <setjmp.h>
#include <stdio.h>

#define DEBUG 1

extern "C" {
#include <stack-c.h>
#include <sciprint.h>
#include <Scierror.h>
#include <MALLOC.h>
}

#include <api_scilab.h>

#include <parameters.hpp>

#include <mpi.h>
#undef PETSC_HAVE_MPI_WIN_CREATE
#include <petscsnes.h>
#include <petscpc.h>
#include <petscversion.h>

#define _(A) (A)
#define MAX(A,B) ((A<B)?B:A)

#define F_IN          1
#define X_IN          2 
#define PARAM_IN      3
#define X_OUT         4
#define INFO_SNES_OUT 5
#define INFO_KSP_OUT  6
#define ITS_OUT       7
#define LITS_OUT      8
#define FAIL_ITS_OUT  9
#define FAIL_LITS_OUT 10
#define FEV_OUT       11
#define EXTRA_OUT     12
#define FUNC_TMP_1    Rhs+1
#define FUNC_TMP_2    Rhs+2

#if PETSC_VERSION_MAJOR == 3
#define USE_PETSC_3 1
#endif

//#define USE_MPI 1

// TODO:
// Add finite difference jacobian approximation via graph colouring:
// see at http://www.mcs.anl.gov/petsc/petsc-as/snapshots/petsc-current/src/snes/examples/tutorials/ex5s.c.html
// MatFD*

///////////////////////////////////////////////////////
// Prototype of objective functions for petsc / snes //
///////////////////////////////////////////////////////

PetscErrorCode FormJacobian(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
PetscErrorCode Monitor(SNES,PetscInt,PetscReal,void *);
PetscErrorCode SciPetscTraceBackErrorHandler(int,const char *,const char *,const char *,PetscErrorCode,int,const char *,void *);

//////////////////////////
// The Scilab interface //
//////////////////////////

static jmp_buf call_f_env; 
static int sci_obj, lhs_obj, rhs_obj;
static int size_x;
static int size_f;

static char sci_snes_help[] = "Scilab SNES Interface.\n\n";

extern "C" int sci_fsolver_snes(char * fname)
{
  int m_x_in,          n_x_in,          * x_in_addr = NULL;
  int m_x_tmp,         n_x_tmp, l_x_tmp;
  int m_x_out     = 1, n_x_out     = 1;
  int m_info_snes = 1, n_info_snes = 1, l_info_snes,  l_list_info_snes;
  int m_info_ksp  = 1, n_info_ksp  = 1, l_info_ksp,   l_list_info_ksp;
  int m_its       = 1, n_its       = 1, l_its,        l_list_its;
  int m_lits      = 1, n_lits      = 1, l_lits,       l_list_lits;
  int m_fail_its  = 1, n_fail_its  = 1, l_fail_its,   l_list_fail_its;
  int m_fail_lits = 1, n_fail_lits = 1, l_fail_lits,  l_list_fail_lits;
  int m_fev       = 1, n_fev       = 1, l_fev,        l_list_fev;
  int m_jac_tmp   = 1, n_jac_tmp   = 1, l_jac_tmp;
  int * func_1_addr = NULL, * func_2_addr = NULL, * extra_addr = NULL;;
  int argc = 0;
  char ** argv = NULL;
  double tmp_dbl[1];

  int      n, i, j, Index;
  // variables for GET macros
  int      tmp_int, tmp_res = -1;
  double   tmp_double, tmp_double_2, tmp_double_3;
  char   * tmp_char = NULL;
  int      is_sparse = 0, var_type;
  // some return variables
  int m_extra       = 8, n_extra       = 1, l_extra; 
  int m_list_labels = 8, n_list_labels = 1;
  static const char * ListLabels [8] = {"plist","info_snes","info_ksp","its","lits","fail_its","fail_lits","fev"};
  double * x_in = NULL, * x_out = NULL, * x_tmp = NULL, * jac_tmp = NULL;
  int * param_addr = NULL;
  SciSparse SpResult;
  // PETSC Variables
  SNES           snes; // nonlinear solver context
  KSP            ksp;  // linear solver context
  PC             pc;   // preconditioner context
  Vec            x,r;  // solution, residual vectors
  Mat            J;    // Jacobian matrix
  PetscErrorCode ierr;
  PetscScalar   *xx;
  SNESConvergedReason snes_reason;
  KSPConvergedReason  ksp_reason;
  PetscReal      snes_rtol, snes_abstol, snes_stol;
  PetscInt       snes_maxit, snes_maxf, its, lits, fail_its, fail_lits, fev;
  PetscReal      ksp_rtol, ksp_abstol, ksp_dtol;
  PetscInt       ksp_maxits;
  PetscMPIInt    size,rank;
  SciErr         _SciErr;

  argv = NULL;
  argc = 0;

  //
  // Input 1 (objective function)
  //

  // Here, we get a Scilab script function

  if (GetType(F_IN) != sci_c_function)
    {
      Scierror(999,_("%s: Wrong type for input argument #%d: A scilab function expected.\n"), fname, 1);
      return 0;
    }
  
  GetRhsVar(F_IN, EXTERNAL_DATATYPE, &lhs_obj, &rhs_obj, &sci_obj);

  //
  // Input 2 (x vector) *
  //

  // Here, we get the vector which will be sent as a parameter to the function
  // The size of the vector is stored in a static variable: size_x to be used later in the solver function
  _SciErr = getVarAddressFromPosition(pvApiCtx, X_IN, &x_in_addr);
  _SciErr = getMatrixOfDouble(pvApiCtx, x_in_addr, &m_x_in, &n_x_in, &x_in);
  size_x = m_x_in * n_x_in;

  if (n_x_in != 1)
    {
      sciprint(_("%s: x must be vector of dimension n x 1"),fname);
      return 0;
    }

  n       = size_x;  // number of variables
  m_x_out = m_x_in;
  n_x_out = n_x_in;

  //
  // Input 3 parameter list
  //

  init_parameters(PARAM_IN, &param_addr);

  // First, we process the solver option to allocate the good solver

  if (check_parameters(param_addr))
    {
      Scierror(999,_("%s: Argument %d is not a plist\n"),fname,PARAM_IN);
      return 0;
    }

  ///////////////////////////
  // Initialize PETSC/SNES //
  ///////////////////////////

  PetscInitialize(&argc,&argv,(char *)0,sci_snes_help);

#ifdef USE_MPI
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
  if (size!=1) 
    {
      Scierror(999,_("%s: This is a uniprocessor solver only!\n"),fname);
      return 0;
    }
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
#endif

  //
  // Scilab output allocation
  //

  // We create the variable which will store the result of our computation
  _SciErr = allocMatrixOfDouble(pvApiCtx, X_OUT, m_x_out, n_x_out, &x_out);
  
  //
  // First call to the function to get the size of the systems of equations
  //

#ifdef DEBUG
  sciprint("DEBUG: calling the function to get size_f\n");
#endif

  n_x_tmp = n_x_in;
  m_x_tmp = m_x_in;
  CreateVar(FUNC_TMP_1, MATRIX_OF_DOUBLE_DATATYPE, &m_x_tmp, &n_x_tmp, &l_x_tmp);
  //_SciErr = allocMatrixOfDouble(pvApiCtx, FUNC_TMP_1, m_x_tmp, n_x_tmp, &x_tmp);
  //_SciErr = getVarAddressFromPosition(pvApiCtx, FUNC_TMP_1, &func_1_addr);

  memcpy(stk(l_x_tmp), x_in, size_x*sizeof(double));

  n_jac_tmp = 1; m_jac_tmp = 1;
  CreateVar(FUNC_TMP_2, MATRIX_OF_DOUBLE_DATATYPE, &m_jac_tmp, &n_jac_tmp, &l_jac_tmp);
  //_SciErr = allocMatrixOfDouble(pvApiCtx, FUNC_TMP_2, m_jac_tmp, n_jac_tmp, &jac_tmp);

  PExecSciFunction(FUNC_TMP_1, &sci_obj, &lhs_obj, &rhs_obj, "fsolver_snes", call_f_env);

  _SciErr = getVarAddressFromPosition(pvApiCtx, FUNC_TMP_2, &func_2_addr);
  _SciErr = getVarType(pvApiCtx, func_2_addr, &var_type);

  if (var_type==sci_sparse) // Sparse matrix
    {
      is_sparse = 1;
#ifdef DEBUG
      sciprint("DEBUG: a sparse jacobian\n");
#endif
      getAllocatedSparseMatrix(pvApiCtx, func_2_addr, &n_jac_tmp, &m_jac_tmp, &SpResult.nel, &SpResult.mnel, &SpResult.icol, &SpResult.R);
      
      size_f = m_jac_tmp;

#ifdef DEBUG
      sciprint("DEBUG: We have %d variables and %d equations (n_jac_tmp = %d m_jac_tmp = %d)\n", size_x, size_f, n_jac_tmp, m_jac_tmp);
#endif

      if (size_x!=n_jac_tmp)
	{
	  Scierror(999,_("%s: error, you must return a jacobian with %d variables\n"),fname,size_x);
	  freeAllocatedSparseMatrix(SpResult.mnel, SpResult.icol, SpResult.R);
	  return 0;
	}
    }
  else
    {
      is_sparse = 0;
#ifdef DEBUG
      sciprint("DEBUG: a dense jacobian\n");
#endif

      _SciErr = getMatrixOfDouble(pvApiCtx, func_2_addr, &m_jac_tmp, &n_jac_tmp, &jac_tmp);
      
      size_f = m_jac_tmp;

#ifdef DEBUG
      sciprint("DEBUG: We have %d variables and %d equations (n_jac_tmp = %d, m_jac_tmp = %d)\n", size_x, size_f, n_jac_tmp, m_jac_tmp);
#endif

      if (size_x!=n_jac_tmp)
	{
	  Scierror(999,_("%s: error, you must return a jacobian with %d variables\n"),fname,size_x);
	  freeAllocatedSparseMatrix(SpResult.mnel, SpResult.icol, SpResult.R);
	  return 0;
	}
    }

  // Create nonlinear solver context
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  ierr = SNESMonitorSet(snes,Monitor,NULL,0); CHKERRQ(ierr);
  ierr = PetscPushErrorHandler(SciPetscTraceBackErrorHandler,NULL); CHKERRQ(ierr);

  // Create matrix and vector data structures; set corresponding routines
  // Create vectors for solution and nonlinear function

  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,size_x);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&r);CHKERRQ(ierr);
  ierr = VecSetSizes(r,PETSC_DECIDE,size_f);CHKERRQ(ierr);
  ierr = VecSetFromOptions(r);CHKERRQ(ierr);

  // Create Jacobian matrix data structure
  ierr = MatCreate(PETSC_COMM_WORLD,&J);CHKERRQ(ierr);
  //ierr = MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,size_x,size_f);CHKERRQ(ierr);
  ierr = MatSetSizes(J,size_x,size_f,size_x,size_f);CHKERRQ(ierr);
  ierr = MatSetFromOptions(J);CHKERRQ(ierr);

  if (is_sparse)
    {
      ierr = MatSetType(J,MATAIJ); CHKERRQ(ierr);
    }
  else
    {
      ierr = MatSetType(J,MATDENSE); CHKERRQ(ierr);
    }
     
  ierr = SNESSetFunction(snes,r,FormFunction,PETSC_NULL);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,J,J,FormJacobian,PETSC_NULL);CHKERRQ(ierr);

  // PetscErrorCode PETSCSNES_DLLEXPORT SNESSetType(SNES snes,const SNESType type)
  // Sets the method for the nonlinear solver. 
  // SNESLS - Newton's method with line search (systems of nonlinear equations)
  // SNESTR - Newton's method with trust region (systems of nonlinear equations) 
  // SNESPICARD - Picard nonlinear solver that uses successive substitutions  YC: NO MORE PICARD IN SNES !!
  get_string_parameter(param_addr, "snes_type", &tmp_char, &tmp_res, SNESLS, CHECK_VALUES, 2, SNESLS, SNESTR);
  if (tmp_res!=-1)
    {
      ierr = SNESSetType(snes,tmp_char);CHKERRQ(ierr);
    }

  // PetscErrorCode PETSCSNES_DLLEXPORT SNESSetTrustRegionTolerance(SNES snes,PetscReal tol)
  // Sets the trust region parameter tolerance. 
  get_double_parameter(param_addr, "snes_tr_tol", &tmp_double, &tmp_res, PETSC_DEFAULT, CHECK_MIN, 0.0);
  if (tmp_res!=-1)
    {
      ierr = SNESSetTrustRegionTolerance(snes,tmp_double); CHKERRQ(ierr);
    }

  // snes_ls_type Sets the line search routine to be used by the method SNESLS. 
  // 1 - "ls_cubic"         - default line search
  // 2 - "ls_quadratic"     - quadratic line search
  // 3 - "ls_searchno"      - the full Newton step (actually not a line search)
  // 4 - "ls_searchnonorms" - the full Newton step (calculating no norms; faster in parallel) 
  get_int_parameter(param_addr, "snes_ls_type", &tmp_int, &tmp_res, 1, CHECK_VALUES, 4, 1, 2, 3, 4);
  if (tmp_res!=-1)
    {
      switch(tmp_int)
	{
	case 1:
	  ierr = SNESLineSearchSet(snes,SNESLineSearchCubic,NULL);CHKERRQ(ierr);
	  break;
	case 2:
	  ierr = SNESLineSearchSet(snes,SNESLineSearchQuadratic,NULL);CHKERRQ(ierr);
	  break;
	case 3:
	  ierr = SNESLineSearchSet(snes,SNESLineSearchNo,NULL);CHKERRQ(ierr);
	  break;
	case 4:
	  ierr = SNESLineSearchSet(snes,SNESLineSearchNoNorms,NULL);CHKERRQ(ierr);
	  break;
	default:
	  Scierror(999,_("%s: error, you must give 1, 2, 3 or 4 as argument for snes_ls_type\n"),fname);
	  freeAllocatedSparseMatrix(SpResult.mnel, SpResult.icol, SpResult.R);
	  return 0;
	}
    }

  // snes_ls_set_alpha, snes_ls_set_maxstep, snes_ls_steptol:
  // Sets the parameters associated with the line search routine in the Newton-based method SNESLS.
  // PetscErrorCode PETSCSNES_DLLEXPORT SNESLineSearchSetParams(SNES snes,PetscReal alpha,PetscReal maxstep, PetscReal steptol)
  // Input Parameters
  // alpha   - The scalar such that .5*f_{n+1} . f_{n+1} <= .5*f_n . f_n - alpha |f_n . J . f_n|
  // maxstep - The maximum norm of the update vector
  // steptol - The minimum norm fraction of the original step after scaling
  // Note: Pass in PETSC_DEFAULT for any parameter you do not wish to change.
  // We are finding the zero of f() so the one dimensional minimization problem we are solving in the line search is minimize
  // .5*f(x_n + lambda*step_direction) . f(x_n + lambda*step_direction)
  get_double_parameter(param_addr, "snes_ls_set_alpha", &tmp_double, &tmp_res, PETSC_DEFAULT, CHECK_NONE);
  get_double_parameter(param_addr, "snes_ls_set_maxstep", &tmp_double_2, &tmp_res, PETSC_DEFAULT, CHECK_NONE);
#ifdef USE_PETSC_3
  ierr = SNESLineSearchSetParams(snes,tmp_double,tmp_double_2); CHKERRQ(ierr);
#else
  get_double_parameter(param_addr, "snes_ls_set_steptol", &tmp_double_3, &tmp_res, PETSC_DEFAULT, CHECK_MIN, 0.0);
  ierr = SNESLineSearchSetParams(snes,tmp_double,tmp_double_2,tmp_double_3); CHKERRQ(ierr);
#endif

  // PetscErrorCode PETSCSNES_DLLEXPORT SNESSetMaxNonlinearStepFailures(SNES snes, PetscInt maxFails)
  // Sets the maximum number of unsuccessful steps attempted by the nonlinear solver before it gives up. 
  get_int_parameter(param_addr, "snes_max_nl_step_fail", &tmp_int, &tmp_res, 1, CHECK_MIN, 0);
  if (tmp_res!=-1)
    {
      ierr = SNESSetMaxNonlinearStepFailures(snes, tmp_int); CHKERRQ(ierr);
    }

  // PetscErrorCode PETSCSNES_DLLEXPORT SNESSetMaxLinearSolveFailures(SNES snes, PetscInt maxFails)
  // the number of failed linear solve attempts allowed before SNES returns with a diverged reason of SNES_DIVERGED_LINEAR_SOLVE 
  get_int_parameter(param_addr, "snes_max_lin_solve_fail", &tmp_int, &tmp_res, 0, CHECK_MIN, 0);
  if (tmp_res!=-1)
    {
      ierr = SNESSetMaxLinearSolveFailures(snes, tmp_int); CHKERRQ(ierr);
    }

#ifdef USE_PETSC_3
  // PetscErrorCode PETSCSNES_DLLEXPORT SNESSetLagPreconditioner(SNES snes,PetscInt lag)
  // Determines when the preconditioner is rebuilt in the nonlinear solve. 
  // lag 	
  //   - -1 indicates NEVER rebuild, 
  //   -  1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 
  //   -  2 means every second time the Jacobian is built etc. 
  //   - -2 indicates rebuild preconditioner at next chance but then never rebuild after that 
  get_int_parameter(param_addr, "snes_lag_precond", &tmp_int, &tmp_res, 1, CHECK_VALUES, 4, -1, -2, 1, 2);
  if (tmp_res!=-1)
    {
      ierr = SNESSetLagPreconditioner(snes,tmp_int); CHKERRQ(ierr);
    }

  // PetscErrorCode PETSCSNES_DLLEXPORT SNESSetLagJacobian(SNES snes,PetscInt lag)
  // Determines when the Jacobian is rebuilt in the nonlinear solve. 
  // See SNESSetLagPreconditioner() for determining how often the preconditioner is rebuilt. 
  get_int_parameter(param_addr, "snes_lag_jac", &tmp_int, &tmp_res, 1, CHECK_MIN, 0);
  if (tmp_res!=-1)
    {
      ierr = SNESSetLagJacobian(snes,tmp_int); CHKERRQ(ierr);
    }
#endif

  // Set linear solver defaults for this problem. By extracting the
  // KSP, KSP, and PC contexts from the SNES context, we can then
  // directly call any KSP, KSP, and PC routines to set various options.

  ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  //ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE); CHKERRQ(ierr);

  // PCSetType can be PCNONE, PCJACOBI, PCILU, or PCBJACOBI
  get_string_parameter(param_addr, "snes_pc_type", &tmp_char, &tmp_res, PCJACOBI, CHECK_VALUES, 4, PCNONE, PCJACOBI, PCILU, PCBJACOBI);
  if (tmp_res!=-1)
    {
      ierr = PCSetType(pc,tmp_char);CHKERRQ(ierr);
    }

  // Builds KSP for a particular solver. 
  // KSPRICHARDSON "richardson"
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPRichardsonSetScale(KSP ksp,PetscReal scale)
  // KSPCHEBYCHEV  "chebychev"
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPChebychevSetEigenvalues(KSP ksp,PetscReal emax,PetscReal emin)
  // KSPCG         "cg"
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPCGSetType(KSP ksp,KSPCGType type)
  //   KSPCGNE       "cgne"
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPCGSetType(KSP ksp,KSPCGType type)
  //   KSPNASH       "nash"
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPNASHSetRadius(KSP ksp, PetscReal radius)
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPNASHGetNormD(KSP ksp, PetscReal *norm_d)
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPNASHGetObjFcn(KSP ksp, PetscReal *o_fcn)
  //   KSPSTCG       "stcg"
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPSTCGSetRadius(KSP ksp, PetscReal radius)
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPSTCGGetNormD(KSP ksp, PetscReal *norm_d)
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPSTCGGetObjFcn(KSP ksp, PetscReal *o_fcn)
  //   KSPGLTR       "gltr"
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPGLTRSetRadius(KSP ksp, PetscReal radius)
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPGLTRGetNormD(KSP ksp, PetscReal *norm_d)
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPGLTRGetObjFcn(KSP ksp, PetscReal *o_fcn)
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPGLTRGetMinEig(KSP ksp, PetscReal *e_min)
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPGLTRGetLambda(KSP ksp, PetscReal *lambda)
  // KSPGMRES      "gmres"
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPGMRESSetRestart(KSP ksp, PetscInt restart) 
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPGMRESSetHapTol(KSP ksp,PetscReal tol)
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPGMRESSetPreAllocateVectors(KSP ksp)
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPGMRESSetCGSRefinementType(KSP ksp,KSPGMRESCGSRefinementType type)
  //   KSPFGMRES     "fgmres" 
  // Idem GMRES
  //   KSPLGMRES     "lgmres"
  // Idem GMRES
  // KSPTCQMR      "tcqmr"
  // KSPBCGS       "bcgs"
  // KSPIBCGS        "ibcgs"
  // KSPBCGSL        "bcgsl"
  // KSPCGS        "cgs"
  // KSPTFQMR      "tfqmr"
  // KSPCR         "cr"
  // KSPLSQR       "lsqr"
  // KSPPREONLY    "preonly"
  // KSPQCG        "qcg"
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPQCGSetTrustRegionRadius(KSP ksp,PetscReal delta)
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPQCGGetTrialStepNorm(KSP ksp,PetscReal *tsnorm)
  // PetscErrorCode PETSCKSP_DLLEXPORT KSPQCGGetQuadratic(KSP ksp,PetscReal *quadratic)
  // KSPBICG       "bicg"
  // KSPMINRES     "minres"
  // KSPSYMMLQ     "symmlq"
  // KSPLCD        "lcd"
  get_string_parameter(param_addr, "snes_ksp_type", &tmp_char, &tmp_res, "gmres", CHECK_VALUES, 24, 
		       "richardson", "chebychev", "cg", "cgne", "nash", "stcg", "gltr", "gmres",
		       "fgmres", "lgmres", "tcqmr", "bcgs", "ibcgs", "bcgsl", "cgs", "tfqmr",
		       "cr", "lsqr", "preonly", "qcg", "bicg", "minres", "symmlq", "lcd");
  if (tmp_res!=-1)
    {
      ierr = KSPSetType(ksp, tmp_char); CHKERRQ(ierr);
    }

  // PetscErrorCode PETSCKSP_DLLEXPORT KSPSetTolerances(KSP ksp,
  //                                                    PetscReal rtol,
  //                                                    PetscReal abstol,
  //                                                    PetscReal dtol,
  //                                                    PetscInt maxits)
  // rtol   - the relative convergence tolerance (relative decrease in the residual norm)
  // abstol - the absolute convergence tolerance (absolute size of the residual norm)
  // dtol   - the divergence tolerance (amount residual can increase before KSPDefaultConverged() concludes that the method is diverging)
  // maxits - maximum number of iterations to use 
  get_double_parameter(param_addr, "snes_ksp_rtol",   &ksp_rtol,   &tmp_res, PETSC_DEFAULT, CHECK_MIN, 0.0);
  get_double_parameter(param_addr, "snes_ksp_abstol", &ksp_abstol, &tmp_res, PETSC_DEFAULT, CHECK_MIN, 0.0);
  get_double_parameter(param_addr, "snes_ksp_dtol",   &ksp_dtol,   &tmp_res, PETSC_DEFAULT, CHECK_MIN, 0.0);
  get_int_parameter(param_addr,    "snes_ksp_maxits", &ksp_maxits, &tmp_res, 20,            CHECK_MIN, 0);

  ierr = KSPSetTolerances(ksp,ksp_rtol,ksp_abstol,ksp_dtol,ksp_maxits);CHKERRQ(ierr);

  // KSPSetComputeSingularValues
  // Sets a flag so that the extreme singular values will be calculated via a Lanczos or Arnoldi process as the linear system is solved.
  get_int_parameter(param_addr, "snes_ksp_comp_sing_val", &tmp_int, &tmp_res, PETSC_TRUE, CHECK_VALUES, 2, PETSC_TRUE, PETSC_FALSE);
  if (tmp_res!=-1)
    {
      ierr = KSPSetComputeSingularValues(ksp, PETSC_TRUE); CHKERRQ(ierr);
    }

  // Set SNES/KSP/KSP/PC runtime options, e.g., -snes_view -snes_monitor -ksp_type <ksp> -pc_type <pc>
  // These options will override those specified above as long as SNESSetFromOptions() is called _after_ any other customization routines.

  // PetscErrorCode PETSCSNES_DLLEXPORT SNESSetTolerances(SNES snes,
  //                                                      PetscReal abstol,
  //                                                      PetscReal rtol,
  //                                                      PetscReal stol,
  //                                                      PetscInt maxit,
  //                                                      PetscInt maxf)
  // abstol - absolute convergence tolerance
  // rtol   - relative convergence tolerance
  // stol   - convergence tolerance in terms of the norm of the change in the solution between steps
  // maxit  - maximum number of iterations
  // maxf   - maximum number of function evaluations
  get_double_parameter(param_addr, "snes_rtol",   &snes_rtol,   &tmp_res, PETSC_DEFAULT, CHECK_MIN, 0.0);
  get_double_parameter(param_addr, "snes_abstol", &snes_abstol, &tmp_res, PETSC_DEFAULT, CHECK_MIN, 0.0);
  get_double_parameter(param_addr, "snes_stol",   &snes_stol,   &tmp_res, PETSC_DEFAULT, CHECK_MIN, 0.0);
  get_int_parameter(param_addr, "snes_maxit", &snes_maxit, &tmp_res, 20, CHECK_MIN, 0);
  get_int_parameter(param_addr, "snes_maxf",  &snes_maxf,  &tmp_res, 20, CHECK_MIN, 0);

  ierr = SNESSetTolerances(snes,snes_abstol,snes_rtol,snes_stol,snes_maxit,snes_maxf);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  // Store the initial solution in x
  VecGetArray(x,&xx);
  for(i=0;i<size_x;i++) xx[i] = x_in[i];
  VecRestoreArray(x,&xx);

  // Now solve
  ierr = SNESSolve(snes,PETSC_NULL,x);CHKERRQ(ierr);
  
  // Get the solution and store it so as to be returned by Scilab
  VecGetArray(x,&xx);
  memcpy(x_out, xx, size_x*sizeof(double));
  VecRestoreArray(x,&xx);

  // Get status informations
  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
  ierr = SNESGetConvergedReason(snes,&snes_reason);CHKERRQ(ierr);
  ierr = KSPGetConvergedReason(ksp,&ksp_reason); CHKERRQ(ierr);
  ierr = SNESGetLinearSolveIterations(snes,&lits); CHKERRQ(ierr);
#ifdef USE_PETSC_3
  ierr = SNESGetNumberFunctionEvals(snes, &fev); CHKERRQ(ierr);
#else
  fev = -1;
#endif
  ierr = SNESGetLinearSolveFailures(snes, &fail_lits); CHKERRQ(ierr);
  ierr = SNESGetNonlinearStepFailures(snes, &fail_its); CHKERRQ(ierr);

#ifdef DEBUG
  sciprint("DEBUG: computing and printing r\n");
  SNESComputeFunction(snes,x,r);
  MatStructure flg;
  SNESComputeJacobian(snes,x,&J,&J,&flg);

  sciprint("DEBUG: After SNESSolve - computing and printing J\n");
  MatView(J,PETSC_VIEWER_STDOUT_SELF);
  sciprint("DEBUG: After SNESSolve - computing and printing r\n");
  VecView(r,PETSC_VIEWER_STDOUT_SELF);
#endif

#ifdef DEBUG
  SNESView(snes,PETSC_VIEWER_STDOUT_SELF);
  KSPView(ksp,PETSC_VIEWER_STDOUT_SELF);
#endif

  ///////////////////////
  // SNES return codes //
  ///////////////////////
  // converged
  // SNES_CONVERGED_FNORM_ABS         =  2, ||F|| < atol
  // SNES_CONVERGED_FNORM_RELATIVE    =  3, ||F|| < rtol*||F_initial||
  // SNES_CONVERGED_PNORM_RELATIVE    =  4, Newton computed step size small; || delta x || < tol
  // SNES_CONVERGED_ITS               =  5, maximum iterations reached
  // SNES_CONVERGED_TR_DELTA          =  7,
  // diverged
  // SNES_DIVERGED_FUNCTION_DOMAIN    = -1,  
  // SNES_DIVERGED_FUNCTION_COUNT     = -2,  
  // SNES_DIVERGED_LINEAR_SOLVE       = -3, 
  // SNES_DIVERGED_FNORM_NAN          = -4, 
  // SNES_DIVERGED_MAX_IT             = -5, means that the solver reached the maximum number of iterations without satisfying any convergence criteria. 
  //                                        SNES_CONVERGED_ITS means that SNESSkipConverged() was chosen as the convergence test; 
  //                                        thus the usual convergence criteria have not been checked and may or may not be satisfied. 
  // SNES_DIVERGED_LS_FAILURE         = -6,
  // SNES_DIVERGED_LOCAL_MIN          = -8  || J^T b || is small, implies converged to local minimum of F()
  //                                        This can only occur when using the line-search variant of SNES. 
  //                                        The line search wants to minimize Q(alpha) = 1/2 || F(x + alpha s) ||^2_2 this occurs at 
  //                                        Q'(alpha) = s^T F'(x+alpha s)^T F(x+alpha s) = 0. If s is the Newton direction - F'(x)^(-1)F(x) 
  //                                        then you get Q'(alpha) = -F(x)^T F'(x)^(-1)^T F'(x+alpha s)F(x+alpha s); 
  //                                        when alpha = 0 Q'(0) = - ||F(x)||^2_2 which is always NEGATIVE if F'(x) is invertible. 
  //                                        This means the Newton direction is a descent direction and the line search should succeed 
  //                                        if alpha is small enough.
  //                                        If F'(x) is NOT invertible AND F'(x)^T F(x) = 0 then Q'(0) = 0 and the Newton direction is NOT a
  //                                        descent direction so the line search will fail. 
  //                                        All one can do at this point is change the initial guess and try again.
  //                                        An alternative explanation: Newton's method can be regarded as replacing the function with its
  //                                        linear approximation and minimizing the 2-norm of that. 
  //                                        That is F(x+s) approx F(x) + F'(x)s so we minimize || F(x) + F'(x) s ||^2_2; do this using Least Squares. 
  //                                        If F'(x) is invertible then s = - F'(x)^(-1)F(x) otherwise F'(x)^T F'(x) s = -F'(x)^T F(x). 
  //                                        If F'(x)^T F(x) is NOT zero then there exists a nontrival (that is F'(x)s != 0) solution to the
  //                                        equation and this direction is s = - [F'(x)^T F'(x)]^(-1) F'(x)^T F(x) 
  //                                        so Q'(0) = - F(x)^T F'(x) [F'(x)^T F'(x)]^(-T) F'(x)^T 
  //                                        F(x) = - (F'(x)^T F(x)) [F'(x)^T F'(x)]^(-T) (F'(x)^T F(x)). 
  //                                        Since we are assuming (F'(x)^T F(x)) != 0 and F'(x)^T F'(x) has no negative eigenvalues Q'(0) < 0 
  //                                        so s is a descent direction and the line search should succeed for small enough alpha.
  //                                        Note that this RARELY happens in practice. Far more likely the linear system is not being 
  //                                        solved (well enough?) or the Jacobian is wrong. 
  // SNES_CONVERGED_ITERATING         =  0

  //////////////////////
  // KSP return codes //
  //////////////////////
  // Converged
  // KSP_CONVERGED_RTOL               =  2,
  // KSP_CONVERGED_ATOL               =  3,
  // KSP_CONVERGED_ITS                =  4,
  // KSP_CONVERGED_CG_NEG_CURVE       =  5,
  // KSP_CONVERGED_CG_CONSTRAINED     =  6,
  // KSP_CONVERGED_STEP_LENGTH        =  7,
  // KSP_CONVERGED_HAPPY_BREAKDOWN    =  8,
  // Diverged
  // KSP_DIVERGED_NULL                = -2,
  // KSP_DIVERGED_ITS                 = -3,
  // KSP_DIVERGED_DTOL                = -4,
  // KSP_DIVERGED_BREAKDOWN           = -5,
  // KSP_DIVERGED_BREAKDOWN_BICG      = -6,
  // KSP_DIVERGED_NONSYMMETRIC        = -7,
  // KSP_DIVERGED_INDEFINITE_PC       = -8,
  // KSP_DIVERGED_NAN                 = -9,
  // KSP_DIVERGED_INDEFINITE_MAT      = -10,
  
  // KSP_CONVERGED_ITERATING          =  0

  // Create the 'status' structure of type plist
  _SciErr = createMList(pvApiCtx, EXTRA_OUT, 8, &extra_addr);
  _SciErr = createMatrixOfStringInList(pvApiCtx, EXTRA_OUT, extra_addr, 1, m_list_labels, n_list_labels, (char **)ListLabels);

  // the status of the snes solver
  tmp_dbl[0] = snes_reason;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, EXTRA_OUT, extra_addr, 2, 1, 1, tmp_dbl);

  // the status of the ksp solver
  tmp_dbl[0] = ksp_reason;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, EXTRA_OUT, extra_addr, 3, 1, 1, tmp_dbl);
  
  // the number of nonlinear iterations performed
  tmp_dbl[0]  = its;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, EXTRA_OUT, extra_addr, 4, 1, 1, tmp_dbl);

  // the number of linear iterations performed
  tmp_dbl[0] = lits;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, EXTRA_OUT, extra_addr, 5, 1, 1, tmp_dbl);

  // the number of failed nonlinear iterations performed
  tmp_dbl[0] = fail_its;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, EXTRA_OUT, extra_addr, 6, 1, 1, tmp_dbl);

  // the number of failed linear iterations performed
  tmp_dbl[0] = fail_lits;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, EXTRA_OUT, extra_addr, 7, 1, 1, tmp_dbl);

  // the number of function evaluations
  tmp_dbl[0] = fev;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, EXTRA_OUT, extra_addr, 8, 1, 1, tmp_dbl);

#ifdef DEBUG
  sciprint("DEBUG: its = %d lits = %d snes_reason = %d ksp_reason = %d\n", its, lits, snes_reason, ksp_reason);
#endif

  // Allocation of the return list

  LhsVar(1) = X_OUT;
  LhsVar(2) = EXTRA_OUT;

  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = VecDestroy(r);CHKERRQ(ierr);
  ierr = MatDestroy(J);CHKERRQ(ierr); 
  ierr = SNESDestroy(snes);CHKERRQ(ierr);
  ierr = PetscPopErrorHandler();CHKERRQ(ierr);

  if (is_sparse) freeAllocatedSparseMatrix(SpResult.mnel, SpResult.icol, SpResult.R);

#ifdef USE_MPI
  ierr = PetscFinalize();CHKERRQ(ierr);
#endif

  return 0;
}

//
// The functions for PETSC / SNES
//

PetscErrorCode FormFunction(SNES snes, Vec x_in, Vec f,void *dummy)
{
  PetscScalar *xx, *ff;
  PetscErrorCode ierr;
  int m_x   = size_x, n_x   = 1, * x_addr = NULL, l_x;
  int m_jac = 1,      n_jac = 1, * jac_addr = NULL, l_jac;
  int rhs_old = Rhs, nbvars_old = Nbvars;
  double * x = NULL, * jac = NULL;
  int i;
  SciErr _SciErr;

  Rhs    = FUNC_TMP_1 + MAX(rhs_obj,lhs_obj);
  Nbvars = FUNC_TMP_1 + MAX(rhs_obj,lhs_obj);

#ifdef DEBUG
  sciprint("DEBUG: calling FormFunction\n");
  PetscInt DbgSizeX, DbgSizeF;
  VecGetSize(x_in,&DbgSizeX);
  VecGetSize(f,&DbgSizeF);
  sciprint("size of x = %d size of f = %d\n", DbgSizeX, DbgSizeF);
#endif

  // Get pointers to vector data.
  // - For default PETSc vectors, VecGetArray() returns a pointer to the data array.  Otherwise, the routine is implementation dependent.
  // - You MUST call VecRestoreArray() when you no longer need access to the array.
  VecGetArray(x_in, &xx);
  VecGetArray(f,    &ff);

  CreateVar(FUNC_TMP_1, MATRIX_OF_DOUBLE_DATATYPE, &m_x, &n_x, &l_x);
  memcpy(stk(l_x), xx, size_x*sizeof(double));

  CreateVar(FUNC_TMP_2, MATRIX_OF_DOUBLE_DATATYPE, &m_jac, &n_jac, &l_jac);

  PExecSciFunction(FUNC_TMP_1, &sci_obj, &lhs_obj, &rhs_obj, "FormFunction", call_f_env);

  GetRhsVar(FUNC_TMP_1, MATRIX_OF_DOUBLE_DATATYPE, &m_x, &n_x, &l_x);

  if (m_x*n_x!=size_f)
    {
      Scierror(999,_("fsolver_snes: error, you must return %d values\n"),size_f);
      Rhs    = rhs_old;
      Nbvars = nbvars_old;
      return 0;
    }

  memcpy(ff, stk(l_x), size_f*sizeof(double));

  VecRestoreArray(x_in, &xx);
  VecRestoreArray(f,    &ff);

#ifdef DEBUG
  sciprint("DEBUG: After FormFunction - computing and printing x\n");
  VecView(x_in,PETSC_VIEWER_STDOUT_SELF);
  sciprint("DEBUG: After FormFunction - computing and printing f\n");
  VecView(f,PETSC_VIEWER_STDOUT_SELF);
#endif

  Rhs    = rhs_old;
  Nbvars = nbvars_old;

  return 0;
}

PetscErrorCode FormJacobian(SNES snes, Vec x_in, Mat *jac_in, Mat *B, MatStructure *flag, void *dummy)
{
  PetscScalar A, *xx;
  PetscErrorCode ierr;
  int m_x   = size_x, n_x   = 1, * x_addr = NULL, l_x;
  int m_jac = 1,      n_jac = 1, * jac_addr = NULL, l_jac;
  double * x = NULL, * jac = NULL;
  int i, j, Index, irow, jcol, var_type, is_sparse = 0;
  int rhs_old = Rhs, nbvars_old = Nbvars;
  SciSparse SpResult;
  SciErr _SciErr;

  Rhs    = FUNC_TMP_1 + MAX(rhs_obj,lhs_obj);
  Nbvars = FUNC_TMP_1 + MAX(rhs_obj,lhs_obj);

#ifdef DEBUG
  sciprint("DEBUG: calling FormJacobian\n");
  PetscInt DbgSizeX;
  PetscInt DbgM, DbgN;
  VecGetSize(x_in,&DbgSizeX);
  MatGetSize(*jac_in,&DbgM,&DbgN);
  sciprint("size of x = %d\n", DbgSizeX);
  sciprint("size of jac: %d %d\n", DbgM, DbgN);
#endif

  VecGetArray(x_in,&xx);

  CreateVar(FUNC_TMP_1, MATRIX_OF_DOUBLE_DATATYPE, &m_x, &n_x, &l_x);
  memcpy(stk(l_x), xx, size_x*sizeof(double));

  CreateVar(FUNC_TMP_2, MATRIX_OF_DOUBLE_DATATYPE, &m_jac, &n_jac, &l_jac);

  PExecSciFunction(FUNC_TMP_1, &sci_obj, &lhs_obj, &rhs_obj, "FormJacobian", call_f_env);

  _SciErr = getVarAddressFromPosition(pvApiCtx, FUNC_TMP_2, &jac_addr);
  _SciErr = getVarType(pvApiCtx, jac_addr, &var_type);

  if (var_type==sci_sparse)
    {
      is_sparse = 1;

      getAllocatedSparseMatrix(pvApiCtx, jac_addr, &m_jac, &n_jac, &SpResult.nel, &SpResult.mnel, &SpResult.icol, &SpResult.R);

      if ((m_jac!=size_f)&&(n_jac!=size_x))
	{
	  Scierror(999,_("fsolver_snes: error, you must return a %d * %d sparse matrix\n"),size_f,size_x);
	  Rhs    = rhs_old;
	  Nbvars = nbvars_old;
	  return 0;
	}

      Index = 0;
      for(i=0;i<SpResult.m;i++)
	{
	  for(j=0;j<SpResult.mnel[i];j++)
	    {
	      A = SpResult.R[Index];
	      irow = i;
	      jcol = SpResult.icol[Index]-1;
	      Index++;
#ifdef DEBUG
	      sciprint("DEBUG: i = %d, j = %d, val = %f\n", irow, jcol, A);
#endif
	      ierr = MatSetValues(*B,1,&irow,1,&jcol,&A,INSERT_VALUES);CHKERRQ(ierr);
	    }
	}
    }
  else
    {
      _SciErr = getMatrixOfDouble(pvApiCtx, jac_addr, &m_jac, &n_jac, &jac); 
      
      if ((m_jac!=size_f)&&(n_jac!=size_x))
	{
	  Scierror(999,_("fsolver_snes: error, you must return a %d * %d full matrix\n"),size_f,size_x);
	  Rhs    = rhs_old;
	  Nbvars = nbvars_old;
	  return 0;
	}
      
      for(i=0;i<m_jac;i++)
	{
	  for(j=0;j<n_jac;j++)
	    {
	      A = *(jac+i+j*m_jac);
	      irow = i;
	      jcol = j;
	  
	      ierr = MatSetValues(*B,1,&irow,1,&jcol,&A,INSERT_VALUES);CHKERRQ(ierr);
	    }
	}
    }


  // SAME_NONZERO_PATTERN
  // DIFFERENT_NONZERO_PATTERN
  // SAME_PRECONDITIONER 
  // SUBSET_NONZERO_PATTERN
  //*flag = DIFFERENT_NONZERO_PATTERN;
  *flag = DIFFERENT_NONZERO_PATTERN;

  // Restore the array
  VecRestoreArray(x_in,&xx);

  // Assemble matrix

  // First: preconditionner
  ierr = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  // Second: if the preconditionner if different from the jacobian, then, we process it
  if (*jac_in!=*B)
    {
      ierr = MatAssemblyBegin(*jac_in,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(*jac_in,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }

#ifdef DEBUG
  sciprint("DEBUG: After FormJacobian - computing and printing jac\n");
  MatView(*jac_in,PETSC_VIEWER_STDOUT_SELF);
  sciprint("DEBUG: After FormJacobian - computing and printing B\n");
  MatView(*B,PETSC_VIEWER_STDOUT_SELF);
#endif

  if (is_sparse) freeAllocatedSparseMatrix(SpResult.mnel, SpResult.icol, SpResult.R);

  Rhs    = rhs_old;
  Nbvars = nbvars_old;

  return 0;
}

PetscErrorCode Monitor(SNES snes,PetscInt its,PetscReal fnorm,void *ctx)
{
  sciprint(_("fsolver_snes: iter = %d, SNES Function norm %f\n"),its,fnorm);
  return 0;
}

PetscErrorCode SciPetscTraceBackErrorHandler(int line,const char *fun,const char* file,const char *dir,PetscErrorCode n,int p,const char *mess,void *ctx)
{
  sciprint(_("fsolver_snes - error: line = %d, function %s, file %s, dir %s. Error %d\n"),line, fun, file, dir, n);
  sciprint(_("fsolver_snes - error: Specific error %d, message %s\n"),p,mess);
  return 0;
}
