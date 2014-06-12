#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <setjmp.h>

#include "ds_common.h"

#include <stack-c.h>
#include <sciprint.h>
#include <Scierror.h>
#include <sci_types.h>
#include <MALLOC.h>

#include <api_scilab.h>
#include <api_parameters.h>
#include <freeArrayOfString.h>

#include "func_utils.h"

// #define DEBUG 1

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

#define FILTERSD_ERROR if(_SciErr.iErr)		\
    {                                           \
      printError(&_SciErr, 0);                  \
      return _SciErr.iErr;                      \
    }

#define FILTERSD_ERROR_NORETURN if(_SciErr.iErr)	\
    {							\
      printError(&_SciErr, 0);				\
      return;						\
    }

#define FILTERSD_ERROR_RETURN_NULL if(_SciErr.iErr)	\
    {							\
      printError(&_SciErr, 0);				\
      return NULL;					\
    }

#define X_IN          1
#define FOBJ_IN       2
#define GRAD_IN       3
#define X_LOWER_IN    4
#define X_UPPER_IN    5
#define PARAMS_IN     6
#define LAST_PARAM_IN 6

#define X_OUT          (Rhs+1)
#define F_OUT          (Rhs+2)
#define IFAIL_OUT      (Rhs+3)
#define LAMBDA_OUT     (Rhs+4)
#define PARAMS_OUT     (Rhs+5)
#define LAST_PARAM_OUT (Rhs+5)

void C2F(dense_functions)(int * n, int * m, double * x, double * f, double * c, double * user, int * iuser);
void C2F(dense_gradients)(int * n, int * m, double * x, double * a, double * user, int * iuser);

struct param_obj {
  int fobj_type; // 0 for a Scilab function, 1 for a C function
  int sci_obj;
  int lhs_obj;
  int rhs_obj;
  int stack_pos;
  int n_x;
  int n_c;
  int is_sparse;
  fobj_filtersd_ptr function;
  grad_filtersd_ptr gradient;
};

/*************************
 * Definition of FTables *
 *************************/

static FTAB FTab_filtersd_function[] = {{(char *) 0, (voidf) 0}};

/************************
 * The Scilab interface *
 ************************/

// static jmp_buf nlopt_function_env; 

static struct param_obj param_fobj;
static struct param_obj param_grad;

// Parameters in:
// - "rho"
// - "htol"
// - "rgtol"
// - "fmin"
// - "maxit"
// - "iprint"
// - "kmax"
// - "maxg"

// Parameters out:
// - "dnorm"
// - "h"
// - "hJt"
// - "hJ"
// - "ipeq"
// - "k"
// - "itn"
// - "nft"
// - "ngt"

int sci_filtersd_dense(char * fname)
{
  int n_x_in = 0,        m_x_in = 0,        * pi_x_in        = NULL; double * x_in       = NULL;
  int n_x_upper_in = 0,  m_x_upper_in = 0,  * pi_x_upper_in  = NULL; double * x_upper_in = NULL;
  int n_x_lower_in = 0,  m_x_lower_in = 0,  * pi_x_lower_in  = NULL; double * x_lower_in = NULL;
  int n_x_out = 0,       m_x_out = 0;      double * x_out      = NULL;
  int n_f_out = 0,       m_f_out = 0;      double * f_out      = NULL;
  int n_lambda_out = 0,  m_lambda_out = 0; double * lambda_out = NULL;
  int n_ifail_out = 0,   m_ifail_out = 0;  int ifail           = 0;
  int n_tmp_var = 0, m_tmp_var = 0, * pi_tmp_var = NULL; double * tmp_var = NULL;
  int  * param_in_addr  = NULL;
  int  * param_out_addr = NULL;
  char * cstype         = NULL;
  double rho    = 100.0;
  double htol   = 1.0e-6;
  double rgtol  = 1.0e-4;
  double fmin   = -1.0e6;
  double ainfty = 1e20;
  double ubd    = 1e5;
  int    mlp    = 50;
  int    mxf    = 50;
  int    mxgr   = 1e6;
  int    mxm1   = 0;
  int    mxws   = 30000;
  int    mxlws  = 30000;
  int    maxit  = 100;
  int    iprint = 0;
  int    kmax   = 1;
  int    maxg   = 5;
  int    Log    = 0;
  double tmp_double = 0.0, f_tmp_out = 0.0;
  int tmp_int = 0, tmp_res = 0, grad_type = 0;
  int sci_obj = 0, lhs_obj = 0, rhs_obj = 0;
  double * ws  = NULL;
  int    * lws = NULL;
  double * v   = NULL; int nv = 0;
  double * x_filtersd = NULL, * x_lower_filtersd = NULL, * x_upper_filtersd = NULL;
  int maxa = 0, maxla = 0, maxu = 0, maxiu = 0, nout = 0, n = 0, m = 0, i = 0, j = 0, k = 0;
  double * tmp_ptr_dbl = NULL;
  fobj_filtersd_ptr func_obj  = NULL;
  grad_filtersd_ptr func_grad = NULL;
  const char * LabelList_out[10] = {"dnorm", "h", "hJt", "hJ", "ipeq", "k", "itn", "nft", "ngt", "cstype"};
  int nbvars_old = Nbvars, nbconstr = 0;
  int is_external_func = 0, is_scilab_func = 0;
  SciErr _SciErr;
  
  nout = 0; // stderr
  //nout = 6; // stdout
  
  if (Rhs<LAST_PARAM_IN)
    {
      Scierror(999,"%s: %d parameters are required\n", fname, LAST_PARAM_IN);
      return 0;
    } /* End If */

  ////////////////////////
  // Initialize commons //
  ////////////////////////

  C2F(defaultc).ainfty = 1.0e+020;
  C2F(defaultc).ubd    = 10000.0;
  C2F(defaultc).mlp    = 50;
  C2F(defaultc).mxf    = 50;

  C2F(wsc).kk          = 0;
  C2F(wsc).ll          = 0;
  C2F(wsc).kkk         = 0;
  C2F(wsc).lll         = 0;
  C2F(wsc).mxws        = 0;
  C2F(wsc).mxlws       = 0;

  C2F(statsc).dnorm    = 0.0;
  C2F(statsc).h        = 0.0;
  C2F(statsc).hJt      = 0.0;
  C2F(statsc).hJ       = 0.0;
  C2F(statsc).ipeq     = 0;
  C2F(statsc).k        = 0;
  C2F(statsc).itn      = 0;
  C2F(statsc).nft      = 0;
  C2F(statsc).ngt      = 0;

  C2F(ngrc).mxgr       = 1000000;

  C2F(mxm1c).mxm1      = 0;

  C2F(epsc).eps        = 1.11e-16;
  C2F(epsc).tol        = 9.9e-13;
  C2F(epsc).emin       = 0.0;

  C2F(repc).sgnf       = 1.0e-8;
  C2F(repc).nrep       = 2;
  C2F(repc).npiv       = 3;
  C2F(repc).nres       = 2;

  ////////////////////////
  // Get the parameters //
  ////////////////////////

  _SciErr = getVarAddressFromPosition(pvApiCtx, X_IN, &pi_x_in); FILTERSD_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, pi_x_in, &n_x_in, &m_x_in, &x_in); FILTERSD_ERROR;
  n = n_x_in*m_x_in;
  
  _SciErr = getVarAddressFromPosition(pvApiCtx, X_LOWER_IN, &pi_x_lower_in); FILTERSD_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, pi_x_lower_in, &n_x_lower_in, &m_x_lower_in, &x_lower_in); FILTERSD_ERROR;

  x_lower_filtersd = (double*)MALLOC(m_x_lower_in*n_x_lower_in*sizeof(double));
  for(i=0; i<m_x_lower_in*n_x_lower_in; i++) x_lower_filtersd[i] = x_lower_in[i];

  _SciErr = getVarAddressFromPosition(pvApiCtx, X_UPPER_IN, &pi_x_upper_in); FILTERSD_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, pi_x_upper_in, &n_x_upper_in, &m_x_upper_in, &x_upper_in); FILTERSD_ERROR;

  x_upper_filtersd = (double*)MALLOC(m_x_upper_in*n_x_upper_in*sizeof(double));
  for(i=0; i<m_x_upper_in*n_x_upper_in; i++) x_upper_filtersd[i] = x_upper_in[i];

#ifdef DEBUG
  for(i=0;i<m_x_upper_in*n_x_upper_in; i++) printf("DEBUG: %d - x_lower = %f, x_upper = %f\n", i, x_lower_filtersd[i], x_upper_filtersd[i]);
#endif

  // Get Fobj
  // Here, we get either a Scilab script function or a string which gives the name of the C/C++/Fortran function
  // loaded via "link".
  func_obj = (fobj_filtersd_ptr)GetFunctionPtr("filtersd_fobj", FOBJ_IN, FTab_filtersd_function, (voidf)C2F(dense_functions), &sci_obj, &lhs_obj, &rhs_obj);

  if ((lhs_obj!=2)&&(rhs_obj!=1))
    {
      sciprint("%s : Error - function must have 2 outputs parameters (f and c) and 1 input parameter (x)", fname);
      FREE(x_lower_filtersd);
      FREE(x_upper_filtersd);
      return 0;
    }

  param_fobj.fobj_type = 0;
  param_fobj.n_x       = n;
  param_fobj.n_c       = -1;
  param_fobj.is_sparse = 0;
  param_fobj.sci_obj   = sci_obj;
  param_fobj.lhs_obj   = lhs_obj;
  param_fobj.rhs_obj   = rhs_obj;
  param_fobj.stack_pos = Rhs + 1; // We preallocate the objective function value on position 3 on the stack
                                  // So, scifunction will create all the needed variables on LAST_PARAM_OUT + 1 to
                                  // avoid the destruction of this preallocated output variable

  if (func_obj==(fobj_filtersd_ptr)0) 
    {
      sciprint("%s : Error - Last argument must be a pointer to a scilab function", fname);
      FREE(x_lower_filtersd);
      FREE(x_upper_filtersd);
      return 0;
    }

  if (func_obj!=(voidf)C2F(dense_functions))
    {
      param_fobj.fobj_type = 1;
      param_fobj.function  = func_obj;
      is_external_func++;
    }
  else
    {
      param_fobj.fobj_type = 0;
      param_fobj.function  = NULL;
      is_scilab_func++;
    }
    
  // Get gradient
  // Here, we get either a Scilab script function or a string which gives the name of the C/C++/Fortran function
  // loaded via "link".
  func_grad = (grad_filtersd_ptr)GetFunctionPtr("filtersd_grad", GRAD_IN, FTab_filtersd_function, (voidf)C2F(dense_gradients), &sci_obj, &lhs_obj, &rhs_obj);
  
  if ((lhs_obj!=1)&&(rhs_obj!=1))
    {
      sciprint("%s : Error - gradient must have 1 output parameters (a) and 1 input parameters (x)", fname);
      FREE(x_lower_filtersd);
      FREE(x_upper_filtersd);
      return 0;
    }

  param_grad.fobj_type = 0;
  param_grad.n_x       = n;
  param_grad.n_c       = -1;
  param_grad.is_sparse = 0;
  param_grad.sci_obj   = sci_obj;
  param_grad.lhs_obj   = lhs_obj;
  param_grad.rhs_obj   = rhs_obj;
  param_grad.stack_pos = Rhs + 1; // We preallocate the objective function value on position 3 on the stack
                                  // So, scifunction will create all the needed variables on LAST_PARAM_OUT + 1 to
                                  // avoid the destruction of this preallocated output variable

  if (func_grad==(grad_filtersd_ptr)0) 
    {
      sciprint("%s : Error - Last argument must be a pointer to a scilab function", fname);
      FREE(x_lower_filtersd);
      FREE(x_upper_filtersd);
      return 0;
    }

  if (func_grad!=(voidf)C2F(dense_gradients))
    {
      param_grad.fobj_type = 1;
      param_grad.gradient  = func_grad;
      is_external_func++;
    }
  else
    {
      param_grad.fobj_type = 0;
      param_grad.gradient  = NULL;
      is_scilab_func++;
    }

  if ((is_external_func!=2) && (is_scilab_func!=2))
    {
      sciprint("%s : Error - if one function is a pointer to an external function the other must be a pointer to an external function.", fname);
      FREE(x_lower_filtersd);
      FREE(x_upper_filtersd);
      return 0;
    }

  initPList(pvApiCtx, PARAMS_IN, &param_in_addr);
  if (!checkPList(pvApiCtx, param_in_addr))
    {
      Scierror(999, "%s: argument n° %d is not a plist\n", fname, PARAMS_IN);
      FREE(x_lower_filtersd);
      FREE(x_upper_filtersd);
      return 0;
    }

  getDoubleInPList(pvApiCtx, param_in_addr, "rho", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) rho = tmp_double;

  getDoubleInPList(pvApiCtx, param_in_addr, "htol", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) htol = tmp_double;

  getDoubleInPList(pvApiCtx, param_in_addr, "rgtol", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) rgtol = tmp_double;

  getDoubleInPList(pvApiCtx, param_in_addr, "ainfty", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) ainfty = tmp_double;

  getDoubleInPList(pvApiCtx, param_in_addr, "ubd", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) ubd = tmp_double;

  getDoubleInPList(pvApiCtx, param_in_addr, "fmin", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) fmin = tmp_double;

  getIntInPList(pvApiCtx, param_in_addr, "maxit", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) maxit = tmp_int;

  getIntInPList(pvApiCtx, param_in_addr, "iprint", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) iprint = tmp_int;

  kmax = n;
  getIntInPList(pvApiCtx, param_in_addr, "kmax", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) kmax = tmp_int;
  kmax = MIN(n, kmax);

  getIntInPList(pvApiCtx, param_in_addr, "maxg", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) maxg = tmp_int;

  getIntInPList(pvApiCtx, param_in_addr, "mlp", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) mlp = tmp_int;

  getIntInPList(pvApiCtx, param_in_addr, "mxf", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) mxf = tmp_int;

  getIntInPList(pvApiCtx, param_in_addr, "mxgr", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) mxgr = tmp_int;

  getIntInPList(pvApiCtx, param_in_addr, "mxws", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) mxws = tmp_int;

  getIntInPList(pvApiCtx, param_in_addr, "mxlws", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) mxlws = tmp_int;

#ifdef DEBUG
  printf("DEBUG:defaultc.ainfty = %f\n", C2F(defaultc).ainfty);
  printf("DEBUG:defaultc.ubd    = %f\n", C2F(defaultc).ubd);
  printf("DEBUG:defaultc.mlp    = %d\n", C2F(defaultc).mlp);
  printf("DEBUG:defaultc.mxf    = %d\n", C2F(defaultc).mxf);

  printf("DEBUG:wsc.kk    = %d\n", C2F(wsc).kk);
  printf("DEBUG:wsc.ll    = %d\n", C2F(wsc).ll);
  printf("DEBUG:wsc.kkk   = %d\n", C2F(wsc).kkk);
  printf("DEBUG:wsc.lll   = %d\n", C2F(wsc).lll);
  printf("DEBUG:wsc.mxws  = %d\n", C2F(wsc).mxws);
  printf("DEBUG:wsc.mxlws = %d\n", C2F(wsc).mxlws);

  printf("DEBUG:statsc.dnorm = %f\n", C2F(statsc).dnorm);
  printf("DEBUG:statsc.h     = %f\n", C2F(statsc).h);
  printf("DEBUG:statsc.hJt   = %f\n", C2F(statsc).hJt);
  printf("DEBUG:statsc.hJ    = %f\n", C2F(statsc).hJ);
  printf("DEBUG:statsc.ipeq  = %d\n", C2F(statsc).ipeq);
  printf("DEBUG:statsc.k     = %d\n", C2F(statsc).k);
  printf("DEBUG:statsc.itn   = %d\n", C2F(statsc).itn);
  printf("DEBUG:statsc.nft   = %d\n", C2F(statsc).nft);
  printf("DEBUG:statsc.ngt   = %d\n", C2F(statsc).ngt);

  printf("DEBUG:ngrc.mxgr = %d\n", C2F(ngrc).mxgr);

  printf("DEBUG:mxm1c.mxml = %d\n", C2F(mxm1c).mxm1);
#endif

  /**********************
   * Pre-Initialisation *
   **********************/

  // maxu    length of workspace user(*) passed through to user subroutines 'functions' and 'gradients'
  // maxiu   length of workspace iuser(*) passed through to user subroutines
  // maxa    maximum number of entries in the Jacobian a(*) set by gradients
  // maxla   number of entries required for sparse matrix indices and pointers
  //         la(0:*) to be set up in lws(*) (maxla>=maxa+m+3).
  //         Set maxla=1 if using dense matrix format 

  maxu  = 0;
  maxiu = 0;

  maxla = 1;

  // set stride in lws(maxiu+1) and constant elements of a(*) in ws(maxu+1) on
  
  // La taille de ws et lws sont reglees dans le common wsc via les variables mxws mxlws

  C2F(wsc).mxws = mxws;
  ws = (double *)MALLOC(C2F(wsc).mxws*sizeof(double));
  for(i=0;i<C2F(wsc).mxws;i++) ws[i] = 0.0;

  C2F(wsc).mxlws = mxlws;
  lws = (int *)MALLOC(C2F(wsc).mxlws*sizeof(int));
  for(i=0;i<C2F(wsc).mxlws;i++) lws[i] = 0;
  lws[maxiu+0] = n; // stride initialization

  /*********************
   * Call the function *
   *********************/

  if (!is_external_func)
    {
      Nbvars = param_fobj.stack_pos + MAX(param_fobj.rhs_obj,param_fobj.lhs_obj) + 3;
      
      _SciErr = createMatrixOfDouble(pvApiCtx, param_fobj.stack_pos+0, n, 1, x_in); FILTERSD_ERROR;

      _SciErr = createMatrixOfDouble(pvApiCtx, param_fobj.stack_pos+1, 0, 0, NULL); FILTERSD_ERROR;
      
      if (!C2F(scifunction)(&param_fobj.stack_pos, &param_fobj.sci_obj, &param_fobj.lhs_obj, &param_fobj.rhs_obj))
	{
	  Scierror(999,"filtersd: error while calling the functions\n");
	  Nbvars = nbvars_old;
	  FREE(x_lower_filtersd);
	  FREE(x_upper_filtersd);
	  FREE(ws);
	  FREE(lws);
	  return 0;
	}

      if (Err>0) 
	{
	  Scierror(999,"filtersd: error when calling objective function\n");
	  Nbvars = nbvars_old;
	  FREE(x_lower_filtersd);
	  FREE(x_upper_filtersd);
	  FREE(ws);
	  FREE(lws);
	  return 0;
	} 
      
      // Get F
      _SciErr = getVarAddressFromPosition(pvApiCtx, param_fobj.stack_pos+0, &pi_tmp_var); FILTERSD_ERROR;
      
      // Get C
      _SciErr = getVarAddressFromPosition(pvApiCtx, param_fobj.stack_pos+1, &pi_tmp_var); FILTERSD_ERROR;
      _SciErr = getMatrixOfDouble(pvApiCtx, pi_tmp_var, &m_tmp_var, &n_tmp_var, &tmp_var);

      nbconstr = m_tmp_var * n_tmp_var;
      m = nbconstr;

      param_grad.n_c = nbconstr;
      param_fobj.n_c = nbconstr;

      /*********************
       * Call the gradient *
       *********************/

      _SciErr = createMatrixOfDouble(pvApiCtx, param_grad.stack_pos+0, n, 1, x_in); FILTERSD_ERROR;
      
      if (!C2F(scifunction)(&param_grad.stack_pos, &param_grad.sci_obj, &param_grad.lhs_obj, &param_grad.rhs_obj))
	{
	  Scierror(999,"filtersd: error while calling the gradient function\n");
	  Nbvars = nbvars_old;
	  FREE(x_lower_filtersd);
	  FREE(x_upper_filtersd);
	  FREE(ws);
	  FREE(lws);
	  return 0;
	}
      
      if (Err>0) 
	{
	  Scierror(999,"filtersd: error when calling gradient function\n");
	  Nbvars = nbvars_old;
	  FREE(x_lower_filtersd);
	  FREE(x_upper_filtersd);
	  FREE(ws);
	  FREE(lws);
	  return 0;
	} 
      
      _SciErr = getVarAddressFromPosition(pvApiCtx, param_grad.stack_pos+0, &pi_tmp_var); FILTERSD_ERROR;
      _SciErr = getVarType(pvApiCtx, pi_tmp_var, &grad_type);
      
      if (grad_type==sci_sparse) param_grad.is_sparse = 1;

      if (!param_grad.is_sparse)
	{
#ifdef DEBUG
	  printf("DEBUG: gradient matrix is full\n");
#endif

	  _SciErr = getMatrixOfDouble(pvApiCtx, pi_tmp_var, &m_tmp_var, &n_tmp_var, &tmp_var); FILTERSD_ERROR;

	  if (m_tmp_var*n_tmp_var != n*(m+1))
	    {
	      Scierror(999,"%s: gradients - a must be of size %d * %d\n", fname, n, m+1);
	      Scierror(999,"%s: current size: %d * %d\n", fname, n_tmp_var, m_tmp_var);
	      Nbvars = nbvars_old;
	      FREE(x_lower_filtersd);
	      FREE(x_upper_filtersd);
	      FREE(ws);
	      FREE(lws);
	      return 0;
	    }
	  
	  k = 0;
	  for(i=0;i<m_tmp_var;i++) 
	    {
	      for(j=0; j<n_tmp_var;j++)
		{
		  ws[maxu+k] = tmp_var[j + i*n_tmp_var]; // store a
		  ws[maxa+maxu+k] = ws[maxu+k]; // replicate a
#ifdef DEBUG
		  printf("ws[%d] = %f  - replicate ws[%d] = %f \n", maxu+k, ws[maxu+k], maxa+maxu+k, ws[maxa+maxu+k]);
#endif
		  k++;
		}
	    }
	}
      else
	{
#ifdef DEBUG
	  printf("DEBUG: gradient matrix is sparse\n");
#endif
	  Scierror(999,"%s: does not manage sparse matrix\n", fname);
	  return 0;
	}

      Nbvars = nbvars_old;
    }
  else
    {
      getIntInPList(pvApiCtx, param_in_addr, "nbconstr", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) nbconstr = tmp_int;

      if (nbconstr==-1)
	{
	  Scierror(999,"%s: call to external C functions - nbconstr must be set in params_in\n", fname);
	  FREE(x_lower_filtersd);
	  FREE(x_upper_filtersd);
	  FREE(ws);
	  FREE(lws);
	  return 0;
	}
    }

  m = nbconstr;
  cstype = (char *)MALLOC((m+1)*sizeof(char));
  for(i=0;i<m;i++) cstype[i] = '0';

  // Now, we can initialize maxa
  maxa = n*(m+1);

  if (m_x_upper_in*n_x_upper_in!=n+m)
    {
      sciprint("%s : Error - x_upper must be of dimension %d", fname, n+m);
      FREE(x_lower_filtersd);
      FREE(x_upper_filtersd);
      FREE(ws);
      FREE(lws);
      FREE(cstype);
      return 0;
    }
  
  if (m_x_lower_in*n_x_lower_in!=n+m)
    {
      sciprint("%s : Error - x_lower must be of dimension %d", fname, n+m);
      FREE(x_lower_filtersd);
      FREE(x_upper_filtersd);
      FREE(ws);
      FREE(lws);
      FREE(cstype);
      return 0;
    }

  lambda_out = (double *)MALLOC((n+m+1)*sizeof(double));
  for(i=0;i<n;i++)       lambda_out[i] = 1.0e-2;
  for(i=n;i<(n+m+1);i++) lambda_out[i] = 0.0;

  ////////////////////
  // Initialization //
  ////////////////////

  // Allouer l'espace pour
  // al  -> n+m
  // bl  -> n+m
  // bu  -> n+m
  // ws  -> ws(mxws)
  // lws -> ws(mxlws)
  // v   -> maxg (then nv = maxg - otherwise, set nv=1 and v(1) = 1)

  v = (double *)MALLOC((maxg+1)*sizeof(double));
  maxg = MIN(n,maxg)+1;

  for(i=0;i<maxg; i++) v[i] = 0.0;
  nv = 0;

  x_filtersd = (double *)MALLOC((n+m+1)*sizeof(double));
  for(i=n;i<(n+m+1);i++) x_filtersd[i] = 0.0; // initialize de x workspace
  for(i=0;i<n;i++)       x_filtersd[i] = x_in[i];

#ifdef DEBUG
  for(i=0;i<(n+m); i++) printf("DEBUG: x_upper[%d] = %f - x_in[%d] = %f - x_lower[%d] = %f\n", i, x_upper_in[i], i, x_in[i], i, x_lower_in[i]);
#endif
    
  // common/defaultc/ainfty,ubd,mlp,mxf -> 1e20, 1e4, 50, 50 (default values)
  // ainfty: 1e20 (default value) represent infinity
  // ubd: 1e4 (default values) upper bound on the allowed constraint violation
  // mlp: 50 (default values) maximum length of arrays used in the degeneracy control
  // mxf: 50 (default values) maximum length of filter arrays

  C2F(defaultc).ainfty = ainfty;
  C2F(defaultc).ubd    = ubd;
  C2F(defaultc).mlp    = mlp;
  C2F(defaultc).mxf    = mxf;

  // common/ngrc/mxgr
  // mxgr: 1000000 (default value) maximum time spent in each call of the LCP solver

  C2F(ngrc).mxgr = mxgr;

  // common/mxm1c/mxm1
  // mxm1: min(m+1,n) (default value) maximum number of general constraints in the active set

  mxm1 = MIN(m+1, n);
  getIntInPList(pvApiCtx, param_in_addr, "mxm1", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) mxm1 = tmp_int;

  C2F(mxm1c).mxm1 = mxm1;

  // common/statsc/dnorm,h,hJt,hJ,ipeq,k,itn,nft,ngt

  // dnorm: final step length
  // h:     final constraint violation
  // hJt:   ditto for 'N' constraints
  // hJ:    ditto for 'A' and 'Z' constraints
  // ipeq:  number of active equations
  // k:     number of free variables
  // itn:   number of iterations
  // nft:   total number of function calls
  // ngt:   total number of gradient calls

#ifdef DEBUG
  sciprint("DEBUG: defaultc common - ainfty = %f\n", C2F(defaultc).ainfty);
  sciprint("DEBUG: defaultc common - ubd    = %f\n", C2F(defaultc).ubd);
  sciprint("DEBUG: defaultc common - mlp    = %d\n", C2F(defaultc).mlp);
  sciprint("DEBUG: defaultc common - mxf    = %d\n", C2F(defaultc).mxf);

  sciprint("DEBUG: statsc common - dnorm = %f\n", C2F(statsc).dnorm);
  sciprint("DEBUG: statsc common - h     = %f\n", C2F(statsc).h);
  sciprint("DEBUG: statsc common - hJt   = %f\n", C2F(statsc).hJt);
  sciprint("DEBUG: statsc common - hJ    = %f\n", C2F(statsc).hJ);
  sciprint("DEBUG: statsc common - ipeq  = %d\n", C2F(statsc).ipeq);
  sciprint("DEBUG: statsc common - k     = %d\n", C2F(statsc).k);
  sciprint("DEBUG: statsc common - itn   = %d\n", C2F(statsc).itn);
  sciprint("DEBUG: statsc common - nft   = %d\n", C2F(statsc).nft);
  sciprint("DEBUG: statsc common - ngt   = %d\n", C2F(statsc).ngt);
#endif

  // Parameters for glcpd:
  // common/epsc/eps,tol,emin
  // common/repc/sgnf,nrep,npiv,nres
  // common/refactorc/mc,mxmc

  // f=x(1)+x(2)+x(3)
  // c(1)=(x(4)+x(6))-1.0
  // c(2)=(x(5)+x(7)-x(4))-1.0
  // c(3)=1e-2*(x(8)-x(5))-1.0
  // c(4)=1e2*x(1)+833.33252*x(4)-x(1)*x(6)-83333.333
  // c(5)=125e1*(x(5)-x(4))+x(2)*(x(4)-x(7))
  // c(6)=125e4-25e2*x(5)+x(3)*(x(5)-x(8))

  // lws(1)=n -> n = 8, m = 6 (stride is equal to 8) with the sparse format, this will be not constant and 
  // not always equal to 8
  // The constant part for the gradient

  // a(1,0) = 1 \
  // a(2,0) = 1  |
  // a(3,0) = 1  |
  // a(4,0) = 0   > Objective function
  // a(5,0) = 0  |
  // a(6,0) = 0  |
  // a(7,0) = 0  |
  // a(8,0) = 0 /

  // a(1,1) = 0 \
  // a(2,1) = 0  |
  // a(3,1) = 0  |
  // a(4,1) = 1   > Constraint c(1)
  // a(5,1) = 0  |
  // a(6,1) = 1  |
  // a(7,1) = 0  |
  // a(8,1) = 0 /

  // a(1,2) = 0 \
  // a(2,2) = 0  |
  // a(3,2) = 0  |
  // a(4,2) = -1  > Constraint c(2)
  // a(5,2) = 1  |
  // a(6,2) = 0  |
  // a(7,2) = 1  |
  // a(8,2) = 0 /

  // a(1,3) = 0     \
  // a(2,3) = 0      |
  // a(3,3) = 0      |
  // a(4,3) = 0       > Constraint c(3)
  // a(5,3) = -1e-2  |
  // a(6,3) = 0      |
  // a(7,3) = 0      |
  // a(8,3) = 1e-2  /

  // a(1,4) = NC        \
  // a(2,4) = 0          |
  // a(3,4) = 0          |
  // a(4,4) = 833.33252   > Constraint c(4)
  // a(5,4) = 0          |
  // a(6,4) = NC         |
  // a(7,4) = 0          |
  // a(8,4) = 0         /

  // a(1,5) = 0    \
  // a(2,5) = NC    |
  // a(3,5) = 0     |
  // a(4,5) = NC     > Constraint c(5)
  // a(5,5) = 125e1 |
  // a(6,5) = 0     |
  // a(7,5) = NC    |
  // a(8,5) = 0    /

  // a(1,6) = 0  \
  // a(2,6) = 0   |
  // a(3,6) = NC  |
  // a(4,6) = 0    > Constraint c(6)
  // a(5,6) = NC  |
  // a(6,6) = 0   |
  // a(7,6) = 0   |
  // a(8,6) = NC /

  // The other part:
  // a(1,4)=1e2-x(6)   -> constraint 4, variable 1
  // a(6,4)=-x(1)      -> constraint 4, variable 6
  // a(2,5)=x(4)-x(7)  -> constraint 5, variable 2
  // a(4,5)=x(2)-125e1 -> constraint 5, variable 4
  // a(7,5)=-x(2)      -> constraint 5, variable 7
  // a(3,6)=x(5)-x(8)  -> constraint 6, variable 3
  // a(5,6)=x(3)-25e2  -> constraint 6, variable 5
  // a(8,6)=-x(3)      -> constraint 6, variable 8
  
  // For sparse format
  // data lws(0:25) /26,1,2,3,4,5,6,7,8,4,6,4,5,7,5,8,1,4,6,2,4,5,7,3,5,8/
  // data lws(26:33)/1,9,11,14,16,19,23,26/
  // data ws(2:26)/3*1.D0,5*0.D0,t,t,tm,t,t,-1.D-2,1.D-2,0.D0,833.33252D0,3*0.D0,125.D1,4*0.D0/

  // lws(0): the length of lws
  // lws(1:8): indexes of coeff of the objective function (non sparse)
  // lws(9:10): indexes of coeff of constraint 1
  // lws(11:13): indexes of coeff of constraint 2
  // lws(14:15): indexes of coeff of constraint 3
  // lws(16:18): indexes of coeff of constraint 4
  // lws(19:22): indexes of coeff of constraint 5
  // lws(23:25): indexes of coeff of constraint 6

  // lws(ws(ws(0)+0):ws(ws(0)+1)-1): indexes of coeff of objective function
  // lws(ws(ws(0)+1):ws(ws(0)+2)-1): indexes of coeff of constraint 1

  //                           n,    m,    x,       al,      f,       fmin,    cstype,bl,      bu,      ws,      lws,  v,       nv,
  // extern void C2F(filtersd)(int *,int *,double *,double *,double *,double *,char *,double *,double *,double *,int *,double *,int *,
  //                           int *,int *,int *,int *,int *,int *,double *,double *,double *,int *,int *, int *,int *);
  //                           maxa, maxla,maxu, maxiu,kmax, maxg, rho,     htol,    rgtol,   maxit,iprint,nout, ifail)

  /*
   * Calling filterSD
   */

  if (!is_external_func)
    {
      C2F(ds_filtersd)(&n,&m,x_filtersd,lambda_out,&f_tmp_out,&fmin,cstype,x_lower_filtersd,x_upper_filtersd,
		       ws,lws,v,&nv,&maxa,&maxla,&maxu,&maxiu,&kmax,&maxg,&rho,&htol,&rgtol,&maxit,&iprint,&nout,&ifail,
		       C2F(dense_functions),C2F(dense_gradients));
    }
  else
    {
      C2F(ds_filtersd)(&n,&m,x_filtersd,lambda_out,&f_tmp_out,&fmin,cstype,x_lower_filtersd,x_upper_filtersd,
		       ws,lws,v,&nv,&maxa,&maxla,&maxu,&maxiu,&kmax,&maxg,&rho,&htol,&rgtol,&maxit,&iprint,&nout,&ifail,
		       (fobj_filtersd_ptr)param_fobj.function,(grad_filtersd_ptr)param_grad.function);
    }

  //////////////////////////
  // return the variables //
  //////////////////////////

  n_x_out = n_x_in; m_x_out = m_x_in;
  _SciErr = allocMatrixOfDouble(pvApiCtx, X_OUT, n_x_out, m_x_out, &x_out); FILTERSD_ERROR;
  for(i=0; i<n_x_out*m_x_out; i++) x_out[i] = x_filtersd[i];

  // f_out
  n_f_out = 1; m_f_out = 1;
  _SciErr = allocMatrixOfDouble(pvApiCtx, F_OUT, n_f_out, m_f_out, &f_out); FILTERSD_ERROR;
  f_out[0] = f_tmp_out;

  // ifail_out

  //  ifail   returns failure indication as follows
  //               0 = successful run
  //               1 = unbounded NLP (f <= fmin at an htol-feasible point)
  //               2 = bounds on x are inconsistent
  //               3 = local minimum of feasibility problem and h > htol
  //                   (nonlinear constraints are locally inconsistent)
  //               4 = initial point x has h > ubd (reset ubd or x and re-enter)
  //               5 = maxit major iterations have been carried out
  //               6 = termination with rho <= htol
  //               7 = not enough workspace in ws or lws (see message)
  //               8 = insufficient space for filter (increase mxf and re-enter)

  //              >9 = unexpected fail in LCP solver (10 has been added to ifail)

  //              10 = solution obtained
  //              11 = unbounded problem (f(x)<fmin has occurred: note grad is not
  //                    evaluated in this case)
  //              12 = bl(i) > bu(i) for some i
  //              13 = infeasible problem detected in Phase 1
  //              14 = line search cannot improve f (possibly increase rgtol)
  //              15 = mxgr gradient calls exceeded (this test is only carried
  //                    out at the start of each iteration)
  //              16 = incorrect setting of m, n, kmax, maxg, mlp, m0de or tol
  //              17 = not enough space in ws or lws
  //              18 = not enough space in lp (increase mlp)
  //              19 = dimension of reduced space too large (increase kmax)
  //              20 = maximum number of unsuccessful restarts taken
  //             >20= possible use by later sparse matrix codes

  n_ifail_out = 1; m_ifail_out = 1;
  _SciErr = allocMatrixOfDouble(pvApiCtx, IFAIL_OUT, n_ifail_out, m_ifail_out, &tmp_ptr_dbl); FILTERSD_ERROR;
  tmp_ptr_dbl[0] = (double)ifail;

  // lambda_out
  n_lambda_out = n+m; m_lambda_out = 1;
  _SciErr = allocMatrixOfDouble(pvApiCtx, LAMBDA_OUT, n_lambda_out, m_lambda_out, &tmp_ptr_dbl); FILTERSD_ERROR;
  for(i=0; i<n_lambda_out*m_lambda_out; i++) tmp_ptr_dbl[i] = lambda_out[i];

#ifdef DEBUG
  sciprint("DEBUG: statsc common - dnorm = %f\n", C2F(statsc).dnorm);
  sciprint("DEBUG: statsc common - h     = %f\n", C2F(statsc).h);
  sciprint("DEBUG: statsc common - hJt   = %f\n", C2F(statsc).hJt);
  sciprint("DEBUG: statsc common - hJ    = %f\n", C2F(statsc).hJ);
  sciprint("DEBUG: statsc common - ipeq  = %d\n", C2F(statsc).ipeq);
  sciprint("DEBUG: statsc common - k     = %d\n", C2F(statsc).k);
  sciprint("DEBUG: statsc common - itn   = %d\n", C2F(statsc).itn);
  sciprint("DEBUG: statsc common - nft   = %d\n", C2F(statsc).nft);
  sciprint("DEBUG: statsc common - ngt   = %d\n", C2F(statsc).ngt);
  sciprint("DEBUG: cstype                = %s\n", cstype);
  sciprint("DEBUG: X_OUT = %d\n", X_OUT);
  sciprint("DEBUG: F_OUT = %d\n", F_OUT);
  sciprint("DEBUG: IFAIL_OUT = %d\n", IFAIL_OUT);
  sciprint("DEBUG: LAMBDA_OUT = %d\n", LAMBDA_OUT);
  sciprint("DEBUG: PARAMS_OUT = %d\n", PARAMS_OUT);
  sciprint("DEBUG: LAST_PARAM_OUT = %d\n", LAST_PARAM_OUT);
#endif

  _SciErr = createPList(pvApiCtx, PARAMS_OUT, &param_out_addr, (char **)LabelList_out, 10); FILTERSD_ERROR;

  _SciErr = createDoubleInPList(pvApiCtx, PARAMS_OUT, param_out_addr, "dnorm", C2F(statsc).dnorm); FILTERSD_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, PARAMS_OUT, param_out_addr, "h",     C2F(statsc).h);     FILTERSD_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, PARAMS_OUT, param_out_addr, "hJt",   C2F(statsc).hJt);   FILTERSD_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, PARAMS_OUT, param_out_addr, "hJ",    C2F(statsc).hJ);    FILTERSD_ERROR;
  
  _SciErr = createIntInPList(pvApiCtx, PARAMS_OUT, param_out_addr, "ipeq", C2F(statsc).ipeq); FILTERSD_ERROR;
  _SciErr = createIntInPList(pvApiCtx, PARAMS_OUT, param_out_addr, "k",    C2F(statsc).k);    FILTERSD_ERROR;
  _SciErr = createIntInPList(pvApiCtx, PARAMS_OUT, param_out_addr, "itn",  C2F(statsc).itn);  FILTERSD_ERROR;
  _SciErr = createIntInPList(pvApiCtx, PARAMS_OUT, param_out_addr, "nft",  C2F(statsc).nft);  FILTERSD_ERROR;
  _SciErr = createIntInPList(pvApiCtx, PARAMS_OUT, param_out_addr, "ngt",  C2F(statsc).ngt);  FILTERSD_ERROR;

  _SciErr = createStringInPList(pvApiCtx, PARAMS_OUT, param_out_addr, "cstype", cstype);

  LhsVar(1) = X_OUT;
  LhsVar(2) = F_OUT;
  LhsVar(3) = IFAIL_OUT;
  LhsVar(4) = LAMBDA_OUT;
  LhsVar(5) = PARAMS_OUT;

  // YC: Plantage ici sinon .... 

  if (x_lower_filtersd) FREE(x_lower_filtersd);
  if (x_upper_filtersd) FREE(x_upper_filtersd);
  if (lambda_out)       FREE(lambda_out);
  if (ws)               FREE(ws);
  if (lws)              FREE(lws);
  if (cstype)           FREE(cstype);
  if (v)                FREE(v);
  if (x_filtersd)       FREE(x_filtersd);

  return 0;
}

/***************************
 * Functions and Gradients *
 ***************************/

// This function is a wrapper which allows to call a Scilab script function like a "normal" C/C++/Fortran function
void C2F(dense_functions)(int * n, int * m, double * x, double * f, double * c, double * user, int * iuser)
{
  int n_x = 0;
  int n_c = 0;
  int * tmp_addr = NULL;
  int sci_obj, lhs_obj, rhs_obj;
  int stack_pos, i;
  int nbvars_old = Nbvars;
  double * tmp_var = NULL, tmp_val = 0.0;
  double f_out = 0.0;
  int n_c_out = 0, m_c_out = 0;
  SciErr _SciErr;

#ifdef DEBUG
  printf("DEBUG: Functions called\n");
#endif

  if (param_fobj.fobj_type==0)
    {
      sci_obj   = param_fobj.sci_obj;
      lhs_obj   = param_fobj.lhs_obj;
      rhs_obj   = param_fobj.rhs_obj;
      stack_pos = param_fobj.stack_pos;
      n_x       = param_fobj.n_x;
      if (param_fobj.n_c!=-1) n_c = param_fobj.n_c;
      else                    n_c = *m;

      Nbvars = stack_pos + MAX(rhs_obj,lhs_obj) + 3;
      
#ifdef DEBUG
      printf("DEBUG: n = %d, m = %d, n_x = %d, n_c = %d\n", *n, *m, n_x, n_c);
      for(i=0; i<n_x; i++) printf("DEBUG: x[%d] = %f\n", i, x[i]);
#endif
      _SciErr = createMatrixOfDouble(pvApiCtx, stack_pos+0, n_x, 1, x); FILTERSD_ERROR_NORETURN;

      // The scilab function return 2 output argument: f and c
      _SciErr = createMatrixOfDouble(pvApiCtx, stack_pos+1, 0, 0, &tmp_val); FILTERSD_ERROR_NORETURN;
      
      _SciErr = createMatrixOfDouble(pvApiCtx, param_fobj.stack_pos+1, 0, 0, NULL); FILTERSD_ERROR_NORETURN;

      if (!C2F(scifunction)(&stack_pos, &sci_obj, &lhs_obj, &rhs_obj))
	{
	  Scierror(999,"filtersd: error when calling objective function\n");
	  Nbvars = nbvars_old;
	  return;
	}
      
      if (Err>0) 
	{
	  Scierror(999,"filtersd: error when calling objective function\n");
	  Nbvars = nbvars_old;
	  return;
	} 
      
      // Get F
      _SciErr = getVarAddressFromPosition(pvApiCtx, stack_pos+0, &tmp_addr); FILTERSD_ERROR_NORETURN;
      getScalarDouble(pvApiCtx, tmp_addr, &f_out);
      *f = f_out;
#ifdef DEBUG
      printf("DEBUG: f = %f\n", *f);
#endif
      
      // Get C
      _SciErr = getVarAddressFromPosition(pvApiCtx, stack_pos+1, &tmp_addr); FILTERSD_ERROR_NORETURN;
      _SciErr = getMatrixOfDouble(pvApiCtx, tmp_addr, &m_c_out, &n_c_out, &tmp_var);
      if (n_c_out*m_c_out!=*m)
	{
	  Scierror(999,"filtersd: objective function - wrong size for the constraints, must be of size %d.\n", *m);
	  Nbvars = nbvars_old;
	  return;
	}
      for(i=0; i<n_c; i++) 
	{
	  c[i] = tmp_var[i];
#ifdef DEBUG
	  printf("DEBUG: c[%d] = %f\n", i, c[i]);
#endif
	}

      Nbvars = nbvars_old;
    }
  else
    {
      (*param_fobj.function)(n, m, x, f, c, NULL, NULL);
      //Scierror(999,"filtersd_dense: call to C functions not implemented yet\n");
      //return;
    }

#ifdef DEBUG
  printf("DEBUG: End of Functions\n");
#endif
}

// This function is a wrapper which allows to call a Scilab script function like a "normal" C/C++/Fortran function
void C2F(dense_gradients)(int * n, int * m, double * x, double * a, double * user, int * iuser)
{
  int n_x = 0;
  int n_c = 0;
  int n_tmp, m_tmp, * tmp_addr = NULL;
  int sci_obj, lhs_obj, rhs_obj;
  int stack_pos, i = 0, j = 0, k = 0;
  int nbvars_old = Nbvars;
  double * tmp_var = NULL;
  SciErr _SciErr;

#ifdef DEBUG
  printf("DEBUG: Gradients called\n");
  printf("DEBUG: n = %d, m = %d\n", *n, *m);
#endif

  if (param_grad.fobj_type==0)
    {
      sci_obj   = param_grad.sci_obj;
      lhs_obj   = param_grad.lhs_obj;
      rhs_obj   = param_grad.rhs_obj;
      stack_pos = param_grad.stack_pos;
      n_x       = param_grad.n_x;
      if (param_grad.n_c!=-1) n_c = param_grad.n_c;
      else                    n_c = *m;
      
      Nbvars = stack_pos + MAX(rhs_obj,lhs_obj) + 3;

#ifdef DEBUG
      for(i=0; i<n_x; i++) printf("DEBUG: x[%d] = %f\n", i, x[i]);
#endif

      _SciErr = createMatrixOfDouble(pvApiCtx, stack_pos+0, n_x, 1, x); FILTERSD_ERROR_NORETURN;
      
      if (!C2F(scifunction)(&stack_pos, &sci_obj, &lhs_obj, &rhs_obj))
	{
	  Scierror(999,"filtersd: error when calling gradient function\n");
	  Nbvars = nbvars_old;
	  return;
	}
      
      if (Err>0) 
	{
	  Scierror(999,"filtersd: error when calling gradient function\n");
	  Nbvars = nbvars_old;
	  return;
	} 
      
      _SciErr = getVarAddressFromPosition(pvApiCtx, stack_pos+0, &tmp_addr); FILTERSD_ERROR_NORETURN;
      _SciErr = getMatrixOfDouble(pvApiCtx, tmp_addr, &m_tmp, &n_tmp, &tmp_var); FILTERSD_ERROR_NORETURN;
      if (m_tmp*n_tmp != (n_x*(n_c+1)))
	{
	  Scierror(999,"filtersd: gradients - a must be of size %d * %d\n", n_x, n_c+1);
	  Scierror(999,"          current size: %d * %d\n", n_tmp, m_tmp);
	  Nbvars = nbvars_old;
	  return;
	}

#ifdef DEBUG
      for(i=0; i<m_tmp*n_tmp; i++) printf("DEBUG: tmp_var[%d] = %f\n", i, tmp_var[i]);
#endif

      if (!param_grad.is_sparse)
	{
#ifdef DEBUG
	  printf("DEBUG: gradient matrix is full\n");
#endif
	  k = 0;
	  for(i=0;i<n_c;i++) 
	    {
	      for(j=0;j<n_x;j++)
		{
#ifdef DEBUG
		  printf("DEBUG a_before[%d] = %f", k, a[k]);
#endif
		  a[k] = tmp_var[j+(i*n_x)];
#ifdef DEBUG
		  printf(" a[%d] = %f\n", k, a[k]);
#endif
		  k++;
		}
	    }
	}
      else
	{
#ifdef DEBUG
	  printf("DEBUG: gradient matrix is sparse\n");
#endif
	  Scierror(999,"filtersd: does not manage sparse matrix\n");
	  return;
	}

      Nbvars = nbvars_old;
    }
  else
    {
      (*param_grad.gradient)(n, m, x, a, NULL, NULL);
      //Scierror(999,"filtersd_dense: call to C functions not implemented yet\n");
      //return;
    }

#ifdef DEBUG
  printf("DEBUG: End of Gradients\n");
#endif
}
