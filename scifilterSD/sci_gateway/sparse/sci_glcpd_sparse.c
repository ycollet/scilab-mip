#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <setjmp.h>

#include "sp_common.h"

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

#define GLCPD_ERROR if(_SciErr.iErr)		\
    {                                           \
      printError(&_SciErr, 0);                  \
      return _SciErr.iErr;                      \
    }

#define GLCPD_ERROR_NORETURN if(_SciErr.iErr)	\
    {							\
      printError(&_SciErr, 0);				\
      return;						\
    }

#define GLCPD_ERROR_RETURN_NULL if(_SciErr.iErr)	\
    {							\
      printError(&_SciErr, 0);				\
      return NULL;					\
    }

#define X_IN          1
#define FOBJ_IN       2
#define GRAD_IN       3
#define A_IN          4
#define LA_IN         5
#define X_LOWER_IN    6
#define X_UPPER_IN    7
#define PARAMS_IN     8
#define LAST_PARAM_IN 8

#define X_OUT          (Rhs+1)
#define F_OUT          (Rhs+2)
#define IFAIL_OUT      (Rhs+3)
#define PARAMS_OUT     (Rhs+4)
#define LAST_PARAM_OUT (Rhs+4)

void C2F(sparse_glcpd_functions)(int * n, double * x, double * f, double * ws, int * lws, char * cws);
void C2F(sparse_glcpd_gradients)(int * n, double * x, double * g, double * ws, int * lws, char * cws);

struct param_obj {
  int fobj_type; // 0 for a Scilab function, 1 for a C function
  int sci_obj;
  int lhs_obj;
  int rhs_obj;
  int stack_pos;
  int n_x;
  fobj_glcpd_ptr function;
  grad_glcpd_ptr gradient;
};

/*************************
 * Definition of FTables *
 *************************/

static FTAB FTab_glcpd_function[] = {{(char *) 0, (voidf) 0}};

/************************
 * The Scilab interface *
 ************************/

// static jmp_buf nlopt_function_env; 

static struct param_obj param_fobj;
static struct param_obj param_grad;

// Parameters in:
// - "rgtol"
// - "ainfty"
// - "ubd"
// - "fmin"
// - "iprint"
// - "kmax"
// - "maxg"
// - "mlp"
// - "mode"
// - "mxgr"
// - "mxws"
// - "mwlws"
// - "nout"

// Parameters out:
// - "ifail"
// - "rgnorm" - infoc
// - "vstep"  - infoc
// - "iter"   - infoc
// - "npv"    - infoc
// - "nfn"    - infoc
// - "ngr"    - infoc

// Workspaces
// - r(m+n)     - double
// - w(m+n)     - double
// - e(m+n)     - double
// - ls(m+n)    - int 
// - alp(mlp)   - double
// - lp(mlp)    - int
// - ws(mxws)   - double
// - lws(mxlws) - int
// - cws(m)     - char
// - v(maxg)    - double
// - g(n)       - double

int sci_glcpd_sparse(char * fname)
{
  int n_x_in = 0,        m_x_in = 0,        * pi_x_in        = NULL; double * x_in       = NULL;
  int n_x_upper_in = 0,  m_x_upper_in = 0,  * pi_x_upper_in  = NULL; double * x_upper_in = NULL;
  int n_x_lower_in = 0,  m_x_lower_in = 0,  * pi_x_lower_in  = NULL; double * x_lower_in = NULL;
  int n_a_in  = 0,       m_a_in  = 0,       * pi_a_in        = NULL; double * a_in       = NULL;
  int n_la_in = 0,       m_la_in = 0,       * pi_la_in       = NULL; int    * la_in      = NULL;
  int n_x_out = 0,       m_x_out = 0;      double * x_out    = NULL;
  int n_f_out = 0,       m_f_out = 0;      double * f_out    = NULL;
  int n_ifail_out = 0,   m_ifail_out = 0;  int ifail         = 0;
  int n_tmp_var = 0, m_tmp_var = 0, * pi_tmp_var = NULL; double * tmp_var = NULL;
  int n_grad_c_start = 0, m_grad_c_start = 0, * grad_c_start_addr = NULL; double * grad_c_start = NULL;
  int n_grad_c_index = 0, m_grad_c_index = 0, * grad_c_index_addr = NULL; double * grad_c_index = NULL;
  int n_grad_c_val   = 0, m_grad_c_val   = 0, * grad_c_val_addr   = NULL; double * grad_c_val   = NULL;
  int  * param_in_addr  = NULL;
  int  * param_out_addr = NULL;
  char * cstype         = NULL;
  double rgtol  = 1.0e-4;
  double fmin   = -1.0e6;
  double ainfty = 1e20;
  double ubd    = 1e5;
  int    mlp    = 50;
  int    mxf    = 50;
  int    mxgr   = 1e6;
  int    mxws   = 30000;
  int    mxlws  = 30000;
  int    iprint = 0;
  int    kmax   = 1;
  int    maxg   = 5;
  int    mode   = 0;
  int    Log    = 0;
  int    peq    = 0;
  int    nnza   = 0;
  double tmp_double = 0.0, f_tmp_out = 0.0;
  int tmp_int = 0, tmp_res = 0, grad_type = 0;
  int sci_obj = 0, lhs_obj = 0, rhs_obj = 0;
  // Workspaces
  // - r(m+n)     - double
  // - w(m+n)     - double
  // - e(m+n)     - double
  // - ls(m+n)    - int 
  // - alp(mlp)   - double
  // - lp(mlp)    - int
  // - ws(mxws)   - double
  // - lws(mxlws) - int
  // - cws(m)     - char
  // - v(maxg)    - double
  // - g(n)       - double
  double * ws  = NULL, * r = NULL, * w = NULL, * e = NULL, * alp = NULL, * g = NULL, * a = NULL;
  int    * lws = NULL, * ls = NULL, * lp = NULL;
  double * v   = NULL; int nv = 0;
  double * x_glcpd = NULL, * x_lower_glcpd = NULL, * x_upper_glcpd = NULL;
  int maxu = 0, maxiu = 0, nout = 0, n = 0, m = 0, i = 0, j = 0, k = 0, k_old = 0;
  double * tmp_ptr_dbl = NULL;
  fobj_glcpd_ptr func_obj  = NULL;
  grad_glcpd_ptr func_grad = NULL;
  // Parameters out:
  // - "ifail"
  // - "rgnorm" - infoc
  // - "vstep"  - infoc
  // - "iter"   - infoc
  // - "npv"    - infoc
  // - "nfn"    - infoc
  // - "ngr"    - infoc
  // - "cstype"
  const char * LabelList_out[6] = {"rgnorm", "vstep", "iter", "npv", "nfn", "ngr"};
  int nbvars_old = Nbvars, nbconstr = 0;
  int is_external_func = 0, is_scilab_func = 0;
  SciErr _SciErr;
  
  nout = 0; // stderr
  
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

  _SciErr = getVarAddressFromPosition(pvApiCtx, X_IN, &pi_x_in); GLCPD_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, pi_x_in, &n_x_in, &m_x_in, &x_in); GLCPD_ERROR;
  n = n_x_in*m_x_in;
  
  _SciErr = getVarAddressFromPosition(pvApiCtx, X_LOWER_IN, &pi_x_lower_in); GLCPD_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, pi_x_lower_in, &n_x_lower_in, &m_x_lower_in, &x_lower_in); GLCPD_ERROR;

  x_lower_glcpd = (double*)MALLOC(m_x_lower_in*n_x_lower_in*sizeof(double));
  for(i=0; i<m_x_lower_in*n_x_lower_in; i++) x_lower_glcpd[i] = x_lower_in[i];

  _SciErr = getVarAddressFromPosition(pvApiCtx, X_UPPER_IN, &pi_x_upper_in); GLCPD_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, pi_x_upper_in, &n_x_upper_in, &m_x_upper_in, &x_upper_in); GLCPD_ERROR;

  x_upper_glcpd = (double*)MALLOC(m_x_upper_in*n_x_upper_in*sizeof(double));
  for(i=0; i<m_x_upper_in*n_x_upper_in; i++) x_upper_glcpd[i] = x_upper_in[i];

#ifdef DEBUG
  for(i=0;i<m_x_upper_in*n_x_upper_in; i++) printf("DEBUG: %d - x_lower = %f, x_upper = %f\n", i, x_lower_glcpd[i], x_upper_glcpd[i]);
#endif

  // *************************************************
  // Specification of A and LA in sparse matrix format
  // *************************************************

  // The matrix A contains gradients of the linear terms in the objective function (column 0) and the general constraints (columns 1:m).
  // No explicit reference to simple bound constraints is required in A.
  // The information is set in the parameters a(*) (double precision real) and la(*) (integer).
  //
  // In this sparse format, these vectors have dimension a(1:maxa) and la(0:maxla-1), where maxa is at least nnza (the number of nonzero elements
  // in A), and maxla is at least  nnza+m+3. la(0) and the last m+2 elements in la are pointers.
  //
  // The vectors a(.) and la(.) must be set as follows:
  //
  // a(j) and la(j) for j=1,nnza are set to the values and row indices (resp.) of all the nonzero elements of A. Entries for each column are grouped
  // together in increasing column order. Within each column group, it is not necessary to have the row indices in increasing order.
  //
  // la(0) is a pointer which points to the start of the pointer information in la. la(0) must be set to nnza+1 (or a larger value if it is desired to
  // allow for future increases to nnza).
  //
  // The last m+2 elements of la(.) contain pointers to the first elements in each of the column groupings. Thus la(la(0)+i)) for i=0,m is set to the
  // location in a(.) containing the first nonzero element for column i of A.
  // Also la(la(0)+m+1)) is set to nnza+1 (the first unused location in a(.)).

  _SciErr = getVarAddressFromPosition(pvApiCtx, A_IN, &pi_a_in); GLCPD_ERROR;
  _SciErr = getListItemNumber(pvApiCtx, pi_a_in, &n_a_in); GLCPD_ERROR;

  _SciErr = getVarAddressFromPosition(pvApiCtx, LA_IN, &pi_la_in); GLCPD_ERROR;
  _SciErr = getListItemNumber(pvApiCtx, pi_la_in, &n_la_in); GLCPD_ERROR;

  nbconstr = n_a_in;
  m = nbconstr;

  nnza = 0;
  for(i=0; i<nbconstr; i++)
    {
      _SciErr = getMatrixOfDoubleInList(pvApiCtx, pi_a_in, i+1, &m_tmp_var, &n_tmp_var, &tmp_var);
      nnza += m_tmp_var*n_tmp_var;
    }

  a_in = (double *)MALLOC(nnza*sizeof(double));

  if (a_in==NULL)
    {
      sciprint("%s : Error while allocating a_in (dimension required %d)", fname, nnza);
      FREE(x_lower_glcpd);
      FREE(x_upper_glcpd);
      return 0;
    }

  la_in = (int *)MALLOC((nnza + nbconstr + 3)*sizeof(int));

  if (la_in==NULL)
    {
      sciprint("%s : Error while allocating la_in (dimension required %d)", fname, nnza + nbconstr + 3);
      FREE(x_lower_glcpd);
      FREE(x_upper_glcpd);
      FREE(a_in);
      return 0;
    }

  // Fill a_in
  k = 0;
  for(i=0; i<nbconstr; i++)
    {
      _SciErr = getMatrixOfDoubleInList(pvApiCtx, pi_a_in,  i+1, &m_tmp_var, &n_tmp_var, &tmp_var);
      for(j=0; j<m_tmp_var*n_tmp_var; j++)
	{
	  a_in[k] = tmp_var[j];
	  k++;
	}
    }

#ifdef DEBUG
  for(i=0; i<nnza; i++)
    {
      printf("DEBUG: a_in[%d] = %f\n", i, a_in[i]);
    }
#endif

  // Fill la
  la_in[0] = nnza + 1;

  k = 0; k_old = 1;
  for(i=0; i<nbconstr; i++)
    {
      _SciErr = getMatrixOfDoubleInList(pvApiCtx, pi_la_in,  i+1, &m_tmp_var, &n_tmp_var, &tmp_var);
      for(j=0; j<m_tmp_var*n_tmp_var; j++)
	{
	  la_in[k+1] = tmp_var[j];
	  k++;
	}
      la_in[nnza+2+i] = k_old;
      k_old = k + 1;
    }

  la_in[nnza+4] = nnza + 1;

  // duplicate la_in[nnza+2]
  la_in[nnza+1] = la_in[nnza+2];

#ifdef DEBUG
  for(i=0; i<nnza + nbconstr + 3; i++)
    {
      printf("DEBUG: la_in[%d] = %d\n", i, la_in[i]);
    }
#endif

  // Check upper and lower bounds
  if (m_x_lower_in*n_x_lower_in!=n+m)
    {
      sciprint("%s : Error - x_lower must be of dimension %d", fname, n+m);
      FREE(x_lower_glcpd);
      FREE(x_upper_glcpd);
      FREE(a_in);
      FREE(la_in);
      return 0;
    }

  if (m_x_upper_in*n_x_upper_in!=n+m)
    {
      sciprint("%s : Error - x_upper must be of dimension %d", fname, n+m);
      FREE(x_lower_glcpd);
      FREE(x_upper_glcpd);
      FREE(a_in);
      FREE(la_in);
      return 0;
    }

  // Get Fobj
  // Here, we get either a Scilab script function or a string which gives the name of the C/C++/Fortran function
  // loaded via "link".
  func_obj = (fobj_glcpd_ptr)GetFunctionPtr("glcpd_fobj", FOBJ_IN, FTab_glcpd_function, (voidf)C2F(sparse_glcpd_functions), &sci_obj, &lhs_obj, &rhs_obj);

  if ((lhs_obj!=2)&&(rhs_obj!=1))
    {
      sciprint("%s : Error - function must have 2 outputs parameters (f and c) and 1 input parameter (x)", fname);
      FREE(x_lower_glcpd);
      FREE(x_upper_glcpd);
      FREE(a_in);
      FREE(la_in);
      return 0;
    }

  param_fobj.fobj_type = 0;
  param_fobj.n_x       = n;
  param_fobj.sci_obj   = sci_obj;
  param_fobj.lhs_obj   = lhs_obj;
  param_fobj.rhs_obj   = rhs_obj;
  param_fobj.stack_pos = Rhs + 1; // We preallocate the objective function value on position 3 on the stack
                                  // So, scifunction will create all the needed variables on LAST_PARAM_OUT + 1 to
                                  // avoid the destruction of this preallocated output variable

  if (func_obj==(fobj_glcpd_ptr)0) 
    {
      sciprint("%s : Error - Last argument must be a pointer to a scilab function", fname);
      FREE(x_lower_glcpd);
      FREE(x_upper_glcpd);
      FREE(a_in);
      FREE(la_in);
      return 0;
    }

  if (func_obj!=(voidf)C2F(sparse_glcpd_functions))
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
  func_grad = (grad_glcpd_ptr)GetFunctionPtr("glcpd_grad", GRAD_IN, FTab_glcpd_function, (voidf)C2F(sparse_glcpd_gradients), &sci_obj, &lhs_obj, &rhs_obj);
  
  if ((lhs_obj!=1)&&(rhs_obj!=1))
    {
      sciprint("%s : Error - gradient must have 1 output parameters (a) and 1 input parameters (x)", fname);
      FREE(x_lower_glcpd);
      FREE(x_upper_glcpd);
      FREE(a_in);
      FREE(la_in);
      return 0;
    }

  param_grad.fobj_type = 0;
  param_grad.n_x       = n;
  param_grad.sci_obj   = sci_obj;
  param_grad.lhs_obj   = lhs_obj;
  param_grad.rhs_obj   = rhs_obj;
  param_grad.stack_pos = Rhs + 1; // We preallocate the objective function value on position 3 on the stack
                                  // So, scifunction will create all the needed variables on LAST_PARAM_OUT + 1 to
                                  // avoid the destruction of this preallocated output variable

  if (func_grad==(grad_glcpd_ptr)0) 
    {
      sciprint("%s : Error - Last argument must be a pointer to a scilab function", fname);
      FREE(x_lower_glcpd);
      FREE(x_upper_glcpd);
      FREE(a_in);
      FREE(la_in);
      return 0;
    }

  if (func_grad!=(voidf)C2F(sparse_glcpd_gradients))
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
      FREE(x_lower_glcpd);
      FREE(x_upper_glcpd);
      FREE(a_in);
      FREE(la_in);
      return 0;
    }

  initPList(pvApiCtx, PARAMS_IN, &param_in_addr);
  if (!checkPList(pvApiCtx, param_in_addr))
    {
      Scierror(999, "%s: argument n° %d is not a plist\n", fname, PARAMS_IN);
      FREE(x_lower_glcpd);
      FREE(x_upper_glcpd);
      FREE(a_in);
      FREE(la_in);
      return 0;
    }

  getDoubleInPList(pvApiCtx, param_in_addr, "rgtol", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) rgtol = tmp_double;

  getDoubleInPList(pvApiCtx, param_in_addr, "ainfty", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) ainfty = tmp_double;

  getDoubleInPList(pvApiCtx, param_in_addr, "ubd", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) ubd = tmp_double;

  getDoubleInPList(pvApiCtx, param_in_addr, "fmin", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) fmin = tmp_double;

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

  getIntInPList(pvApiCtx, param_in_addr, "mode", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) mode = tmp_int;

  getIntInPList(pvApiCtx, param_in_addr, "mxgr", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) mxgr = tmp_int;

  getIntInPList(pvApiCtx, param_in_addr, "mxws", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) mxws = tmp_int;

  getIntInPList(pvApiCtx, param_in_addr, "mxlws", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) mxlws = tmp_int;

  getIntInPList(pvApiCtx, param_in_addr, "nout", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) nout = tmp_int;

  getIntInPList(pvApiCtx, param_in_addr, "peq", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) peq = tmp_int;

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
#endif

  /**********************
   * Pre-Initialisation *
   **********************/

  // set stride in lws(maxiu+1) and constant elements of a(*) in ws(maxu+1) on
  
  // La taille de ws et lws sont reglees dans le common wsc via les variables mxws mxlws

  C2F(wsc).mxws = mxws;
  ws = (double *)MALLOC(C2F(wsc).mxws*sizeof(double));
  for(i=0;i<C2F(wsc).mxws;i++) ws[i] = 0.0;

  C2F(wsc).mxlws = mxlws;
  lws = (int *)MALLOC(C2F(wsc).mxlws*sizeof(int));
  for(i=0;i<C2F(wsc).mxlws;i++) lws[i] = 0;
  lws[maxiu+0] = n; // stride initialization

  cstype = (char *)MALLOC((m+1)*sizeof(char));
  for(i=0;i<m;i++) cstype[i] = '\0';
  
  r = (double *)MALLOC((m+n)*sizeof(double));
  for(i=0;i<m+n;i++) r[i] = 0.0;

  w = (double *)MALLOC((m+n)*sizeof(double));
  for(i=0;i<m+n;i++) w[i] = 0.0;

  e = (double *)MALLOC((m+n)*sizeof(double));
  for(i=0;i<m+n;i++) e[i] = 0.0;

  alp = (double *)MALLOC((m+n)*sizeof(double));
  for(i=0;i<m+n;i++) alp[i] = 0.0;

  g = (double *)MALLOC(n*sizeof(double));
  for(i=0;i<n;i++) g[i] = 0.0;

  ls = (int *)MALLOC((m+n)*sizeof(int));
  for(i=0;i<m+n;i++) ls[i] = 0.0;

  lp = (int *)MALLOC((m+n)*sizeof(int));
  for(i=0;i<m+n;i++) lp[i] = 0.0;

  ////////////////////
  // Initialization //
  ////////////////////

  // Allouer l'espace pour
  // ws  -> ws(mxws)
  // lws -> ws(mxlws)
  // v   -> maxg (then nv = maxg - otherwise, set nv=1 and v(1) = 1)

  v = (double *)MALLOC((maxg+1)*sizeof(double));
  maxg = MIN(n,maxg)+1;

  for(i=0;i<maxg; i++) v[i] = 0.0;
  nv = 0;

  x_glcpd = (double *)MALLOC((n+m+1)*sizeof(double));
  for(i=n;i<(n+m+1);i++) x_glcpd[i] = 0.0; // initialize de x workspace
  for(i=0;i<n;i++)       x_glcpd[i] = x_in[i];

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

#ifdef DEBUG
  sciprint("DEBUG: defaultc common - ainfty = %f\n", C2F(defaultc).ainfty);
  sciprint("DEBUG: defaultc common - ubd    = %f\n", C2F(defaultc).ubd);
  sciprint("DEBUG: defaultc common - mlp    = %d\n", C2F(defaultc).mlp);
  sciprint("DEBUG: defaultc common - mxf    = %d\n", C2F(defaultc).mxf);
#endif

  // Parameters for glcpd:
  // common/epsc/eps,tol,emin
  // common/repc/sgnf,nrep,npiv,nres
  // common/refactorc/mc,mxmc

  /*
   * Calling glcpd
   */

  if (!is_external_func)
    {
      C2F(sci_sp_glcpd)(&n,&m,&k,&kmax,&maxg,a_in,la_in,x_glcpd,x_lower_glcpd,x_upper_glcpd,&f_tmp_out,&fmin,g,r,w,e,ls,alp,lp,&mlp,&peq,ws,lws,cstype,v,&nv,&rgtol,&mode,&ifail,
			&mxgr,&iprint,&nout,C2F(sparse_glcpd_functions),C2F(sparse_glcpd_gradients));
    }
  else
    {
      C2F(sci_sp_glcpd)(&n,&m,&k,&kmax,&maxg,a_in,la_in,x_glcpd,x_lower_glcpd,x_upper_glcpd,&f_tmp_out,&fmin,g,r,w,e,ls,alp,lp,&mlp,&peq,ws,lws,cstype,v,&nv,&rgtol,&mode,&ifail,
			&mxgr,&iprint,&nout,(fobj_glcpd_ptr)param_fobj.function,(grad_glcpd_ptr)param_grad.function);
    }

  //////////////////////////
  // return the variables //
  //////////////////////////

  n_x_out = n_x_in; m_x_out = m_x_in;
  _SciErr = allocMatrixOfDouble(pvApiCtx, X_OUT, n_x_out, m_x_out, &x_out); GLCPD_ERROR;
  for(i=0; i<n_x_out*m_x_out; i++) x_out[i] = x_glcpd[i];

  // f_out
  n_f_out = 1; m_f_out = 1;
  _SciErr = allocMatrixOfDouble(pvApiCtx, F_OUT, n_f_out, m_f_out, &f_out); GLCPD_ERROR;
  f_out[0] = f_tmp_out;

  // ifail_out   outcome of the process
  // 0   = solution obtained
  // 1   = unbounded problem (f(x)<fmin has occurred: note grad is not evaluated in this case)
  // 2   = bl(i) > bu(i) for some i
  // 3   = infeasible problem detected in Phase 1
  // 4   = line search cannot improve f (possibly increase rgtol)
  // 5   = mxgr gradient calls exceeded (this test is only carried out at the start of each iteration)
  // 6   = incorrect setting of m, n, kmax, maxg, mlp, mode or tol
  // 7   = not enough space in ws or lws
  // 8   = not enough space in lp (increase mlp)
  // 9   = dimension of reduced space too large (increase kmax)
  // 10  = maximum number of unsuccessful restarts taken
  // >10 = possible use by later sparse matrix codes

  n_ifail_out = 1; m_ifail_out = 1;
  _SciErr = allocMatrixOfDouble(pvApiCtx, IFAIL_OUT, n_ifail_out, m_ifail_out, &tmp_ptr_dbl); GLCPD_ERROR;
  tmp_ptr_dbl[0] = (double)ifail;

#ifdef DEBUG
  sciprint("DEBUG: X_OUT          = %d\n", X_OUT);
  sciprint("DEBUG: F_OUT          = %d\n", F_OUT);
  sciprint("DEBUG: IFAIL_OUT      = %d\n", IFAIL_OUT);
  sciprint("DEBUG: PARAMS_OUT     = %d\n", PARAMS_OUT);
  sciprint("DEBUG: LAST_PARAM_OUT = %d\n", LAST_PARAM_OUT);
#endif

  // LabelList_out[6] = {"rgnorm", "vstep", "iter", "npv", "nfn", "ngr"};

  _SciErr = createPList(pvApiCtx, PARAMS_OUT, &param_out_addr, (char **)LabelList_out, 6); GLCPD_ERROR;

  _SciErr = createDoubleInPList(pvApiCtx, PARAMS_OUT, param_out_addr, "rgnorm", C2F(infoc).rgnorm); GLCPD_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, PARAMS_OUT, param_out_addr, "vstep",  C2F(infoc).vstep);  GLCPD_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, PARAMS_OUT, param_out_addr, "iter",   C2F(infoc).iter);   GLCPD_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, PARAMS_OUT, param_out_addr, "npv",    C2F(infoc).npv);    GLCPD_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, PARAMS_OUT, param_out_addr, "nfn",    C2F(infoc).nfn);    GLCPD_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, PARAMS_OUT, param_out_addr, "ngr",    C2F(infoc).ngr);    GLCPD_ERROR;

  LhsVar(1) = X_OUT;
  LhsVar(2) = F_OUT;
  LhsVar(3) = IFAIL_OUT;
  LhsVar(4) = PARAMS_OUT;

  // YC: Plantage ici sinon .... 

  if (x_lower_glcpd) FREE(x_lower_glcpd);
  if (x_upper_glcpd) FREE(x_upper_glcpd);
  if (ws)            FREE(ws);
  if (lws)           FREE(lws);
  if (cstype)        FREE(cstype);
  if (v)             FREE(v);
  if (x_glcpd)       FREE(x_glcpd);
  if (r)             FREE(r);
  if (w)             FREE(w);
  if (e)             FREE(e);
  if (alp)           FREE(alp);
  if (ls)            FREE(ls);
  if (lp)            FREE(lp);
  if (g)             FREE(g);
  if (a)             FREE(a);
  if (a_in)          FREE(a_in);
  if (la_in)         FREE(la_in);

  return 0;
}

/***************************
 * Functions and Gradients *
 ***************************/

// This function is a wrapper which allows to call a Scilab script function like a "normal" C/C++/Fortran function
void C2F(sparse_glcpd_functions)(int * n, double * x, double * f, double * ws, int * lws, char * cws)
{
  int n_x = 0;
  int * tmp_addr = NULL;
  int sci_obj, lhs_obj, rhs_obj;
  int stack_pos, i;
  int nbvars_old = Nbvars;
  double * tmp_var = NULL, tmp_val = 0.0;
  double f_out = 0.0;
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

      Nbvars = stack_pos + MAX(rhs_obj,lhs_obj) + 1;
      
#ifdef DEBUG
      printf("DEBUG: n = %d, m = %d, n_x = %d\n", *n, *m, n_x);
      for(i=0; i<n_x; i++) printf("DEBUG: x[%d] = %f\n", i, x[i]);
#endif
      _SciErr = createMatrixOfDouble(pvApiCtx, stack_pos+0, n_x, 1, x); GLCPD_ERROR_NORETURN;

      if (!C2F(scifunction)(&stack_pos, &sci_obj, &lhs_obj, &rhs_obj))
	{
	  Scierror(999,"sparse glcpd: error when calling objective function\n");
	  Nbvars = nbvars_old;
	  return;
	}
      
      if (Err>0) 
	{
	  Scierror(999,"sparse glcpd: error when calling objective function\n");
	  Nbvars = nbvars_old;
	  return;
	} 
      
      // Get F
      _SciErr = getVarAddressFromPosition(pvApiCtx, stack_pos+0, &tmp_addr); GLCPD_ERROR_NORETURN;
      getScalarDouble(pvApiCtx, tmp_addr, &f_out);
      *f = f_out;
#ifdef DEBUG
      printf("DEBUG: f = %f\n", *f);
#endif

      Nbvars = nbvars_old;
    }
  else
    {
      (*param_fobj.function)(n, x, f, NULL, NULL, NULL);
      //Scierror(999,"sparse glcpd: call to C functions not implemented yet\n");
      //return;
    }

#ifdef DEBUG
  printf("DEBUG: End of Functions\n");
#endif
}

// This function is a wrapper which allows to call a Scilab script function like a "normal" C/C++/Fortran function
void C2F(sparse_glcpd_gradients)(int * n, double * x, double * g, double * ws, int * lws, char * cws)
{
  int n_grad_f = 0, m_grad_f = 0, * grad_f_addr = NULL; double * grad_f = NULL;
  int n_x = 0;
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
      
      Nbvars = stack_pos + MAX(rhs_obj,lhs_obj) + 1;

#ifdef DEBUG
      for(i=0; i<n_x; i++) printf("DEBUG: x[%d] = %f\n", i, x[i]);
#endif

      _SciErr = createMatrixOfDouble(pvApiCtx, stack_pos+0, n_x, 1, x); GLCPD_ERROR_NORETURN;
      
      if (!C2F(scifunction)(&stack_pos, &sci_obj, &lhs_obj, &rhs_obj))
	{
	  Scierror(999,"sparse glcpd: error when calling gradient function\n");
	  Nbvars = nbvars_old;
	  return;
	}
      
      if (Err>0) 
	{
	  Scierror(999,"sparse glcpd: error when calling gradient function\n");
	  Nbvars = nbvars_old;
	  return;
	} 

      _SciErr = getVarAddressFromPosition(pvApiCtx, stack_pos+0, &grad_f_addr); GLCPD_ERROR_NORETURN;
      _SciErr = getMatrixOfDouble(pvApiCtx, grad_f_addr, &m_grad_f, &n_grad_f, &grad_f); GLCPD_ERROR_NORETURN;
#ifdef DEBUG
      printf("DEBUG: m_grad_f = %d, n_grad_f = %d\n", m_grad_f, n_grad_f);
#endif

      for(i=0;i<n_x;i++)
	{
	  g[i] = grad_f[i];
	}

      Nbvars = nbvars_old;
    }
  else
    {
      (*param_grad.gradient)(n, x, g, NULL, NULL, NULL);
      //Scierror(999,"sparse glcpd: call to C functions not implemented yet\n");
      //return;
    }

#ifdef DEBUG
  printf("DEBUG: End of Gradients\n");
#endif
}
