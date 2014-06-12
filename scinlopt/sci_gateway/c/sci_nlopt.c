#include <stack-c.h>
#include <MALLOC.h>
#include <api_scilab.h>
#include <api_parameters.h>

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <setjmp.h>

#include <nlopt.h>

//#define DEBUG 1

#define X_IN          1
#define F_IN          2
#define CINEQ_IN      3
#define CEQ_IN        4
#define LOWER_IN      5
#define UPPER_IN      6
#define ITER_IN       7
#define PARAMS_IN     8
#define LAST_ARG      8
#define X_OPT_OUT     Rhs+1
#define F_OPT_OUT     Rhs+2
#define STATUS_OUT    Rhs+3
#define METH_NAME_OUT Rhs+4

#define OK   1
#define FAIL 0

#define _(A) A
#define MAX(A,B) ((A<B)?B:A)

#define NLOPT_ERROR if(_SciErr.iErr)		\
    {						\
      printError(&_SciErr, 0);			\
      return 0;					\
    }

/******************************
 * Work around for the        *
 * AddFunctionInTable problem *
 ******************************/

typedef void (*voidf)();
typedef struct {
  char *name;
  voidf f;
} FTAB;

extern voidf AddFunctionInTable(char *name, int *rep, FTAB *table);  

/******************************************
 * General functions for scilab interface *
 ******************************************/

// SearchInDynLinks: a scilab function which tries to find a C/C++/Fortran function loaded via link
// SearchComp: a function which tries to find a C/C++/Fortran function stored in the FTAB table (ours is empty)
// SetFunction: this function looks inside the table or inside Scilab for dynamically loaded functions
// sciobj: a wrapper function which has the same prototype as a dynamically loaded C/C++/Fortran function and which is used
//         to call a Scilab script function like a C/C++/Fortran function
// Emptyfunc: an empty function which is used when no dynamically functions has been found
// GetFunctionPtr: the main function. Get the parameter on the stack. If it's a Scilab function, then it returns a pointer to
//                 sciobj. Otherwise, returns a pointer to the dynamically loaded function

extern int   SearchInDynLinks(char *op, void (**realop) ());
static int   SearchComp(FTAB *Ftab, char *op, void (**realop) ( ));  
static voidf SetFunction(char *name, int *rep, FTAB *table);  
static int   sciobj(int size_x, int size_constr, double * x, double * f, double * con, void * param);
static void  Emptyfunc(void) {} ;
voidf        GetFunctionPtr(char *, int, FTAB *, voidf, int *, int*, int*);

struct param_obj {
  int sci_obj;
  int lhs_obj;
  int rhs_obj;
  int stack_pos;
};

/***********************
 *Definition of FTables*
 ***********************/

static FTAB FTab_nlopt_function[] = {{(char *) 0, (voidf) 0}};

/************************
 * The Scilab interface *
 ************************/

static jmp_buf nlopt_function_env; 

double sci_fobj(int size_x, const double * x, double * gradient, double * con, void * param);
void sci_constr_ineq(unsigned int size_constr, double * result, unsigned int size_x, const double * x, double * gradient, void * param);
void sci_constr_eq(unsigned int size_constr, double * result, unsigned int size_x, const double * x, double * gradient, void * param);

static struct param_obj param_fobj;
static struct param_obj param_cineq;
static struct param_obj param_ceq;

int sci_nlopt(char * fname)
{
  int m_x,         n_x,        * x_addr = NULL, l_x;
  int m_fobj  = 1, n_fobj  = 1,                 l_fobj;
  int m_cineq = 1, n_cineq = 1,                 l_cineq;
  int m_ceq   = 1, n_ceq   = 1,                 l_ceq;
  double * pdbl_x = NULL, * pdbl_xout = NULL;
  double * pdbl_nbiter = NULL;
  double * pdbl_lower = NULL;
  double * pdbl_upper = NULL;
  double * ceq_tol = NULL, * cineq_tol = NULL, * tmp_ptr_double = NULL;
  double f_out = 0, tmp_double;
  int * nbiter_addr = NULL, m_nbiter, n_nbiter, * param_in_addr = NULL;
  int * upper_addr = NULL, * lower_addr = NULL, m_upper, n_upper, m_lower, n_lower;
  int i, size_x, maxiter, status_out, tmp_int, tmp_res, Log = 0;
  int sci_obj, lhs_obj, rhs_obj, status;
  int real_size, method_number, nb_cineq, nb_ceq;
  nlopt_opt Method;
  nlopt_func func_obj;
  nlopt_mfunc func_cineq, func_ceq;
  nlopt_result res_nlopt;
  SciErr _SciErr;

  CheckRhs(LAST_ARG,LAST_ARG);
  CheckLhs(2,4);

  // Get X
  _SciErr = getVarAddressFromPosition(pvApiCtx, X_IN, &x_addr); NLOPT_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, x_addr, &m_x, &n_x, &pdbl_x); NLOPT_ERROR;
  size_x = m_x*n_x;

  // Get Fobj
  // Here, we get either a Scilab script function or a string which gives the name of the C/C++/Fortran function
  // loaded via "link".
  func_obj = (nlopt_func)GetFunctionPtr("nlopt_fobj", F_IN, FTab_nlopt_function, (voidf)sci_fobj, &sci_obj, &lhs_obj, &rhs_obj);

  param_fobj.sci_obj   = sci_obj;
  param_fobj.lhs_obj   = lhs_obj;
  param_fobj.rhs_obj   = rhs_obj;
  param_fobj.stack_pos = Rhs+2; // We preallocate the objective function value on position 3 on the stack
                                // So, scifunction will create all the needed variables on Rhs+2 == 4 to
                                // avoid the destruction of this preallocated output variable

  if (func_obj==(nlopt_func)0) 
    {
      sciprint("%s : Error - Last argument must be a pointer to a scilab function", fname);
      return 0;
    }

  // Get Cineq
  // Here, we get either a Scilab script function or a string which gives the name of the C/C++/Fortran function
  // loaded via "link".
  if (GetType(CINEQ_IN)==sci_matrix)
    {
      func_cineq = NULL;
      param_cineq.sci_obj   = -1;
      param_cineq.lhs_obj   = 0;
      param_cineq.rhs_obj   = 0;
      param_cineq.stack_pos = Rhs+2; // We preallocate the objective function value on position 3 on the stack
    }
  else
    {
      func_cineq = (nlopt_mfunc)GetFunctionPtr("nlopt_cineq", CINEQ_IN, FTab_nlopt_function, (voidf)sci_constr_ineq, &sci_obj, &lhs_obj, &rhs_obj);
      
      param_cineq.sci_obj   = sci_obj;
      param_cineq.lhs_obj   = lhs_obj;
      param_cineq.rhs_obj   = rhs_obj;
      param_cineq.stack_pos = Rhs+2; // We preallocate the objective function value on position 3 on the stack
                                     // So, scifunction will create all the needed variables on Rhs+2 == 4 to
                                     // avoid the destruction of this preallocated output variable
    }

  // Get Ceq
  // Here, we get either a Scilab script function or a string which gives the name of the C/C++/Fortran function
  // loaded via "link".
  if (GetType(CEQ_IN)==sci_matrix)
    {
      func_ceq = NULL;
      param_ceq.sci_obj   = -1;
      param_ceq.lhs_obj   = 0;
      param_ceq.rhs_obj   = 0;
      param_ceq.stack_pos = Rhs+2; // We preallocate the objective function value on position 3 on the stack
    }
  else
    {
      func_ceq = (nlopt_mfunc)GetFunctionPtr("nlopt_ceq", CEQ_IN, FTab_nlopt_function, (voidf)sci_constr_eq, &sci_obj, &lhs_obj, &rhs_obj);
      
      param_ceq.sci_obj   = sci_obj;
      param_ceq.lhs_obj   = lhs_obj;
      param_ceq.rhs_obj   = rhs_obj;
      param_ceq.stack_pos = Rhs+2; // We preallocate the objective function value on position 3 on the stack
                                   // So, scifunction will create all the needed variables on Rhs+2 == 4 to
                                   // avoid the destruction of this preallocated output variable
    }

  // Lower
  _SciErr = getVarAddressFromPosition(pvApiCtx, LOWER_IN, &lower_addr); NLOPT_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, lower_addr, &m_lower, &n_lower, &pdbl_lower); NLOPT_ERROR;
  if (m_lower*n_lower!=size_x)
    {
      sciprint("%s : Error - lower vector must be of the same size as x", fname);
      return 0;
    }

  // Upper
  _SciErr = getVarAddressFromPosition(pvApiCtx, UPPER_IN, &upper_addr); NLOPT_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, upper_addr, &m_upper, &n_upper, &pdbl_upper); NLOPT_ERROR;
  if (m_lower*n_lower!=size_x)
    {
      sciprint("%s : Error - upper vector must be of the same size as x", fname);
      return 0;
    }

  // Nb iters
  _SciErr = getVarAddressFromPosition(pvApiCtx, ITER_IN, &nbiter_addr); NLOPT_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, nbiter_addr, &m_nbiter, &n_nbiter, &pdbl_nbiter); NLOPT_ERROR;
  maxiter = (int)pdbl_nbiter[0];

  initPList(pvApiCtx, PARAMS_IN, &param_in_addr);
  if (!checkPList(pvApiCtx, param_in_addr))
    {
      Scierror(999, "%s: argument n° %d is not a plist\n", fname, PARAMS_IN);

      return 0;
    }

  getIntInPList(pvApiCtx, param_in_addr, "srand", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) nlopt_srand(tmp_int);

  getIntInPList(pvApiCtx, param_in_addr, "nb_ceq", &nb_ceq, &tmp_res, 0, Log, CHECK_NONE);
  getIntInPList(pvApiCtx, param_in_addr, "nb_cineq", &nb_cineq, &tmp_res, 0, Log, CHECK_NONE);

  getIntInPList(pvApiCtx, param_in_addr, "method", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) 
    {
      Method = nlopt_create((nlopt_algorithm)tmp_int, (unsigned int)size_x);

      if (Method==NULL)
        {
          Scierror(999,"%s: Problem while allocating the method. Maybe the selected method is not available\n", fname);

          return 0;
        }

      method_number = tmp_int;
    }
  else
    {
      Method = nlopt_create((nlopt_algorithm)0, (unsigned int)size_x);
      method_number = 0;
    }

  nlopt_set_lower_bounds(Method, pdbl_lower);
  nlopt_set_upper_bounds(Method, pdbl_upper);

  getIntInPList(pvApiCtx, param_in_addr, "max", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if ((tmp_res!=-1)&&(tmp_int!=0)) 
    {
      res_nlopt = nlopt_set_max_objective(Method, func_obj, NULL);
    }
  else
    {
      res_nlopt = nlopt_set_min_objective(Method, func_obj, NULL);
    }

  if ((func_cineq!=NULL)&&(nb_cineq!=0))
    {
      tmp_ptr_double = (double *)MALLOC(nb_cineq*sizeof(double));
      getColVectorOfDoubleInPList(pvApiCtx, param_in_addr, "cineq_tol", tmp_ptr_double, &tmp_res, 1e-6, nb_cineq, &real_size, Log, CHECK_NONE);
      res_nlopt = nlopt_remove_inequality_constraints(Method);
      res_nlopt = nlopt_add_inequality_mconstraint(Method, nb_cineq, func_cineq, NULL, tmp_ptr_double);
      FREE(tmp_ptr_double);

      if (res_nlopt<0)
        {
          createEmptyMatrix(pvApiCtx, Rhs+1);
          LhsVar(1) = Rhs+1;
          createEmptyMatrix(pvApiCtx, Rhs+2);
          LhsVar(2) = Rhs+2;
          if (Lhs>=3)
            {
              createScalarDouble(pvApiCtx, Rhs+3, (double)res_nlopt); NLOPT_ERROR;
              LhsVar(3) = Rhs+3;
            }

          if (Lhs>=4)
            {
              createSingleString(pvApiCtx, Rhs+4, (char *)nlopt_algorithm_name((nlopt_algorithm)method_number));
              LhsVar(4) = Rhs+4;
            }

          // Now release the structure
          nlopt_destroy(Method);
          
          return 0;
        }
    }

  if ((func_ceq!=NULL)&&(nb_ceq!=0))
    {
      tmp_ptr_double = (double *)MALLOC(nb_ceq*sizeof(double));
      getColVectorOfDoubleInPList(pvApiCtx, param_in_addr, "ceq_tol", tmp_ptr_double, &tmp_res, 1e-6, nb_ceq, &real_size, Log, CHECK_NONE);
      res_nlopt = nlopt_remove_equality_constraints(Method);
      res_nlopt = nlopt_add_equality_mconstraint(Method, nb_ceq, func_ceq, NULL, tmp_ptr_double);
      FREE(tmp_ptr_double);

      if (res_nlopt<0)
        {
          createEmptyMatrix(pvApiCtx, Rhs+1);
          LhsVar(1) = Rhs+1;
          createEmptyMatrix(pvApiCtx, Rhs+2);
          LhsVar(2) = Rhs+2;
          if (Lhs>=3)
            {
              createScalarDouble(pvApiCtx, Rhs+3, (double)res_nlopt); NLOPT_ERROR;
              LhsVar(3) = Rhs+3;
            }

          if (Lhs>=4)
            {
              createSingleString(pvApiCtx, Rhs+4, (char *)nlopt_algorithm_name((nlopt_algorithm)method_number));
              LhsVar(4) = Rhs+4;
            }

          // Now release the structure
          nlopt_destroy(Method);
          
          return 0;
        }
    }

  nlopt_set_maxeval(Method, maxiter);

  getDoubleInPList(pvApiCtx, param_in_addr, "stopval", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) nlopt_set_stopval(Method, tmp_double);

  getDoubleInPList(pvApiCtx, param_in_addr, "ftol_rel", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) nlopt_set_ftol_rel(Method, tmp_double);

  getDoubleInPList(pvApiCtx, param_in_addr, "ftol_abs", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) nlopt_set_ftol_abs(Method, tmp_double);

  getDoubleInPList(pvApiCtx, param_in_addr, "xtol_rel", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) nlopt_set_xtol_rel(Method, tmp_double);

  getDoubleInPList(pvApiCtx, param_in_addr, "xtol_abs1", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) nlopt_set_xtol_abs1(Method, tmp_double);

  tmp_ptr_double = (double *)MALLOC(size_x*sizeof(double));
  getColVectorOfDoubleInPList(pvApiCtx, param_in_addr, "xtol_abs", tmp_ptr_double, &tmp_res, 0, size_x, &real_size, Log, CHECK_SIZE);
  if (tmp_res!=-1) nlopt_set_xtol_abs(Method, tmp_ptr_double);
  FREE(tmp_ptr_double);

  getDoubleInPList(pvApiCtx, param_in_addr, "maxtime", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) nlopt_set_maxtime(Method, tmp_double);

  getIntInPList(pvApiCtx, param_in_addr, "force_stop", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) nlopt_set_force_stop(Method, tmp_int);

  getIntInPList(pvApiCtx, param_in_addr, "population", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) nlopt_set_population(Method, tmp_int);

  tmp_ptr_double = (double *)MALLOC(size_x*sizeof(double));
  getColVectorOfDoubleInPList(pvApiCtx, param_in_addr, "default_initial_step", tmp_ptr_double, &tmp_res, 0, size_x, &real_size, Log, CHECK_SIZE);
  if (tmp_res!=-1) nlopt_set_default_initial_step(Method, tmp_ptr_double);
  FREE(tmp_ptr_double);

  tmp_ptr_double = (double *)MALLOC(size_x*sizeof(double));
  getColVectorOfDoubleInPList(pvApiCtx, param_in_addr, "initial_step", tmp_ptr_double, &tmp_res, 0, size_x, &real_size, Log, CHECK_SIZE);
  if (tmp_res!=-1) nlopt_set_initial_step(Method, tmp_ptr_double);
  FREE(tmp_ptr_double);

  getDoubleInPList(pvApiCtx, param_in_addr, "initial_step1", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) nlopt_set_initial_step1(Method, tmp_double);

  // YC: for some methods, we need to specify the local optimizer ...
  // YC: TO BE DONE
  // nlopt_result nlopt_set_local_optimizer(nlopt_opt opt, const nlopt_opt local_opt);

  _SciErr = allocMatrixOfDouble(pvApiCtx, X_OPT_OUT, size_x, 1, &pdbl_xout); NLOPT_ERROR;
  for(i=0;i<size_x;i++) pdbl_xout[i] = pdbl_x[i];

  status = nlopt_optimize(Method, pdbl_xout, &f_out);

  createScalarDouble(pvApiCtx, F_OPT_OUT, f_out); NLOPT_ERROR;
  if (Lhs>=3)
    {
      createScalarDouble(pvApiCtx, STATUS_OUT, status); NLOPT_ERROR;
    }

  if (Lhs>=4)
    {
      createSingleString(pvApiCtx, METH_NAME_OUT, (char *)nlopt_algorithm_name((nlopt_algorithm)method_number));
    }

  LhsVar(1) = F_OPT_OUT;
  LhsVar(2) = X_OPT_OUT;
  if (Lhs>=3)
    {
      LhsVar(3) = STATUS_OUT;
    }
  if (Lhs>=4)
    {
      LhsVar(4) = METH_NAME_OUT;
    }

  // Now release the structure
  nlopt_destroy(Method);

  return 0;
}

// This function is a wrapper which allows to call a Scilab script function like a "normal" 
// C/C++/Fortran function
double sci_fobj(int size_x, const double * x, double * gradient, double * con, void * param)
{
  int n_x = size_x, m_x = 1;
  int n_tmp, m_tmp, * tmp_addr = NULL;
  int sci_obj, lhs_obj, rhs_obj;
  int stack_pos, i, output_state;
  int rhs_old = Rhs, nbvars_old = Nbvars, Index = 0;
  double * x_sci = NULL;
  double * tmp_var = NULL;
  double tmp_val, f_out;
  SciErr _SciErr;

  sci_obj   = param_fobj.sci_obj;
  lhs_obj   = param_fobj.lhs_obj;
  rhs_obj   = param_fobj.rhs_obj;
  stack_pos = param_fobj.stack_pos;

  Nbvars = stack_pos + MAX(rhs_obj,lhs_obj);

  _SciErr = allocMatrixOfDouble(pvApiCtx, stack_pos+0, n_x, m_x, &x_sci); NLOPT_ERROR;
  for(i=0;i<size_x;i++) x_sci[i] = x[i];

  _SciErr = allocMatrixOfDouble(pvApiCtx, stack_pos+1, 1, 1, &tmp_var); NLOPT_ERROR;

  SciFunction(&stack_pos, &sci_obj, &lhs_obj, &rhs_obj);

  if (Err>0) 
    {
      Scierror(999,"nlopt: error when calling objective function\n");
      Nbvars = nbvars_old;
      return;
    } 

  _SciErr = getVarAddressFromPosition(pvApiCtx, stack_pos+0, &tmp_addr); NLOPT_ERROR;
  getScalarDouble(pvApiCtx, tmp_addr, &f_out);

  if (gradient!=NULL)
    {
      _SciErr = getVarAddressFromPosition(pvApiCtx, stack_pos+1, &tmp_addr); NLOPT_ERROR;
      _SciErr = getMatrixOfDouble(pvApiCtx, tmp_addr, &n_tmp, &m_tmp, &tmp_var); NLOPT_ERROR;
      // YC: check gradient size
      for(i=0;i<size_x;i++) gradient[i] = tmp_var[i];
    }

  Nbvars = nbvars_old;

  return f_out;
}

// This function is a wrapper which allows to call a Scilab script function like a "normal" 
// C/C++/Fortran function
void sci_constr_ineq(unsigned int size_constr, double * result, unsigned int size_x, const double * x, double * gradient, void * param)
{
  int n_x = size_x, m_x = 1;
  int n_tmp, m_tmp, * tmp_addr = NULL;
  int sci_obj, lhs_obj, rhs_obj;
  int stack_pos, i, output_state;
  int rhs_old = Rhs, nbvars_old = Nbvars, Index = 0;
  double * x_sci = NULL;
  double * tmp_var = NULL;
  double tmp_val, f_out;
  SciErr _SciErr;

  sci_obj   = param_cineq.sci_obj;
  lhs_obj   = param_cineq.lhs_obj;
  rhs_obj   = param_cineq.rhs_obj;
  stack_pos = param_cineq.stack_pos;

  Nbvars = stack_pos + MAX(rhs_obj,lhs_obj) + 5;

  _SciErr = allocMatrixOfDouble(pvApiCtx, stack_pos+0, n_x, m_x, &x_sci); NLOPT_ERROR;
  for(i=0;i<size_x;i++) 
    {
      x_sci[i] = x[i];
    }

  _SciErr = allocMatrixOfDouble(pvApiCtx, stack_pos+1, 1, 1, &tmp_var);

  SciFunction(&stack_pos, &sci_obj, &lhs_obj, &rhs_obj);

  if (Err>0) 
    {
      Scierror(999,"nlopt: error when calling inequality constraint function\n");
      Nbvars = nbvars_old;
      return;
    } 

  _SciErr = getVarAddressFromPosition(pvApiCtx, stack_pos+0, &tmp_addr); NLOPT_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, tmp_addr, &m_tmp, &n_tmp, &tmp_var); NLOPT_ERROR;

  for(i=0;i<m_tmp*n_tmp;i++) 
    {
      result[i] = tmp_var[i];
    }

  if (gradient!=NULL)
    {
      _SciErr = getVarAddressFromPosition(pvApiCtx, stack_pos+1, &tmp_addr); NLOPT_ERROR;
      _SciErr = getMatrixOfDouble(pvApiCtx, tmp_addr, &n_tmp, &m_tmp, &tmp_var); NLOPT_ERROR;
      // YC: check gradient size
      for(i=0;i<size_x*size_constr;i++) gradient[i] = tmp_var[i];
    }

  Nbvars = nbvars_old;
}

// This function is a wrapper which allows to call a Scilab script function like a "normal" 
// C/C++/Fortran function
void sci_constr_eq(unsigned int size_constr, double * result, unsigned int size_x, const double * x, double * gradient, void * param)
{
  int n_x = size_x, m_x = 1;
  int n_tmp, m_tmp, * tmp_addr = NULL;
  int sci_obj, lhs_obj, rhs_obj;
  int stack_pos, i, output_state;
  int rhs_old = Rhs, nbvars_old = Nbvars, Index = 0;
  double * x_sci = NULL;
  double * tmp_var = NULL;
  double tmp_val, f_out;
  SciErr _SciErr;

  sci_obj   = param_ceq.sci_obj;
  lhs_obj   = param_ceq.lhs_obj;
  rhs_obj   = param_ceq.rhs_obj;
  stack_pos = param_ceq.stack_pos;

  Nbvars = stack_pos + MAX(rhs_obj,lhs_obj);

  _SciErr = allocMatrixOfDouble(pvApiCtx, stack_pos+0, n_x, m_x, &x_sci); NLOPT_ERROR;
  for(i=0;i<size_x;i++) x_sci[i] = x[i];

  _SciErr = allocMatrixOfDouble(pvApiCtx, stack_pos+1, 1, 1, &tmp_var);

  SciFunction(&stack_pos, &sci_obj, &lhs_obj, &rhs_obj);

  if (Err>0) 
    {
      Scierror(999,"nlopt: error when calling equality constraint function\n");
      Nbvars = nbvars_old;
      return;
    } 

  _SciErr = getVarAddressFromPosition(pvApiCtx, stack_pos+0, &tmp_addr); NLOPT_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, tmp_addr, &m_tmp, &n_tmp, &tmp_var); NLOPT_ERROR;

  for(i=0;i<m_tmp*n_tmp;i++) result[i] = tmp_var[i];

  if (gradient!=NULL)
    {
      _SciErr = getVarAddressFromPosition(pvApiCtx, stack_pos+1, &tmp_addr); NLOPT_ERROR;
      _SciErr = getMatrixOfDouble(pvApiCtx, tmp_addr, &n_tmp, &m_tmp, &tmp_var); NLOPT_ERROR;
      // YC: check gradient size
      for(i=0;i<size_x*size_constr;i++) gradient[i] = tmp_var[i];
    }

  Nbvars = nbvars_old;
}

// This function returns a pointer to a function or a Scilab pointer to a function
// The choice is performed on the type of the given Scilab parameter
voidf GetFunctionPtr(char * name, int n, FTAB Table[], voidf scifun, int * ifunc, int * lhs, int * rhs) 
{
  int type, rep, m_tmp, n_tmp, i;
  int * tmp_addr = NULL;
  char * tmp_char = NULL;
  int * pi_len = NULL;
  char ** pst_strings = NULL;
  voidf f;
  SciErr _SciErr;

  _SciErr = getVarAddressFromPosition(pvApiCtx, n, &tmp_addr);
  _SciErr = getVarType(pvApiCtx, tmp_addr, &type);

  switch(type) 
    {
    case sci_strings: 
      _SciErr = getMatrixOfString(pvApiCtx, tmp_addr, &n_tmp, &m_tmp, NULL, NULL);
      pi_len = (int *)MALLOC(n_tmp*m_tmp*sizeof(int));
      pst_strings = (char **)MALLOC(n_tmp*m_tmp*sizeof(char *));
      _SciErr = getMatrixOfString(pvApiCtx, tmp_addr, &n_tmp, &m_tmp, pi_len, NULL);
      for(i=0;i<n_tmp*m_tmp;i++) pst_strings[i] = (char *)MALLOC((pi_len[i]+1)*sizeof(char));
      _SciErr = getMatrixOfString(pvApiCtx, tmp_addr, &n_tmp, &m_tmp, pi_len, pst_strings);

      f = SetFunction(pst_strings[0], &rep, Table);

      if (pst_strings) freeArrayOfString(pst_strings, n_tmp*m_tmp);
      if (pi_len) FREE(pi_len);

      if (rep==1)
	{
	  Scierror(999,"Function not found is %s\n", name);
	  return (voidf)0;
	}

      return f ;

    case sci_c_function: 
      GetRhsVar(n, "f", lhs, rhs, ifunc);
      return (voidf)scifun ;
      
    default: 
      Scierror(999,"Wrong parameter in %s ! (number %d)\r\n",name,n);
      return (voidf)0;
    }
}

// This function searches in the FTAB or in Scilab for corresponding function
voidf SetFunction(char * name, int * rep, FTAB table[]) 
{
  voidf loc;
  char * s = NULL;
  char buf[csiz];

  strncpy(buf,name,csiz);
  s = buf;
  while((*s!=' ')&&(*s != '\0')) {s++;};
  *s =  '\0';

  if (SearchComp(table,buf,&loc)==OK) 
    {
      *rep = 0;
      return(loc);
    }

  if (SearchInDynLinks(buf,&loc)>=0)
    {
      *rep = 0;
      return(loc);
    }

  loc = Emptyfunc;
  *rep = 1;

  return(loc);
}

// This function search in FTAB (here, we will use FTab_cobyla_function) for the corresponding name of the function
int SearchComp(FTAB Ftab[], char * op, void (**realop)()) 
{
  int i=0;

  while(Ftab[i].name!=(char *)0) 
    {
      int j;

      j = strcmp(op,Ftab[i].name);
      if ( j == 0 )
	{
	  *realop = Ftab[i].f;
	  return(OK);
	}
      else
	{ 
	  if ( j <= 0)
	    {
	      return(FAIL);
	    }
	  else i++;
	}
    }

  return(FAIL);
}
