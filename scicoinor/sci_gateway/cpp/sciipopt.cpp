//////////////////////////////////////////////////////////////////////////////////////////
// sciipopt: a scilab interface to the ipopt non linear constrained optimization solver //
//////////////////////////////////////////////////////////////////////////////////////////

//  Copyright (C) 2008-2010 Yann Collette.
//
//  SCIIPOPT is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2, or (at your option) any
//  later version.
//
//  This part of code is distributed with the FURTHER condition that it 
//  can be compiled and linked with the Matlab libraries and it can be 
//  used within the Matlab environment.
//
//  SCIIPOPT is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with SciCoinOR; see the file COPYING.  If not, write to the Free
//  Software Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <IpIpoptApplication.hpp>
#include <IpSolveStatistics.hpp>
#include <IpTNLP.hpp>
#include <IpIpoptCalculatedQuantities.hpp>
#include <IpSmartPtr.hpp>

#include <scilabjournal.hpp>
#include <call_function.hpp>

#include <iomanip>
#include <fstream>
#include <string>
#include <new>
#include <exception>
#include <string>
#include <algorithm>

#undef min
#undef max

extern "C"
{
#include <stack-c.h>
#include <sciprint.h>
#include <Scierror.h>
#include <MALLOC.h>
#include <api_parameters.h>
}

#include <api_scilab.h>

#include <def_solvers.hpp>
#include <manage_params.hpp>

// Prototypes for objective functions, gradient and constraints ...

// int objective(double * fobj, double * x, int n_size_x, double x_new);
// int objective_grad(double * fobj, double * x, int n_size_x, double x_new);
// int constraints(double * gobj, int n_size_constr, double * x, int n_size_x, double x_new);
// int constraint_jac(double * x, int n_size_x, double new_x, int n_size_g, int nnz_jac, 
//                    int * iRow, int * jCol, double * values, void * param);
// int hessian(double * x, int n_size_x, double new_x, double obj_factor, int n_size_g, double * lambda,
//             double new_lambda, int nele_hess, int * iRow, int * jCol, double * values, void * param);

typedef int (*call_f_objective)(double *, double *, int, double, void *);
typedef int (*call_f_objective_grad)(double *, double *, int, double, void *);
typedef int (*call_f_constraints)(double *, int, double *, int, double, void *);
typedef int (*call_f_constraints_jac)(double *, int, double, int, int, int *, int *, double *, void *);
typedef int (*call_f_hessian)(double *, int, double, double, int, double *,
			      double, int, int *, int *, double *, void *);

static FTAB FTab_ipopt_call_f[] = {{(char *) 0, (voidf) 0}};

int sci_ipopt_objective(double * x, double * f, int n_size_x, double x_new, void * param);
int sci_ipopt_objective_grad(double * x, double * f, int n_size_x, double x_new, void * param);
int sci_ipopt_constraints(double * x, int n_size_x, double * g, int n_size_g, double x_new, void * param);
int sci_ipopt_constraints_jac(double * x, int n_size_x, double new_x, int n_size_g, int nnz_jac, 
			      int * iRow, int * jCol, double * values, void * param);
int sci_ipopt_hessian(double * x, int n_size_x, double new_x, double obj_factor, int n_size_g, double * lambda,
		      double new_lambda, int nele_hess, int * iRow, int * jCol, double * values, void * param);

//#define DEBUG 1
//#define DBGPRINTF printf
#define DBGPRINTF sciprint

#include <string.h>

#include <helper.hpp>

using namespace Ipopt;

#define X_POS                   1
#define FOBJ_POS                2
#define DFOBJ_POS               3
#define GOBJ_POS                4
#define DGOBJ_POS               5
#define SPARSE_DGOBJ_POS        6
#define DHOBJ_POS               7
#define SPARSE_DHOBJ_POS        8
#define VAR_LIN_TYPE_POS        9
#define CONSTR_LIN_TYPE_POS     10
#define CONSTR_RHS_POS          11
#define CONSTR_LHS_POS          12
#define X_LOWER_POS             13
#define X_UPPER_POS             14
#define INTCB_POS               15
#define PARAMS_POS              16
#define LAST_PARAMS             PARAMS_POS
#define X_SOL_OUT_POS           Rhs+1
#define F_SOL_OUT_POS           Rhs+2
#define EXTRA_OUT_POS           Rhs+3

// A structure which handles informations related to objective functions
struct sci_ipopt_info
{
  int fobj_lhs,  fobj_rhs,  l_fobj,  fobj_is_list,  * fobj_list_addr;
  int dfobj_lhs, dfobj_rhs, l_dfobj, dfobj_is_list, * dfobj_list_addr;
  int gobj_lhs,  gobj_rhs,  l_gobj,  gobj_is_list,  * gobj_list_addr;
  int dgobj_lhs, dgobj_rhs, l_dgobj, dgobj_is_list, * dgobj_list_addr;
  int dhobj_lhs, dhobj_rhs, l_dhobj, dhobj_is_list, * dhobj_list_addr;
  int intcb_lhs, intcb_rhs, l_intcb;
  int m_sparse_dgobj, n_sparse_dgobj; // the sparsity structure of the jacobian
  int m_sparse_dhobj, n_sparse_dhobj; // the sparsity structure of the Hessian
  double * sparse_dgobj, * sparse_dhobj;
  call_f_objective       objective;
  call_f_objective_grad  objective_grad;
  call_f_constraints     constraints;
  call_f_constraints_jac constraints_jac;
  call_f_hessian         hessian;
  std::vector<param_fobj> fobj_parameter_list;  // List of parameters for the function
  std::vector<param_fobj> dfobj_parameter_list; // List of parameters for the function
  std::vector<param_fobj> gobj_parameter_list;  // List of parameters for the function
  std::vector<param_fobj> dgobj_parameter_list; // List of parameters for the function
  std::vector<param_fobj> dhobj_parameter_list; // List of parameters for the function
  int error;       // ==1 if an error occured
  int ibegin;      // the position of the top of the stack
  int nbvar;       // number of variables
  int nbconstr;    // number of constraints
  int nnz_jac_g;   // number of non zeroes in Jacobian
  int nnz_h_lag;   // number of non zeroes in Hessian of Lagrangean
  int index_style; // if the index starts from 0 (C) or 1 (FORTRAN)
  double * var_lin_type;    // position on the stack of the type of variables linearity (int)
  double * constr_lin_type; // position on the stack of the type of constraints linearity (int)
  double * constr_rhs;      // position on the stack of the upper bound for constraints (double)
  double * constr_lhs;      // position on the stack of the lower bound for constraints (double)
  double * x_lower;         // position on the stack of the lower bound for the variables (double)
  double * x_upper;         // position on the stack of the upper bound for the variables (double)
  double * x_0;             // position on the stack of the starting point (double)
  SmartPtr<IpoptApplication> ipopt_app; // a pointer to the current solver to retrieve the final statistics
  bool HessianPresent;
  double * pdbl_x_sol;
  double * pdbl_lambda;
  double f_sol;
  double status;
  double it_count;
  double cpu_time;
  double fobj_eval;
  double fobj_grad_eval;
  double constr_eval;
  double constr_jac_eval;
  double hess_eval;
  double dual_inf;
  double constr_viol;
  double complementarity;
  double kkt_error;
};


class IpoptTNLP : public TNLP
{
public:
  IpoptTNLP() : printSol_(false), sci_parameters(NULL) {}
  virtual ~IpoptTNLP() {}
  IpoptTNLP(const IpoptTNLP &other) : printSol_(other.printSol_), sci_parameters(other.sci_parameters) {}
  virtual bool get_variables_linearity(Index n, TNLP::LinearityType* var_types);
  virtual bool get_constraints_linearity(Index m, TNLP::LinearityType* const_types);
  virtual bool get_nlp_info(Index& n, Index&m, Index& nnz_jac_g,Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style);
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,Index m, Number* g_l, Number* g_u);
  virtual bool get_starting_point(Index n, bool init_x, Number* x,bool init_z, Number* z_L, Number* z_U,Index m, bool init_lambda,Number* lambda);
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,Index m, Index nele_jac, Index* iRow, Index *jCol,Number* values);
  virtual bool eval_h(Index n, const Number* x, bool new_x,Number obj_factor, Index m, const Number* lambda,bool new_lambda, 
		      Index nele_hess, Index* iRow, Index* jCol, Number* values);
  virtual Index get_number_of_nonlinear_variables();
  virtual bool  get_list_of_nonlinear_variables(Index num_nonlin_vars, Index* pos_nonlin_vars);

  virtual void finalize_solution(SolverReturn status,Index n, const Number *x, const Number *z_L, const Number *z_U, Index m, const Number *g, 
				 const Number *lambda, Number obj_value, const IpoptData *ip_data, IpoptCalculatedQuantities *ip_cq);
  virtual bool intermediate_callback(AlgorithmMode mode, Index iter, Number obj_value, Number inf_pr, Number inf_du, Number mu, 
				     Number d_norm, Number regularization_size, Number alpha_du, Number alpha_pr, Index ls_trials, 
				     const IpoptData *ip_data, IpoptCalculatedQuantities *ip_cq);
  void set_scilab_parameters(struct sci_ipopt_info * param) {sci_parameters = param;}
  struct sci_ipopt_info * get_scilab_parameters() {return sci_parameters;}
private:
  bool printSol_;
  struct sci_ipopt_info * sci_parameters;
};


bool IpoptTNLP::get_variables_linearity(Index n, TNLP::LinearityType* var_lin_types)
{
#ifdef DEBUG
  DBGPRINTF("sciipopt: in get_variables_linearity\n");
#endif
  int i;

  printf("in get_variables_linearity - n = %d\n", n);

  if (sci_parameters->var_lin_type)
    {
      for(i=0;i<n; i++) 
        {
          sciprint("var_lin_type[%d] = %d\n", i, (int)*(sci_parameters->var_lin_type+i));
          switch((int)*(sci_parameters->var_lin_type+i))
            {
            case 0:
              var_lin_types[i] = TNLP::LINEAR;
              break;
            default:
              var_lin_types[i] = TNLP::NON_LINEAR;
              break;
            }
        }
      return true;
    }
  else
    {
      return false;
    }
}

Index IpoptTNLP::get_number_of_nonlinear_variables()
{
#ifdef DEBUG
  DBGPRINTF("sciipopt: in get_number_of_nonlinear_variables\n");
#endif
  int i, _Index = 0;

  if (sci_parameters->var_lin_type)
    {
      for(i=0;i<sci_parameters->nbvar; i++) 
        {
          if ((int)*(sci_parameters->var_lin_type+i)!=0) _Index++;
        }
      return _Index;
    }
  else
    {
      return (Index)sci_parameters->nbvar;
    }
}

bool  IpoptTNLP::get_list_of_nonlinear_variables(Index num_nonlin_vars, Index* pos_nonlin_vars)
{
#ifdef DEBUG
  DBGPRINTF("sciipopt: in get_list_of_nonlinear_variables\n");
#endif
  int i, _Index = 0;
  
  if (sci_parameters->var_lin_type)
    {
      for(i=0;i<sci_parameters->nbvar; i++) 
        {
          if ((int)*(sci_parameters->var_lin_type+i)!=0) 
            {
              pos_nonlin_vars[_Index] = i + (int)(sci_parameters->index_style==0);
              _Index++;
            }
        }
      return true;
    }
  else
    {
      return false;
    }
}

bool IpoptTNLP::get_constraints_linearity(Index m, TNLP::LinearityType* const_lin_types)
{
#ifdef DEBUG
  DBGPRINTF("sciipopt: in get_constraints_linearity\n");
#endif
  int i;

  printf("in get_constraints_linearity - m = %d\n", m);

  if (sci_parameters->constr_lin_type)
    {
      for(i=0;i<m;i++) 
        {
          switch((int)*(sci_parameters->constr_lin_type+i))
            {
            case 0:
              const_lin_types[i] = TNLP::LINEAR;
              break;
            default:
              const_lin_types[i] = TNLP::NON_LINEAR;
              break;
            }
        }
      return true;
    }
  else
    {
      return false;
    }
}

 bool IpoptTNLP::get_nlp_info(Index& n, Index&m, Index& nnz_jac_g, Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
{
#ifdef DEBUG
  DBGPRINTF("sciipopt: in get_nlp_info\n");
#endif

  n = sci_parameters->nbvar;    // number of variable
  m = sci_parameters->nbconstr; // number of constraints
  nnz_jac_g = sci_parameters->nnz_jac_g; // number of non zeroes in Jacobian
  nnz_h_lag = sci_parameters->nnz_h_lag; // number of non zeroes in Hessian of Lagrangean
  index_style = (sci_parameters->index_style==0) ? TNLP::C_STYLE : TNLP::FORTRAN_STYLE;

  return true;
}

bool IpoptTNLP::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
{
#ifdef DEBUG
  DBGPRINTF("sciipopt: in get_bounds_info\n");
#endif

  int i;
  for(i=0;i<n;i++)
    {
      x_l[i] = *(sci_parameters->x_lower+i);
      x_u[i] = *(sci_parameters->x_upper+i);
    }

  for(i=0;i<m;i++)
    {
      g_l[i] = *(sci_parameters->constr_lhs+i);
      g_u[i] = *(sci_parameters->constr_rhs+i);
    }

  return true;
}

bool IpoptTNLP::get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda)
{
#ifdef DEBUG
  DBGPRINTF("sciipopt: in get_starting_point\n");
#endif

  int i;
  if (sci_parameters->x_0)
    {
      for(i=0;i<n;i++) x[i] = *(sci_parameters->x_0+i);
      
      return true;
    }
  else
    {
      return false;
    }
}

bool IpoptTNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
#ifdef DEBUG
  DBGPRINTF("sciipopt: in eval_f\n");
#endif

  if (sci_parameters->objective)
    {
      return (bool)(*sci_parameters->objective)((double *)x, (double *)&obj_value, (int)n, (double)new_x, (void *)sci_parameters);
    }
  else
    {
      return false;
    }
}

bool IpoptTNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
#ifdef DEBUG
  DBGPRINTF("sciipopt: in eval_grad_f\n");
#endif

  if (sci_parameters->objective_grad)
    {
      return (bool)(*sci_parameters->objective_grad)((double *)x, (double *)grad_f, (int)n, (double)new_x, (void *)sci_parameters);
    }
  else
    {
      return false;
    }
}

bool IpoptTNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
#ifdef DEBUG
  DBGPRINTF("sciipopt: in eval_g\n");
#endif

  if (sci_parameters->constraints)
    {
      return (bool)(*sci_parameters->constraints)((double *) x, (int)n, (double *)g, (int)m, (double)new_x, (void *)sci_parameters);
    }
  else
    {
      return false;
    }
}

bool IpoptTNLP::eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nnz_jac, Index* iRow, Index *jCol, Number* values)
{
#ifdef DEBUG
  DBGPRINTF("sciipopt: in eval_jac_g\n");
#endif

  if (sci_parameters->constraints_jac)
    {
      return (bool)(*sci_parameters->constraints_jac)((double *) x, (int)n, (double)new_x, (int)m, 
						      (int)nnz_jac, (int *)iRow, (int *)jCol, (double *)values, (void *)sci_parameters);
    }
  else
    {
      return false;
    }
}

bool IpoptTNLP::eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda,
		       bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values)
{
#ifdef DEBUG
  DBGPRINTF("sciipopt: in eval_h\n");
#endif

  if (sci_parameters->hessian)
    {
      return (int)(*sci_parameters->hessian)((double *)x, (int)n, (double)new_x, (double)obj_factor, (int)m, (double *)lambda,
					     (double)new_lambda, (int)nele_hess, (int *)iRow, (int *)jCol, (double *)values,
					     (void *)sci_parameters);
    }
  else
    {
      return false;
    }
}

void IpoptTNLP::finalize_solution(SolverReturn status,Index n, const Number *x, 
				  const Number *z_L, const Number *z_U, Index m, const Number *g, 
				  const Number *lambda, Number obj_value, 
				  const IpoptData *ip_data, IpoptCalculatedQuantities *ip_cq)
{
#ifdef DEBUG
  DBGPRINTF("sciipopt: in finalize_solution\n");
#endif

  int int_fobj_eval, int_constr_eval, int_fobj_grad_eval, int_constr_jac_eval, int_hess_eval;

  if (sci_parameters->error==1) status = (SolverReturn)Unrecoverable_Exception;

  //if (status == Solve_Succeeded)
  if (status>=0)
    {
      memcpy(sci_parameters->pdbl_x_sol,x,n*sizeof(double));
      sci_parameters->f_sol = obj_value;
    }
  else
    {
      memcpy(sci_parameters->pdbl_x_sol,sci_parameters->x_0,n*sizeof(double));
      sci_parameters->f_sol = 0;
    }

#ifdef DEBUG
  DBGPRINTF("sciipopt: finalize_solution - allocate status\n");
#endif

  // Value of status:

  // Solve_Succeeded                    0
  // Solved_To_Acceptable_Level         1
  // Infeasible_Problem_Detected        2
  // Search_Direction_Becomes_Too_Small 3
  // Diverging_Iterates                 4
  // User_Requested_Stop                5
  // Feasible_Point_Found               6
  // Maximum_Iterations_Exceeded        -1
  // Restoration_Failed                 -2
  // Error_In_Step_Computation          -3
  // Not_Enough_Degrees_Of_Freedom      -10
  // Invalid_Problem_Definition         -11
  // Invalid_Option                     -12
  // Invalid_Number_Detected            -13
  // Unrecoverable_Exception            -100
  // NonIpopt_Exception_Thrown          -101
  // Insufficient_Memory                -102
  // Internal_Error                     -199 

  sci_parameters->status = (double)status;

  // if (status == Solve_Succeeded)
  if (status>=0)
    {
      memcpy(sci_parameters->pdbl_lambda,lambda,m*sizeof(double));

      sci_parameters->it_count = (double)sci_parameters->ipopt_app->Statistics()->IterationCount();
      sci_parameters->cpu_time = sci_parameters->ipopt_app->Statistics()->TotalCPUTime();

      sci_parameters->ipopt_app->Statistics()->NumberOfEvaluations(int_fobj_eval, 
								   int_constr_eval, 
								   int_fobj_grad_eval,
								   int_constr_jac_eval, 
								   int_hess_eval);

      sci_parameters->fobj_eval       = (double)int_fobj_eval;
      sci_parameters->constr_eval     = (double)int_constr_eval;
      sci_parameters->fobj_grad_eval  = (double)int_fobj_grad_eval;
      sci_parameters->constr_jac_eval = (double)int_constr_jac_eval;
      sci_parameters->hess_eval       = (double)int_hess_eval;
  

      sci_parameters->ipopt_app->Statistics()->Infeasibilities(sci_parameters->dual_inf, 
							       sci_parameters->constr_viol, 
							       sci_parameters->complementarity,
							       sci_parameters->kkt_error);
    }
  else
    {
      memset(sci_parameters->pdbl_lambda,0,m*sizeof(double));

      sci_parameters->it_count        = -1;
      sci_parameters->fobj_eval       = -1;
      sci_parameters->constr_eval     = -1;
      sci_parameters->fobj_grad_eval  = -1;
      sci_parameters->constr_jac_eval = -1;
      sci_parameters->hess_eval       = -1;

      sci_parameters->cpu_time        = -1.0;
      sci_parameters->dual_inf        = -1.0;
      sci_parameters->constr_viol     = -1.0;
      sci_parameters->complementarity = -1.0;
      sci_parameters->kkt_error       = -1.0;
    }

#ifdef DEBUG
  DBGPRINTF("sciipopt: finalize_solution - create mlist\n");
  DBGPRINTF("it count = %f\n", sci_parameters->it_count);
  DBGPRINTF("cpu time = %f\n", sci_parameters->cpu_time);

  DBGPRINTF("fobj eval       = %f\n", sci_parameters->fobj_eval);
  DBGPRINTF("constr eval     = %f\n", sci_parameters->constr_eval);
  DBGPRINTF("fobj grad eval  = %f\n", sci_parameters->fobj_grad_eval);
  DBGPRINTF("constr jac eval = %f\n", sci_parameters->constr_jac_eval);
  DBGPRINTF("hess eval       = %f\n", sci_parameters->hess_eval);

  DBGPRINTF("dual inf        = %f\n", sci_parameters->dual_inf);
  DBGPRINTF("constr viol     = %f\n", sci_parameters->constr_viol);
  DBGPRINTF("complementarity = %f\n", sci_parameters->complementarity);
  DBGPRINTF("kkt error       = %f\n", sci_parameters->kkt_error);
#endif

#ifdef DEBUG
  DBGPRINTF("sciipopt: finalize_solution\n");
#endif
}

bool IpoptTNLP::intermediate_callback(AlgorithmMode mode, Index iter, Number obj_value, Number inf_pr, Number inf_du, Number mu, 
				      Number d_norm, Number regularization_size, Number alpha_du, Number alpha_pr, Index ls_trials, 
				      const IpoptData *ip_data, IpoptCalculatedQuantities *ip_cq)
{
  int * param_addr = NULL;
  int n_tmp, m_tmp, l_tmp, * tmp_addr = NULL;
  double * dbl_tmp = NULL;
  int * int_tmp = NULL;
  int nbvars_old = Nbvars;
  SciErr _SciErr;
  double * test_data = NULL;
  char * LabelList[11] = {"algorithm_mode", "iter", "obj_value", "inf_pr", "inf_du", "mu",
   			  "d_norm", "regularization_size", "alpha_du", "alpha_pr", "ls_trials"};

  if (sci_parameters->l_intcb==-1) return true;

#ifdef _MSC_VER
  Nbvars = sci_parameters->ibegin + max(sci_parameters->intcb_rhs, sci_parameters->intcb_lhs);
#else
  Nbvars = sci_parameters->ibegin + std::max(sci_parameters->intcb_rhs, sci_parameters->intcb_lhs);
#endif

  // Scilab prototype: [bool] = intermediate_callback(param)

  _SciErr = createPList(pvApiCtx, sci_parameters->ibegin+0, &param_addr, LabelList, 11); SCICOINOR_ERROR;
  
  _SciErr = createIntInPList(pvApiCtx,    sci_parameters->ibegin+0, param_addr, "algorithm_mode", mode); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx,    sci_parameters->ibegin+0, param_addr, "iter", iter); SCICOINOR_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, sci_parameters->ibegin+0, param_addr, "obj_value", obj_value); SCICOINOR_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, sci_parameters->ibegin+0, param_addr, "inf_pr", inf_pr); SCICOINOR_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, sci_parameters->ibegin+0, param_addr, "inf_du", inf_du); SCICOINOR_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, sci_parameters->ibegin+0, param_addr, "mu", mu); SCICOINOR_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, sci_parameters->ibegin+0, param_addr, "d_norm", d_norm); SCICOINOR_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, sci_parameters->ibegin+0, param_addr, "regularization_size", regularization_size); SCICOINOR_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, sci_parameters->ibegin+0, param_addr, "alpha_du", alpha_du); SCICOINOR_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, sci_parameters->ibegin+0, param_addr, "alpha_pr", alpha_pr); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx, sci_parameters->ibegin+0, param_addr, "ls_trials", ls_trials); SCICOINOR_ERROR;

  //////////////////////////
  // Infos from IpoptData //
  //////////////////////////

  // General algorithmic data

  // Index iter_count () const
  // bool MuInitialized () const
  // Number curr_tau () const
  // bool TauInitialized () const
  // bool FreeMuMode () const
  // bool tiny_step_flag ()
  // Number tol () const - Obtain the tolerance.

  // Information gathered for iteration output

  // Number info_regu_x () const
  // Number info_alpha_primal () const
  // Number info_alpha_dual () const
  // Index info_ls_count () const

  /////////////////////////////////////////
  // Info from IpoptCalculatedQuantities //
  /////////////////////////////////////////

  // SmartPtr< const Vector > curr_slack_x_L () - Slacks for x_L (at current iterate).
  // SmartPtr< const Vector > curr_slack_x_U () - Slacks for x_U (at current iterate).
  // SmartPtr< const Vector > curr_slack_s_L () - Slacks for s_L (at current iterate).
  // SmartPtr< const Vector > curr_slack_s_U () - Slacks for s_U (at current iterate).
  // SmartPtr< const Vector > trial_slack_x_L () - Slacks for x_L (at trial point).
  // SmartPtr< const Vector > trial_slack_x_U () - Slacks for x_U (at trial point).
  // SmartPtr< const Vector > trial_slack_s_L () - Slacks for s_L (at trial point).
  // SmartPtr< const Vector > trial_slack_s_U () - Slacks for s_U (at trial point).
  // Index AdjustedTrialSlacks () - Indicating whether or not we "fudged" the slacks.

  // Objective function

  // virtual Number curr_f () - Value of objective function (at current point).
  // virtual Number unscaled_curr_f () - Unscaled value of the objective function (at the current point).
  // virtual Number trial_f () - Value of objective function (at trial point).
  // virtual Number unscaled_trial_f () - Unscaled value of the objective function (at the trial point).
  // SmartPtr< const Vector > curr_grad_f () - Gradient of objective function (at current point).
  // SmartPtr< const Vector > trial_grad_f () - Gradient of objective function (at trial point).

  // Barrier Objective Function

  // virtual Number curr_barrier_obj () - Barrier Objective Function Value (at current iterate with current mu).
  // virtual Number trial_barrier_obj () - Barrier Objective Function Value (at trial point with current mu).
  // SmartPtr< const Vector > curr_grad_barrier_obj_x () - Gradient of barrier objective function with respect to x (at current point with current mu).
  // SmartPtr< const Vector > curr_grad_barrier_obj_s () - Gradient of barrier objective function with respect to s (at current point with current mu).
  // SmartPtr< const Vector > grad_kappa_times_damping_x () - Gradient of the damping term with respect to x (times kappa_d).
  // SmartPtr< const Vector > grad_kappa_times_damping_s () - Gradient of the damping term with respect to s (times kappa_d).

  // Constraints

  // SmartPtr< const Vector > curr_c () - c(x) (at current point)
  // SmartPtr< const Vector > unscaled_curr_c () - unscaled c(x) (at current point)
  // SmartPtr< const Vector > trial_c () - c(x) (at trial point)
  // SmartPtr< const Vector > curr_d () - d(x) (at current point)
  // SmartPtr< const Vector > unscaled_curr_d () - unscaled d(x) (at current point)
  // SmartPtr< const Vector > trial_d () - d(x) (at trial point)
  // SmartPtr< const Vector > curr_d_minus_s () - d(x) - s (at current point)
  // SmartPtr< const Vector > trial_d_minus_s () - d(x) - s (at trial point)
  // SmartPtr< const Matrix > curr_jac_c () - Jacobian of c (at current point).
  // SmartPtr< const Matrix > trial_jac_c () - Jacobian of c (at trial point).
  // SmartPtr< const Matrix > curr_jac_d () - Jacobian of d (at current point).
  // SmartPtr< const Matrix > trial_jac_d () - Jacobian of d (at trial point).
  // SmartPtr< const Vector > curr_jac_cT_times_vec (const Vector &vec) - Product of Jacobian (evaluated at current point) of C transpose with general vector.
  // SmartPtr< const Vector > trial_jac_cT_times_vec (const Vector &vec) - Product of Jacobian (evaluated at trial point) of C transpose with general vector.
  // SmartPtr< const Vector > curr_jac_dT_times_vec (const Vector &vec) - Product of Jacobian (evaluated at current point) of D transpose with general vector.
  // SmartPtr< const Vector > trial_jac_dT_times_vec (const Vector &vec) - Product of Jacobian (evaluated at trial point) of D transpose with general vector.
  // SmartPtr< const Vector > curr_jac_cT_times_curr_y_c () - Product of Jacobian (evaluated at current point) of C transpose with current y_c.
  // SmartPtr< const Vector > trial_jac_cT_times_trial_y_c () - Product of Jacobian (evaluated at trial point) of C transpose with trial y_c.
  // SmartPtr< const Vector > curr_jac_dT_times_curr_y_d () - Product of Jacobian (evaluated at current point) of D transpose with current y_d.
  // SmartPtr< const Vector > trial_jac_dT_times_trial_y_d () - Product of Jacobian (evaluated at trial point) of D transpose with trial y_d.
  // SmartPtr< const Vector > curr_jac_c_times_vec (const Vector &vec) - Product of Jacobian (evaluated at current point) of C with general vector.
  // SmartPtr< const Vector > curr_jac_d_times_vec (const Vector &vec) - Product of Jacobian (evaluated at current point) of D with general vector.
  // virtual Number curr_constraint_violation () - Constraint Violation (at current iterate).
  // virtual Number trial_constraint_violation () - Constraint Violation (at trial point).
  // virtual Number curr_nlp_constraint_violation (ENormType NormType) - Real constraint violation in a given norm (at current iterate).
  // virtual Number unscaled_curr_nlp_constraint_violation (ENormType NormType) - Unscaled real constraint violation in a given norm (at current iterate).

  // Call to the scilab function (intermediate_callback function)
  try
    {
      SciFunction(&sci_parameters->ibegin,&sci_parameters->l_intcb,&sci_parameters->intcb_lhs,&sci_parameters->intcb_rhs);
    }
  catch(...)
    {
      Scierror(999,"sciipopt: error when calling Intermediate callback function\n");
      Nbvars = nbvars_old;
      sci_parameters->error = 1;
      return false;
    }
  
  if (Err>0) 
    {
      Scierror(999,"sciipopt: error when calling Intermediate callback function\n");
      Nbvars = nbvars_old;
      sci_parameters->error = 1;
      return false;
    } /* End If */

  // We get the resulting vector
  _SciErr = getVarAddressFromPosition(pvApiCtx, sci_parameters->ibegin, &tmp_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, tmp_addr, &n_tmp, &m_tmp, &dbl_tmp); SCICOINOR_ERROR;

  Nbvars = nbvars_old;
  
  return bool(dbl_tmp[0]);
}

extern "C" int sciipopt(char * fname)
{
  // x0 fobj, dfobj, gobj, dgobj, params ...
  int m_x,               n_x,               * x_addr = NULL;
  int m_fobj,            n_fobj,            l_fobj,  fobj_is_list  = 0;
  int m_dfobj,           n_dfobj,           l_dfobj, dfobj_is_list = 0;
  int m_gobj,            n_gobj,            l_gobj,  gobj_is_list  = 0;
  int m_dgobj,           n_dgobj,           l_dgobj, dgobj_is_list = 0;
  int m_dhobj,           n_dhobj,           l_dhobj, dhobj_is_list = 0;
  int m_intcb,           n_intcb,           l_intcb, * intcb_addr = NULL;
  int m_sparse_dgobj,    n_sparse_dgobj,    * sparse_dgobj_addr = NULL; // sparsity structure of dgobj
  int m_sparse_dhobj,    n_sparse_dhobj,    * sparse_dhobj_addr = NULL; // sparsity structure of dgobj
  int m_var_lin_type,    n_var_lin_type,    * var_lin_type_addr = NULL;
  int m_constr_lin_type, n_constr_lin_type, * constr_lin_type_addr = NULL;
  int m_constr_rhs,      n_constr_rhs,      * constr_rhs_addr = NULL;
  int m_constr_lhs,      n_constr_lhs,      * constr_lhs_addr = NULL;
  int m_x_lower,         n_x_lower,         * x_lower_addr = NULL;
  int m_x_upper,         n_x_upper,         * x_upper_addr = NULL;
  double * x = NULL, * x_lower = NULL, * x_upper = NULL, * constr_rhs = NULL, * constr_lhs = NULL;
  double * var_lin_type = NULL, * constr_lin_type = NULL;
  double * sparse_dgobj = NULL, * sparse_dhobj = NULL;
  int Log = 0;
  bool nnz_jac_is_needed = false, nnz_hess_is_needed = false;
  int ipopt_status;
  int * param_in_addr = NULL;
  struct sci_ipopt_info * sci_parameters = new struct sci_ipopt_info;
  SciErr _SciErr;

  CheckRhs(LAST_PARAMS,LAST_PARAMS);
  CheckLhs(3,3);

  if (Rhs<LAST_PARAMS) 
    {
      Scierror(999,"%s: %d inputs required in call to %s. Bug in ipopt.sci ?...\n",fname, fname, LAST_PARAMS);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, X_POS, &x_addr); SCICOINOR_ERROR;
  if (~isEmptyMatrix(pvApiCtx, x_addr))
    {
      _SciErr = getMatrixOfDouble(pvApiCtx, x_addr, &m_x, &n_x, &x); SCICOINOR_ERROR;
    }
  else
    {
      x = NULL;
    }

  sci_parameters->error = 0;

  sci_parameters->fobj_parameter_list.resize(0);
  sci_parameters->dfobj_parameter_list.resize(0);
  sci_parameters->gobj_parameter_list.resize(0);
  sci_parameters->dgobj_parameter_list.resize(0);
  sci_parameters->dhobj_parameter_list.resize(0);

  sci_parameters->objective = 
    (call_f_objective)GetFunctionPtr("ipopt: objective", FOBJ_POS, FTab_ipopt_call_f, 
				     (voidf)sci_ipopt_objective, 
				     &l_fobj, &m_fobj, &n_fobj, &fobj_is_list,
				     sci_parameters->fobj_parameter_list);

  if (GetType(DFOBJ_POS)==sci_matrix)
    {
      m_dfobj = 0;
      n_dfobj = 0;
    }
  else
    {
      sci_parameters->objective_grad = 
	(call_f_objective_grad)GetFunctionPtr("ipopt: objective_grad", DFOBJ_POS, FTab_ipopt_call_f, 
					      (voidf)sci_ipopt_objective_grad, 
					      &l_dfobj, &m_dfobj, &n_dfobj, &dfobj_is_list,
					      sci_parameters->dfobj_parameter_list);
    }
  if (GetType(GOBJ_POS)==sci_matrix)
    {
      m_gobj = 0;
      n_gobj = 0;
      sci_parameters->constraints = NULL;
    }
  else
    {
      sci_parameters->constraints = 
	(call_f_constraints)GetFunctionPtr("ipopt: constraints", GOBJ_POS, FTab_ipopt_call_f, 
					   (voidf)sci_ipopt_constraints, 
					   &l_gobj, &m_gobj, &n_gobj, &gobj_is_list,
					   sci_parameters->gobj_parameter_list);
    }
  if (GetType(DGOBJ_POS)==sci_matrix)
    {
      m_dgobj = 0;
      n_dgobj = 0;
      sci_parameters->constraints_jac = NULL;
    }
  else
    {
      sci_parameters->constraints_jac = 
	(call_f_constraints_jac)GetFunctionPtr("ipopt: constraints_jac", DGOBJ_POS, FTab_ipopt_call_f, 
					       (voidf)sci_ipopt_constraints_jac, 
					       &l_dgobj, &m_dgobj, &n_dgobj, &dgobj_is_list,
					       sci_parameters->dgobj_parameter_list);

      // If we pass a C function here, we need the 'nnz_jac' option which
      // specifies the number of non zeros elements in the Jacobian
      if (sci_parameters->constraints_jac!=sci_ipopt_constraints_jac)
        {
          nnz_jac_is_needed = true;
        }
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, SPARSE_DGOBJ_POS, &sparse_dgobj_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, sparse_dgobj_addr, &m_sparse_dgobj, &n_sparse_dgobj, &sparse_dgobj); SCICOINOR_ERROR;

  if (GetType(DHOBJ_POS)==sci_matrix)
    {
      m_dhobj = 0;
      n_dhobj = 0;
      sci_parameters->hessian = NULL;
    }
  else
    {
      sci_parameters->hessian = 
	(call_f_hessian)GetFunctionPtr("ipopt: hessian", DHOBJ_POS, FTab_ipopt_call_f, 
				       (voidf)sci_ipopt_hessian, 
				       &l_dhobj, &m_dhobj, &n_dhobj, &dhobj_is_list,
				       sci_parameters->dhobj_parameter_list);
      if ((m_dhobj!=1)&&(n_dhobj!=3))
	{
	  Scierror(999,"%s: dhobj must take 3 input parameters and return 1 output parameter.\n",fname);
	  return 0;
	}

      // If we pass a C function here, we need the 'nnz_hess' option which
      // specifies the number of non zeros elements in the Hessian
      if (sci_parameters->hessian!=sci_ipopt_hessian)
        {
          nnz_hess_is_needed = true;
        }
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, INTCB_POS, &intcb_addr); SCICOINOR_ERROR;
  
  if (!isEmptyMatrix(pvApiCtx, intcb_addr))
    {
      GetRhsVar(INTCB_POS, EXTERNAL_DATATYPE, &m_intcb, &n_intcb, &l_intcb);
    }
  else
    {
      m_intcb = -1;
      n_intcb = -1;
      l_intcb = -1;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, SPARSE_DHOBJ_POS, &sparse_dhobj_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, sparse_dhobj_addr, &m_sparse_dhobj, &n_sparse_dhobj, &sparse_dhobj); SCICOINOR_ERROR;

  _SciErr = getVarAddressFromPosition(pvApiCtx, VAR_LIN_TYPE_POS, &var_lin_type_addr); SCICOINOR_ERROR;
  if (!isEmptyMatrix(pvApiCtx, var_lin_type_addr))
    {
      _SciErr = getMatrixOfDouble(pvApiCtx, var_lin_type_addr, &m_var_lin_type, &n_var_lin_type, &var_lin_type); SCICOINOR_ERROR;
    }
  else
    {
      var_lin_type = NULL;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, CONSTR_LIN_TYPE_POS, &constr_lin_type_addr); SCICOINOR_ERROR;
  if (!isEmptyMatrix(pvApiCtx, constr_lin_type_addr))
    {
      _SciErr = getMatrixOfDouble(pvApiCtx, constr_lin_type_addr, &m_constr_lin_type, &n_constr_lin_type, &constr_lin_type); SCICOINOR_ERROR;
    }
  else
    {
      constr_lin_type = NULL;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, CONSTR_RHS_POS, &constr_rhs_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, constr_rhs_addr, &m_constr_rhs, &n_constr_rhs, &constr_rhs); SCICOINOR_ERROR;

  _SciErr = getVarAddressFromPosition(pvApiCtx, CONSTR_LHS_POS, &constr_lhs_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, constr_lhs_addr, &m_constr_lhs, &n_constr_lhs, &constr_lhs); SCICOINOR_ERROR;
  
  _SciErr = getVarAddressFromPosition(pvApiCtx, X_LOWER_POS, &x_lower_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, x_lower_addr, &m_x_lower, &n_x_lower, &x_lower); SCICOINOR_ERROR;
  
  _SciErr = getVarAddressFromPosition(pvApiCtx, X_UPPER_POS, &x_upper_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, x_upper_addr, &m_x_upper, &n_x_upper, &x_upper); SCICOINOR_ERROR;

  if ((m_constr_rhs!=m_constr_lhs)&&(n_constr_rhs!=n_constr_lhs))
    {
      Scierror(999,"%s: constr_lhs and constr_rhs must have the same size.\n",fname);
      return 0;
    }

  if ((m_constr_rhs!=m_constr_lin_type)&&(n_constr_rhs!=n_constr_lin_type))
    {
      Scierror(999,"%s: constr_rhs and constr_lin_type must have the same size.\n",fname);
      return 0;
    }

  if ((m_x_upper!=m_x_lower)&&(n_x_upper!=n_x_lower))
    {
      Scierror(999,"%s: x_lower and x_upper must have the same size.\n",fname);
      return 0;
    }

  if ((m_x_upper!=m_var_lin_type)&&(n_x_upper!=n_var_lin_type))
    {
      Scierror(999,"%s: x_upper and var_lin_type must have the same size.\n",fname);
      return 0;
    }

  // Store the functions and the sparsity structure
  sci_parameters->fobj_lhs      = m_fobj;  
  sci_parameters->fobj_rhs      = n_fobj;  
  sci_parameters->l_fobj        = l_fobj;
  sci_parameters->fobj_is_list  = fobj_is_list;

  sci_parameters->dfobj_lhs     = m_dfobj;
  sci_parameters->dfobj_rhs     = n_dfobj;
  sci_parameters->l_dfobj       = l_dfobj;
  sci_parameters->dfobj_is_list = dfobj_is_list;

  sci_parameters->gobj_lhs      = m_gobj;
  sci_parameters->gobj_rhs      = n_gobj;
  sci_parameters->l_gobj        = l_gobj;
  sci_parameters->gobj_is_list  = gobj_is_list;

  sci_parameters->dgobj_lhs     = m_dgobj;
  sci_parameters->dgobj_rhs     = n_dgobj;
  sci_parameters->l_dgobj       = l_dgobj;
  sci_parameters->dgobj_is_list = dgobj_is_list;

  sci_parameters->dhobj_lhs     = m_dhobj;
  sci_parameters->dhobj_rhs     = n_dhobj;
  sci_parameters->l_dhobj       = l_dhobj;
  sci_parameters->dhobj_is_list = dhobj_is_list;

  sci_parameters->intcb_lhs     = m_intcb;
  sci_parameters->intcb_rhs     = n_intcb;
  sci_parameters->l_intcb       = l_intcb;

  if ((n_sparse_dgobj!=2)&&(n_sparse_dgobj!=0))
    {
      Scierror(999,"%s: sparse_dgobj must be a mx2 matrix or an empty value\n",fname);
      return 0;
    }

  if (sci_parameters->HessianPresent)
    {
      if ((n_sparse_dhobj!=2)&&(n_sparse_dhobj!=0))
	{
	  Scierror(999,"%s: sparse_dhobj must be a mx2 matrix or an empty value\n",fname);
	  return 0;
	}
    }

  // Store informations related to the sparsity structure of the Jacobian
  sci_parameters->sparse_dgobj   = sparse_dgobj;
  sci_parameters->n_sparse_dgobj = 2;
  sci_parameters->m_sparse_dgobj = m_sparse_dgobj;
  sci_parameters->nnz_jac_g      = sci_parameters->m_sparse_dgobj;

  // Store informations related to the sparsity structure of the Hessian
  sci_parameters->sparse_dhobj   = sparse_dhobj;
  sci_parameters->n_sparse_dhobj = 2;
  sci_parameters->m_sparse_dhobj = m_sparse_dhobj;

  if (GetType(DHOBJ_POS)!=sci_matrix)
    {
      sci_parameters->HessianPresent = true;
      sci_parameters->nnz_h_lag      = sci_parameters->m_sparse_dhobj;
    }
  else
    {
      sci_parameters->HessianPresent = false;
      sci_parameters->nnz_h_lag      = 0;
    }

  sci_parameters->ibegin            = Rhs + 1;
  sci_parameters->nbvar             = m_x*n_x;

  if (sci_parameters->constraints)
    {
      sci_parameters->nbconstr = m_constr_rhs * n_constr_rhs;
    }
  else
    sci_parameters->nbconstr        = 0;

  sci_parameters->index_style       = 1; // C_STYLE=0, FORTRAN_STYLE=1
  sci_parameters->var_lin_type      = var_lin_type;
  sci_parameters->constr_lin_type   = constr_lin_type;
  sci_parameters->constr_rhs        = constr_rhs;
  sci_parameters->constr_lhs        = constr_lhs;
  sci_parameters->x_lower           = x_lower;
  sci_parameters->x_upper           = x_upper;
  sci_parameters->x_0               = x;

  // Create the parameters for finalize_solution
  int m_x_sol           = 1,  n_x_sol           = 1;
  int m_f_sol           = 1,  n_f_sol           = 1;
  // Create the 'extra' structure of type plist
  int m_lambda          = 1,  n_lambda          = 1;
  int m_status          = 1,  n_status          = 1;
  int m_it_count        = 1,  n_it_count        = 1;
  int m_cpu_time        = 1,  n_cpu_time        = 1;
  int m_fobj_eval       = 1,  n_fobj_eval       = 1;
  int m_fobj_grad_eval  = 1,  n_fobj_grad_eval  = 1;
  int m_constr_eval     = 1,  n_constr_eval     = 1;
  int m_constr_jac_eval = 1,  n_constr_jac_eval = 1;
  int m_hess_eval       = 1,  n_hess_eval       = 1;
  int m_dual_inf        = 1,  n_dual_inf        = 1;
  int m_constr_viol     = 1,  n_constr_viol     = 1;
  int m_complementarity = 1,  n_complementarity = 1;
  int m_kkt_error       = 1,  n_kkt_error       = 1;
  int m_list_labels     = 1,  n_list_labels     = 14;

  m_x_sol  = sci_parameters->nbvar;
  n_lambda = sci_parameters->nbconstr;

  sci_parameters->pdbl_x_sol  = (double *)MALLOC(m_x_sol*sizeof(double));
  sci_parameters->pdbl_lambda = (double *)MALLOC(n_lambda*sizeof(double));
  sci_parameters->f_sol = 0;

  sci_parameters->status          = 0;
  sci_parameters->it_count        = 0;
  sci_parameters->cpu_time        = 0;
  sci_parameters->fobj_eval       = 0;
  sci_parameters->fobj_grad_eval  = 0;
  sci_parameters->constr_eval     = 0;
  sci_parameters->constr_jac_eval = 0;
  sci_parameters->hess_eval       = 0;
  sci_parameters->dual_inf        = 0;
  sci_parameters->constr_viol     = 0;
  sci_parameters->complementarity = 0;
  sci_parameters->kkt_error       = 0;

#ifdef DEBUG
  DBGPRINTF("sparse_dhobj   = %x\n", sci_parameters->sparse_dhobj);
  DBGPRINTF("m_sparse_dhobj = %d\n", sci_parameters->m_sparse_dhobj);
  DBGPRINTF("n_sparse_dhobj = %d\n", sci_parameters->n_sparse_dhobj);
  DBGPRINTF("nnz_h_lag      = %d\n", sci_parameters->nnz_h_lag);

  DBGPRINTF("sparse_dgobj   = %x\n", sci_parameters->sparse_dgobj);
  DBGPRINTF("m_sparse_dgobj = %d\n", sci_parameters->m_sparse_dgobj);
  DBGPRINTF("n_sparse_dgobj = %d\n", sci_parameters->n_sparse_dgobj);
  DBGPRINTF("nnz_jac_g      = %d\n", sci_parameters->nnz_jac_g);
#endif

  // Create an instance of your nlp...
  SmartPtr<IpoptTNLP> mynlp = new IpoptTNLP();

  // Create an instance of the IpoptApplication
  SmartPtr<IpoptApplication> ipopt_app = new IpoptApplication(false); // false because we don't want messages to be printed in the bash console

  sci_parameters->ipopt_app = ipopt_app;

  mynlp->set_scilab_parameters(sci_parameters);

  // Get the parameters stored in the plist
  int tmp_res, tmp_int;

  _SciErr = initPList(pvApiCtx, PARAMS_POS, &param_in_addr); SCICOINOR_ERROR;
  if (!checkPList(pvApiCtx, param_in_addr))
    {
      Scierror(999, "%s: argument n° %d is not a plist\n", fname, PARAMS_POS);
      return 0;
    }

#ifdef DEBUG
  DBGPRINTF("sciipopt: processing options\n");
#endif

  /////////////
  // Journal //
  /////////////

  // EJournalLevel:

  // J_INSUPPRESSIBLE -1
  // J_NONE 	       0
  // J_ERROR 	       1
  // J_STRONGWARNING   2
  // J_SUMMARY 	       3
  // J_WARNING 	       4
  // J_ITERSUMMARY     5
  // J_DETAILED        6
  // J_MOREDETAILED    7
  // J_VECTOR 	       8
  // J_MOREVECTOR      9
  // J_MATRIX 	       10
  // J_MOREMATRIX      11
  // J_ALL 	       12
  // J_LAST_LEVEL      13

  // Add a new journal to print all the ipopt informations in the scilab console
  _SciErr = getIntInPList(pvApiCtx, param_in_addr, "journal_level", &tmp_int, &tmp_res, 0, Log, CHECK_BOTH,-1,12); SCICOINOR_ERROR;
  SmartPtr<Journal> console = new ScilabJournal((EJournalLevel)tmp_int);
  ipopt_app->Jnlst()->AddJournal(console);

  if (nnz_jac_is_needed)
    {
      _SciErr = getIntInPList(pvApiCtx, param_in_addr, "nnz_jac", &tmp_int, &tmp_res, -1, Log, CHECK_NONE); SCICOINOR_ERROR;
      if (tmp_res==-1) 
        {
          Scierror(999,"%s: the Jacobian function is a C function, you need to specify the nnz_jac option.\n", fname);
          return 0;
        }
      sci_parameters->nnz_jac_g = tmp_int;
    }

  if (nnz_hess_is_needed)
    {
      _SciErr = getIntInPList(pvApiCtx, param_in_addr, "nnz_hess", &tmp_int, &tmp_res, -1, Log, CHECK_NONE); SCICOINOR_ERROR;
      if (tmp_res==-1) 
        {
          Scierror(999,"%s: the Hessian function is a C function, you need to specify the nnz_hess option.\n", fname);
          return 0;
        }
      sci_parameters->nnz_h_lag = tmp_int;
    }

  // Manage all the remaining parameters of Ipopt
  manage_ipopt_params(ipopt_app->Options(), param_in_addr, Log);

  /////////////////////////////////////////////////////////////
  // Initialize the IpoptApplication and process the options //
  /////////////////////////////////////////////////////////////

  ApplicationReturnStatus status;

  status = ipopt_app->Initialize();
  if (status != Solve_Succeeded) 
    {
      Scierror(999,"%s: Error during initialization!\n",fname);
      return 0;
    }

#ifdef DEBUG
  DBGPRINTF("sciipopt: launching optimization\n");
#endif

  // Value of status:
  // Solve_Succeeded                    = 0;
  // Solved_To_Acceptable_Level         = 1;
  // Infeasible_Problem_Detected        = 2;
  // Search_Direction_Becomes_Too_Small = 3;
  // Diverging_Iterates                 = 4;
  // User_Requested_Stop                = 5;
  // Feasible_Point_Found               = 6;
  // Maximum_Iterations_Exceeded        = -1;
  // Restoration_Failed                 = -2;
  // Error_In_Step_Computation          = -3;
  // Not_Enough_Degrees_Of_Freedom      = -10;
  // Invalid_Problem_Definition         = -11;
  // Invalid_Option                     = -12;
  // Invalid_Number_Detected            = -13;
  // Unrecoverable_Exception            = -100;
  // NonIpopt_Exception_Thrown          = -101;
  // Insufficient_Memory                = -102;
  // Internal_Error                     = -199 

  try 
    {
      ipopt_status = 0;
      ipopt_status = ipopt_app->OptimizeTNLP((SmartPtr<TNLP>&)mynlp);
    }
  catch(IpoptException &E) 
    {
      //There has been a failure to solve a problem with Ipopt.
      sciprint("sci_ipopt: Exception - %s\n",E.Message().c_str());
      if (sci_parameters) delete sci_parameters;
      return 0;
    }
  catch(std::bad_alloc &E) 
    {
      sciprint("sci_ipopt: allocation problem - %s\n", E.what());
      if (sci_parameters) delete sci_parameters;
      return 0;
    }

#ifdef DEBUG
  DBGPRINTF("sciipopt: exiting optimization: status = %d\n", status);
#endif

  int * extra_addr = NULL, found = 0;
  char * ListLabels [] = {"lambda","status","it_count","cpu_time",
                          "fobj_eval","fobj_grad_eval","constr_eval","constr_jac_eval","hess_eval",
                          "dual_inf","constr_viol","complementarity","kkt_error"};

  _SciErr = createMatrixOfDouble(pvApiCtx, X_SOL_OUT_POS, m_x_sol, n_x_sol, sci_parameters->pdbl_x_sol); SCICOINOR_ERROR;
  _SciErr = createMatrixOfDouble(pvApiCtx, F_SOL_OUT_POS, m_f_sol, n_f_sol, &sci_parameters->f_sol); SCICOINOR_ERROR;

  _SciErr = createPList(pvApiCtx, EXTRA_OUT_POS, &extra_addr, (char **)ListLabels, 13); SCICOINOR_ERROR;

  _SciErr = createColVectorOfDoubleInPList(pvApiCtx, EXTRA_OUT_POS, extra_addr, "lambda", m_lambda*n_lambda, sci_parameters->pdbl_lambda); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx, EXTRA_OUT_POS,    extra_addr, "status",           sci_parameters->status); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx, EXTRA_OUT_POS,    extra_addr, "it_count",         sci_parameters->it_count); SCICOINOR_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, EXTRA_OUT_POS, extra_addr, "cpu_time",         sci_parameters->cpu_time); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx, EXTRA_OUT_POS,    extra_addr, "fobj_eval",        sci_parameters->fobj_eval); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx, EXTRA_OUT_POS,    extra_addr, "fobj_grad_eval",   sci_parameters->fobj_grad_eval); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx, EXTRA_OUT_POS,    extra_addr, "constr_eval",      sci_parameters->constr_eval); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx, EXTRA_OUT_POS,    extra_addr, "constr_jac_eval",  sci_parameters->constr_jac_eval); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx, EXTRA_OUT_POS,    extra_addr, "hess_eval",        sci_parameters->hess_eval); SCICOINOR_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, EXTRA_OUT_POS, extra_addr, "dual_inf",         sci_parameters->dual_inf); SCICOINOR_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, EXTRA_OUT_POS, extra_addr, "constr_viol",      sci_parameters->constr_viol); SCICOINOR_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, EXTRA_OUT_POS, extra_addr, "complementarity",  sci_parameters->complementarity); SCICOINOR_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, EXTRA_OUT_POS, extra_addr, "kkt_error",        sci_parameters->kkt_error); SCICOINOR_ERROR;

  LhsVar(1) = X_SOL_OUT_POS;
  LhsVar(2) = F_SOL_OUT_POS;
  LhsVar(3) = EXTRA_OUT_POS;

  if (sci_parameters->pdbl_x_sol)  FREE(sci_parameters->pdbl_x_sol);
  if (sci_parameters->pdbl_lambda) FREE(sci_parameters->pdbl_lambda);
  if (sci_parameters) delete sci_parameters;

  return 0;
}

//
// Definition of the "fake" objective and constraints.
// These functions allow to be able to deal with scilab scripts or C functions
//

int sci_ipopt_objective(double * x, double * f, int n_size_x, double x_new, void * param)
{
  int n_tmp_1, m_tmp_1, * tmp_1_addr = NULL;
  int n_tmp_2, m_tmp_2;
  double * tmp_1 = NULL, * tmp_2 = NULL;
  int nbvars_old = Nbvars;
  int Index = 0;
  struct sci_ipopt_info * sci_parameters = (struct sci_ipopt_info *)param;
  SciErr _SciErr;

  ////////////////////////////////
  // Call to objective function //
  ////////////////////////////////

#ifdef _MSC_VER
  Nbvars = sci_parameters->ibegin + max(sci_parameters->fobj_rhs,sci_parameters->fobj_lhs);
#else
  Nbvars = sci_parameters->ibegin + std::max(sci_parameters->fobj_rhs,sci_parameters->fobj_lhs);
#endif

  // Store the current point into tmp_1: x
  n_tmp_1 = n_size_x;
  m_tmp_1 = 1;

  _SciErr = allocMatrixOfDouble(pvApiCtx, sci_parameters->ibegin, n_tmp_1, m_tmp_1, &tmp_1); SCICOINOR_ERROR;
  memcpy(tmp_1, x, sizeof(double)*n_size_x);

  if (sci_parameters->fobj_rhs>=2)
    {
      // Store x_new into tmp_2
      n_tmp_2 = 1;
      m_tmp_2 = 1;
      
      _SciErr = allocMatrixOfDouble(pvApiCtx, sci_parameters->ibegin+1, n_tmp_2, m_tmp_2, &tmp_2); SCICOINOR_ERROR;
      *tmp_2 = x_new;
    }

  // Call to the scilab function (objective function)
  try
    {
      SciFunction(&sci_parameters->ibegin,&sci_parameters->l_fobj,&sci_parameters->fobj_lhs,&sci_parameters->fobj_rhs);
    }
  catch(...)
    {
      Scierror(999,"sci_ipopt: error when calling objective function\n");
      sci_parameters->error = 1;
      Nbvars = nbvars_old;
      return 0;
    }

  if (Err>0) 
    {
      Scierror(999,"sci_ipopt: error when calling objective function\n");
      sci_parameters->error = 1;
      Nbvars = nbvars_old;
      return 0;
    } /* End If */
  
  // Get fobj
  _SciErr = getVarAddressFromPosition(pvApiCtx, sci_parameters->ibegin, &tmp_1_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, tmp_1_addr, &n_tmp_1, &m_tmp_1, &tmp_1); SCICOINOR_ERROR;

  *f = *tmp_1;

  Nbvars = nbvars_old;

  return (int)true;
}

int sci_ipopt_objective_grad(double * x, double * df, int n_size_x, double x_new, void * param)
{
  int n_tmp_1, m_tmp_1, * tmp_1_addr = NULL;
  int n_tmp_2, m_tmp_2;
  double * tmp_1 = NULL, * tmp_2 = NULL;
  int nbvars_old = Nbvars;
  struct sci_ipopt_info * sci_parameters = (struct sci_ipopt_info *)param;
  SciErr _SciErr;

  ////////////////////////////////////////////////
  // Call to the gradient of objective function //
  ////////////////////////////////////////////////

#ifdef _MSC_VER
  Nbvars = sci_parameters->ibegin + max(sci_parameters->dfobj_rhs,sci_parameters->dfobj_lhs);
#else
  Nbvars = sci_parameters->ibegin + std::max(sci_parameters->dfobj_rhs,sci_parameters->dfobj_lhs);
#endif

  // Store the current point into tmp_1: x
  n_tmp_1 = n_size_x;
  m_tmp_1 = 1;

  _SciErr = allocMatrixOfDouble(pvApiCtx, sci_parameters->ibegin, n_tmp_1, m_tmp_1, &tmp_1); SCICOINOR_ERROR;
  memcpy(tmp_1, x, sizeof(double)*n_size_x);

  if (sci_parameters->dfobj_rhs>=2)
    {
      // Store x_new into tmp_2
      n_tmp_2 = 1;
      m_tmp_2 = 1;
      
      _SciErr = allocMatrixOfDouble(pvApiCtx, sci_parameters->ibegin+1, n_tmp_2, m_tmp_2, &tmp_2); SCICOINOR_ERROR;
      *tmp_2 = x_new;
    }

  // Call to the scilab function (gradient of objective function)
  try
    {
      SciFunction(&sci_parameters->ibegin,&sci_parameters->l_dfobj,&sci_parameters->dfobj_lhs,&sci_parameters->dfobj_rhs);
    }
  catch(...)
    {
      Scierror(999,"sci_ipopt: error when calling gradient of objective function\n");
      sci_parameters->error = 1;
      Nbvars = nbvars_old;
      return 0;
    }

  if (Err>0) 
    {
      Scierror(999,"sci_ipopt: error when calling gradient of objective function\n");
      sci_parameters->error = 1;
      Nbvars = nbvars_old;
      return 0;
    } /* End If */

  if (GetType(sci_parameters->ibegin)==sci_sparse)
    {
      Scierror(999,"sci_ipopt: the gradient of the objective function must be not sparse\n");
      sci_parameters->error = 1;
      Nbvars = nbvars_old;
      return 0;
    }

  // Get dfobj
  _SciErr = getVarAddressFromPosition(pvApiCtx, sci_parameters->ibegin, &tmp_1_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, tmp_1_addr, &n_tmp_1, &m_tmp_1, &tmp_1); SCICOINOR_ERROR;

  memcpy(df,tmp_1,n_size_x*sizeof(double));

#ifdef DEBUG
  DBGPRINTF("sciipopt: leaving eval_grad_f\n");
#endif

  Nbvars = nbvars_old;

  return (int)true;
}

int sci_ipopt_constraints(double * x, int n_size_x, double * g, int n_size_g, double x_new, void * param)
{
#ifdef DEBUG
  DBGPRINTF("sciipopt: in eval_g\n");
#endif

  int n_tmp_1, m_tmp_1, * tmp_1_addr = NULL;
  int n_tmp_2, m_tmp_2, * tmp_2_addr = NULL;
  double * tmp_1 = NULL, * tmp_2 = NULL;
  int nbvars_old = Nbvars;
  struct sci_ipopt_info * sci_parameters = (struct sci_ipopt_info *)param;
  SciErr _SciErr;
  
  if (sci_parameters->nbconstr==0) return (int)false;

  /////////////////////////////////////
  // Call to the constraint function //
  /////////////////////////////////////

#ifdef _MSC_VER
  Nbvars = sci_parameters->ibegin + max(sci_parameters->gobj_rhs,sci_parameters->gobj_lhs);
#else
  Nbvars = sci_parameters->ibegin + std::max(sci_parameters->gobj_rhs,sci_parameters->gobj_lhs);
#endif

  // Store the current point into tmp_1: x
  n_tmp_1 = n_size_x;
  m_tmp_1 = 1;

  _SciErr = allocMatrixOfDouble(pvApiCtx, sci_parameters->ibegin, n_tmp_1, m_tmp_1, &tmp_1); SCICOINOR_ERROR;
  memcpy(tmp_1, x, sizeof(double)*n_size_x);

  if (sci_parameters->gobj_rhs>=2)
    {
      // Store new_x into tmp_2
      n_tmp_2 = 1;
      m_tmp_2 = 1;
      
      _SciErr = allocMatrixOfDouble(pvApiCtx, sci_parameters->ibegin+1, n_tmp_2, m_tmp_2, &tmp_2); SCICOINOR_ERROR;
      *tmp_2 = x_new;
    }

  // Call to the scilab function (constraint function)
  try
    {
      SciFunction(&sci_parameters->ibegin,&sci_parameters->l_gobj,&sci_parameters->gobj_lhs,&sci_parameters->gobj_rhs);
    }
  catch(...)
    {
      Scierror(999,"sci_ipopt: error when calling constraint function\n");
      sci_parameters->error = 1;
      Nbvars = nbvars_old;
      return 0;
    }

  if (Err>0) 
    {
      Scierror(999,"sci_ipopt: error when calling constraint function\n");
      sci_parameters->error = 1;
      Nbvars = nbvars_old;
      return 0;
    } /* End If */

  // Get gobj
  _SciErr = getVarAddressFromPosition(pvApiCtx, sci_parameters->ibegin, &tmp_1_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, tmp_1_addr, &n_tmp_1, &m_tmp_1, &tmp_1); SCICOINOR_ERROR;

  memcpy(g,tmp_1,n_size_g*sizeof(double));

  Nbvars = nbvars_old;

  return (int)true;
}

int sci_ipopt_constraints_jac(double * x, int n_size_x, double new_x, int n_size_g, int nnz_jac, 
			      int * iRow, int * jCol, double * values, void * param)
{
#ifdef DEBUG
  DBGPRINTF("sciipopt: in eval_jac_g\n");
#endif

  int n_tmp_1, m_tmp_1, * tmp_1_addr = NULL;
  int n_tmp_2, m_tmp_2, * tmp_2_addr = NULL;
  double * tmp_1 = NULL, * tmp_2 = NULL;
  int i, Index;
  int nbvars_old = Nbvars;
  struct sci_ipopt_info * sci_parameters = (struct sci_ipopt_info *)param;
  SciErr _SciErr;

  if (sci_parameters->nbconstr==0) return (int)false;

#ifdef _MSC_VER
  Nbvars = sci_parameters->ibegin + max(sci_parameters->dgobj_rhs,sci_parameters->dgobj_lhs);
#else
  Nbvars = sci_parameters->ibegin + std::max(sci_parameters->dgobj_rhs,sci_parameters->dgobj_lhs);
#endif

  // Fill iRow and jCol for the sparsity structure of the hessian
  if ((iRow!=NULL)&&(jCol!=NULL)) 
    {
      Index = 0;
      for(i=0; i<sci_parameters->m_sparse_dgobj; i++)
	{
	  iRow[Index] = (int)*(sci_parameters->sparse_dgobj + i + 0 * sci_parameters->m_sparse_dgobj);
	  jCol[Index] = (int)*(sci_parameters->sparse_dgobj + i + 1 * sci_parameters->m_sparse_dgobj);
	  Index++;
	}
    }

  // Fill the values of the jacobian
  // If values!=NULL, then the current point x has been shipped
  if (values!=NULL)
    {
#ifdef DEBUG
      DBGPRINTF("sciipopt: call to the jacobian\n");
#endif

      //////////////////////////////////////////////////
      // Call to the gradient of constraints function //
      //////////////////////////////////////////////////

      // Store the current point into tmp_1: x
      n_tmp_1 = n_size_x;
      m_tmp_1 = 1;
      
      _SciErr = allocMatrixOfDouble(pvApiCtx, sci_parameters->ibegin, n_tmp_1, m_tmp_1, &tmp_1); SCICOINOR_ERROR;

      memcpy(tmp_1, x, sizeof(double)*n_size_x);
  
      if (sci_parameters->dgobj_rhs>=2)
	{
	  // Store new_x into tmp_2
	  n_tmp_2 = 1;
	  m_tmp_2 = 1;
	  
	  _SciErr = allocMatrixOfDouble(pvApiCtx, sci_parameters->ibegin+1, n_tmp_2, m_tmp_2, &tmp_2); SCICOINOR_ERROR;
	  *tmp_2 = new_x;
	}

#ifdef DEBUG
      DBGPRINTF("sciipopt: call to the function\n");
#endif

      // Call to the scilab function (constraint function)
      try
	{
	  SciFunction(&sci_parameters->ibegin,&sci_parameters->l_dgobj,&sci_parameters->dgobj_lhs,&sci_parameters->dgobj_rhs);
	}
      catch(...)
	{
	  Scierror(999,"sci_ipopt: error when calling Jacobian of the constraints function\n");
	  sci_parameters->error = 1;
	  Nbvars = nbvars_old;
	  return 0;
	}

      if (Err>0) 
	{
	  Scierror(999,"sci_ipopt: error when calling Jacobian of the constraints function\n");
	  sci_parameters->error = 1;
	  return 0;
	} /* End If */
      
      if (GetType(sci_parameters->ibegin)==sci_sparse)
	{
	  Scierror(999,"sciipopt: the Jacobian of the constraints must be not sparse\n");
	  sci_parameters->error = 1;
	  Nbvars = nbvars_old;
	  return 0;
	}

      // Get dG
      _SciErr = getVarAddressFromPosition(pvApiCtx, sci_parameters->ibegin, &tmp_1_addr); SCICOINOR_ERROR;
      _SciErr = getMatrixOfDouble(pvApiCtx, tmp_1_addr, &n_tmp_1, &m_tmp_1, &tmp_1); SCICOINOR_ERROR;

      if (n_tmp_1*m_tmp_1!=nnz_jac)
	{
	  Scierror(999,"sciipopt: the Jacobian of the constraint function must return %d values\n",nnz_jac);
	  sci_parameters->error = 1;
	  Nbvars = nbvars_old;
	  return 0;
	}

      // We must return a vector of values. Do we check ??
      memcpy(values, tmp_1, nnz_jac*sizeof(double));
    }

  Nbvars = nbvars_old;

  return (int)true;
}

int sci_ipopt_hessian(double * x, int n_size_x, double new_x, double obj_factor, int n_size_g, double * lambda,
		      double new_lambda, int nele_hess, int * iRow, int * jCol, double * values, void * param)
{
#ifdef DEBUG
  DBGPRINTF("sciipopt: in eval_h\n");
#endif

  int n_tmp_1, m_tmp_1, * tmp_1_addr = NULL;
  int n_tmp_2, m_tmp_2;
  int n_tmp_3, m_tmp_3;
  int n_tmp_4, m_tmp_4;
  int n_tmp_5, m_tmp_5;
  double * tmp_1 = NULL, * tmp_2 = NULL, * tmp_3 = NULL;
  double * tmp_4 = NULL, * tmp_5 = NULL;
  int i, Index, var_type;
  int nbvars_old = Nbvars;
  struct sci_ipopt_info * sci_parameters = (struct sci_ipopt_info *)param;
  SciErr _SciErr;

  if (sci_parameters->HessianPresent)
    {
#ifdef _MSC_VER
      Nbvars = sci_parameters->ibegin + max(sci_parameters->dhobj_rhs,sci_parameters->dhobj_lhs);
#else
      Nbvars = sci_parameters->ibegin + std::max(sci_parameters->dhobj_rhs,sci_parameters->dhobj_lhs);
#endif

      // Fill iRow and jCol for the sparsity structure of the hessian
      if ((iRow!=NULL)&&(jCol!=NULL)) 
	{
#ifdef DEBUG
	  DBGPRINTF("sciipopt: in eval_h - sparsity\n");
#endif
	  Index = 0;
	  for(i=0; i<sci_parameters->m_sparse_dhobj; i++)
	    {
	      iRow[Index] = (int)*(sci_parameters->sparse_dhobj + i + 0 * sci_parameters->m_sparse_dhobj);
	      jCol[Index] = (int)*(sci_parameters->sparse_dhobj + i + 1 * sci_parameters->m_sparse_dhobj);
	      Index++;
	    }
	}
      
      // If values!=NULL, then the current point x has been shipped
      if (values!=NULL)
	{
#ifdef DEBUG
	  DBGPRINTF("sciipopt: in eval_h - hessian\n");
#endif
	  //////////////////////////////////
	  // Call to the Hessian function //
	  //////////////////////////////////
	  
	  // Store the current point into tmp_1
	  n_tmp_1 = n_size_x;
	  m_tmp_1 = 1;
	  
	  _SciErr = allocMatrixOfDouble(pvApiCtx, sci_parameters->ibegin, n_tmp_1, m_tmp_1, &tmp_1); SCICOINOR_ERROR;
	  
	  memcpy(tmp_1, x, sizeof(double)*n_size_x);

	  // Store the lambda vector into tmp_2
	  n_tmp_2 = n_size_g;
	  m_tmp_2 = 1;
	  
	  _SciErr = allocMatrixOfDouble(pvApiCtx, sci_parameters->ibegin+1, n_tmp_2, m_tmp_2, &tmp_2); SCICOINOR_ERROR;
	  
	  // If values!=NULL, then the current lambda vector has been shipped
	  memcpy(tmp_2, lambda, sizeof(double)*n_size_g);
	  
	  // Store the objective function weight into tmp_3
	  n_tmp_3 = 1;
	  m_tmp_3 = 1;
	  
	  _SciErr = allocMatrixOfDouble(pvApiCtx, sci_parameters->ibegin+2, n_tmp_3, m_tmp_3, &tmp_3); SCICOINOR_ERROR;
	  
	  // If values!=NULL, then the current objective function weight has been shipped
	  *tmp_3 = obj_factor;
	  
	  if (sci_parameters->dhobj_rhs>=4)
	    {
	      // Store new_x into tmp_4
	      n_tmp_4 = 1;
	      m_tmp_4 = 1;
	      
	      _SciErr = allocMatrixOfDouble(pvApiCtx, sci_parameters->ibegin+3, n_tmp_4, m_tmp_4, &tmp_4); SCICOINOR_ERROR;
	      *tmp_4 = new_x;
	    }

	  if (sci_parameters->dhobj_rhs>=5)
	    {
	      // Store new_lambda into tmp_5
	      n_tmp_5 = 1;
	      m_tmp_5 = 1;
	      
	      _SciErr = allocMatrixOfDouble(pvApiCtx, sci_parameters->ibegin+4, n_tmp_5, m_tmp_5, &tmp_5); SCICOINOR_ERROR;
	      *tmp_5 = new_lambda;
	    }

#ifdef DEBUG
	  DBGPRINTF("sciipopt: in eval_h - call to the function\n");
	  DBGPRINTF("sciipopt: in eval_h - Hessian = %d\n",sci_parameters->HessianPresent);
#endif

	  // Call to the scilab function (Hessian function)
	  try
	    {
	      SciFunction(&sci_parameters->ibegin,&sci_parameters->l_dhobj,&sci_parameters->dhobj_lhs,&sci_parameters->dhobj_rhs);
	    }
	  catch(...)
	    {
	      Scierror(999,"sciipopt: error when calling Hessian function\n");
	      sci_parameters->error = 1;
	      Nbvars = nbvars_old;
	      return 0;
	    }

	  if (Err>0) 
	    {
	      Scierror(999,"sciipopt: error when calling Hessian function\n");
	      sci_parameters->error = 1;
	      Nbvars = nbvars_old;
	      return 0;
	    } /* End If */
	  
	  _SciErr = getVarAddressFromPosition(pvApiCtx, sci_parameters->ibegin, &tmp_1_addr); SCICOINOR_ERROR;
	  _SciErr = getVarType(pvApiCtx, tmp_1_addr, &var_type); SCICOINOR_ERROR;
	  if (var_type==sci_sparse)
	    {
	      Scierror(999,"sciipopt: the Hessian function must be not sparse\n");
	      Nbvars = nbvars_old;
	      return 0;
	    }
	  
#ifdef DEBUG
	  DBGPRINTF("sciipopt: in eval_h - function called\n");
#endif
	  
	  // We get the resulting vector
	  _SciErr = getMatrixOfDouble(pvApiCtx, tmp_1_addr, &n_tmp_1, &m_tmp_1, &tmp_1); SCICOINOR_ERROR;
	  
	  if (n_tmp_1*m_tmp_1!=nele_hess)
	    {
	      Scierror(999,"sciipopt: the Hessian function must return %d values\n",nele_hess);
	      sci_parameters->error = 1;
	      Nbvars = nbvars_old;
	      return 0;
	    }
	  
	  memcpy(values, tmp_1, nele_hess*sizeof(double));
	}
      
      Nbvars = nbvars_old;

      return (int)true;
    }
  else
    {
      Nbvars = nbvars_old;

      return (int)false;
    }
}
