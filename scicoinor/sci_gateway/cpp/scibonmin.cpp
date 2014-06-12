///////////////////////////////////////////////////////////////////////////////////////////////////////////
// scibonmin: a scilab interface to Bonmin, a tool for non linear constrained mixed integer optimization //
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//  Copyright (C) 2008-2010 Yann Collette.
//
//  SCIBONMIN is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2, or (at your option) any
//  later version.
//
//  This part of code is distributed with the FURTHER condition that it 
//  can be compiled and linked with the Matlab libraries and it can be 
//  used within the Matlab environment.
//
//  SCIBONMIN is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with SciCoinOR; see the file COPYING.  If not, write to the Free
//  Software Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

/////////////////////////////
//  List of bonmin options //
/////////////////////////////

// From the documentation:

////////////////////////////////////////////////////////////////////////////////////////
//  Option 	                type 	default 	B-BB 	B-OA 	B-QG 	B-Hyb //
////////////////////////////////////////////////////////////////////////////////////////
//                              Algorithm choice                                      //
////////////////////////////////////////////////////////////////////////////////////////
//  algorithm 	                S 	B-BB 	        + 	+ 	+ 	+     // 
////////////////////////////////////////////////////////////////////////////////////////
//  Bonmin ecp based strong branching                                                 //
////////////////////////////////////////////////////////////////////////////////////////
//  ecp_abs_tol_strong 	        F 	1e-06 	        + 	+ 	+ 	+     //
//  ecp_max_rounds_strong 	I 	0 	        + 	+ 	+ 	+     //
//  ecp_rel_tol_strong 	        F 	0.1 	        + 	+ 	+ 	+     //
//  lp_strong_warmstart_method 	S       Basis   	+ 	+ 	+ 	+     //
////////////////////////////////////////////////////////////////////////////////////////
//  Branch-and-bound options                                                          //
////////////////////////////////////////////////////////////////////////////////////////
//  allowable_fraction_gap 	F 	0 	        + 	+ 	+ 	+     //
//  allowable_gap 	        F 	0 	        + 	+ 	+ 	+     //
//  cutoff 	                F 	1e+100 	        + 	+ 	+ 	+     //
//  cutoff_decr 	        F 	1e-05 	        + 	+ 	+ 	+     //
//  integer_tolerance 	        F 	1e-06 	        + 	+ 	+ 	+     //
//  iteration_limit 	        I 	INT_MAX 	+ 	+ 	+ 	+     //
//  nlp_failure_behavior 	S 	stop 	        + 	+ 	+ 	+     //
//  node_comparison 	        S 	dynamic 	+ 	+ 	+ 	+     //
//  node_limit 	                I       INT_MAX 	+ 	+ 	+ 	+     //
//  num_cut_passes 	        I 	1 	        - 	+ 	+ 	+     //
//  num_cut_passes_at_root 	I 	20 	        - 	+ 	+ 	+     //
//  number_before_trust 	I 	8 	        + 	+ 	+ 	+     //
//  number_strong_branch 	I 	20 	        + 	+ 	+ 	+     //
//  solution_limit 	        I 	INT_MAX 	+ 	+ 	+ 	+     //
//  sos_constraints 	        S 	enable 	        + 	- 	+ 	+     //
//  time_limit 	                F 	1e+10 	        + 	+ 	+ 	+     //
//  tree_search_strategy 	S 	top-node 	+ 	+ 	+ 	+     //
//  variable_selection 	        S 	strong-branching+ 	+ 	+ 	+     //
////////////////////////////////////////////////////////////////////////////////////////
//  Diving options                                                                    //
////////////////////////////////////////////////////////////////////////////////////////
//  max_backtracks_in_dive 	I 	5 	        + 	+ 	+ 	+     //
//  max_dive_depth 	        I 	INT_MAX 	+ 	+ 	+ 	+     //
//  stop_diving_on_cutoff 	S 	no 	        + 	+ 	+ 	+     //
////////////////////////////////////////////////////////////////////////////////////////
//  Feasibility pump heuristic                                                        //
////////////////////////////////////////////////////////////////////////////////////////
//  feasibility_pump_objective_norm I 	1 	        + 	+ 	+ 	+     //
//  heuristic_feasibility_pump 	S 	no 	        + 	+ 	+ 	+     //
////////////////////////////////////////////////////////////////////////////////////////
//  Fractional diving MIP heuristic                                                   //
////////////////////////////////////////////////////////////////////////////////////////
//  heuristic_dive_MIP_fractional S 	no 	        + 	+ 	+ 	+     //
////////////////////////////////////////////////////////////////////////////////////////
//  Fractional diving heuristic                                                       //
////////////////////////////////////////////////////////////////////////////////////////
//  heuristic_dive_fractional 	S 	no 	        + 	+ 	+ 	+     //
////////////////////////////////////////////////////////////////////////////////////////
//  Local search based heuristics                                                     //
////////////////////////////////////////////////////////////////////////////////////////
//  dummy_pump_heuristic 	S 	no 	        + 	+ 	+ 	+     //
//  fix_and_solve_heuristic 	S 	no 	        + 	+ 	+ 	+     //
//  heuristic_RINS 	        S 	no 	        + 	+ 	+ 	+     //
//  heuristic_local_branching 	S 	no 	        + 	+ 	+ 	+     //
//  local_search_node_limit 	I 	1000 	        + 	+ 	+ 	+     //
//  local_search_solution_limit I 	5 	        + 	+ 	+ 	+     //
//  local_search_time_limit 	F 	60 	        + 	+ 	+ 	+     //
////////////////////////////////////////////////////////////////////////////////////////
//  MILP cutting planes in hybrid                                                     //
////////////////////////////////////////////////////////////////////////////////////////
//  2mir_cuts 	                I 	0 	        - 	+ 	- 	+     //
//  Gomory_cuts 	        I 	-5 	        - 	+ 	- 	+     //
//  clique_cuts 	        I 	-5 	        - 	+ 	- 	+     //
//  cover_cuts 	                I 	-5 	        - 	+ 	- 	+     //
//  flow_covers_cuts 	        I 	-5 	        - 	+ 	- 	+     //
//  lift_and_project_cuts 	I 	0 	        - 	+ 	- 	+     //
//  mir_cuts 	                I 	-5 	        - 	+ 	- 	+     //
//  probing_cuts 	        I 	-5 	        - 	+ 	- 	+     //
//  reduce_and_split_cuts 	I 	0 	        - 	+ 	- 	+     //
////////////////////////////////////////////////////////////////////////////////////////
//  Nlp solution robustness                                                           //
////////////////////////////////////////////////////////////////////////////////////////
//  max_consecutive_failures 	I 	10 	        + 	+ 	+ 	+     //
//  max_random_point_radius 	F 	100000 	        + 	- 	- 	-     //
//  num_iterations_suspect 	I 	-1 	        + 	+ 	+ 	+     //
//  num_retry_unsolved_random_point I 	 0 	        + 	+ 	+ 	+     //
//  random_point_perturbation_interval F 1 	        + 	- 	- 	-     //
//  random_point_type 	        S 	Jon 	        + 	- 	- 	-     //
////////////////////////////////////////////////////////////////////////////////////////
//  Nlp solve options in B-Hyb                                                        //
////////////////////////////////////////////////////////////////////////////////////////
//  nlp_solve_frequency 	I 	10 	        - 	- 	- 	+     //
//  nlp_solve_max_depth 	I 	10 	        - 	- 	- 	+     //
//  nlp_solves_per_depth 	F 	1e+30 	        - 	- 	- 	+     //
////////////////////////////////////////////////////////////////////////////////////////
//  Options for MILP subsolver in OA decomposition                                    //
////////////////////////////////////////////////////////////////////////////////////////
//  milp_log_level 	        I 	0 	        - 	+ 	- 	+     //
//  milp_subsolver 	        S 	Cbc_D 	        - 	+ 	- 	+     //
////////////////////////////////////////////////////////////////////////////////////////
//  Options for OA decomposition                                                      //
////////////////////////////////////////////////////////////////////////////////////////
//  oa_log_frequency 	        F 	100 	        + 	+ 	+ 	+     //
//  oa_log_level 	        I 	1 	        + 	+ 	+ 	+     //
////////////////////////////////////////////////////////////////////////////////////////
//  Options for ecp cuts generation                                                   //
////////////////////////////////////////////////////////////////////////////////////////
//  ecp_abs_tol 	        F 	1e-06 	        - 	- 	- 	+     //
//  ecp_max_rounds 	        I 	5 	        - 	- 	- 	+     //
//  ecp_propability_factor 	F 	1000 	        - 	- 	- 	+     //
//  ecp_rel_tol 	        F 	0 	        - 	- 	- 	+     //
//  filmint_ecp_cuts 	        I 	0 	        - 	- 	- 	+     //
////////////////////////////////////////////////////////////////////////////////////////
//  Options for non-convex problems                                                   //
////////////////////////////////////////////////////////////////////////////////////////
//  max_consecutive_infeasible 	I 	0 	        + 	+ 	+ 	+     //
//  num_resolve_at_infeasibles 	I 	0 	        + 	- 	- 	-     //
//  num_resolve_at_node 	I 	0 	        + 	- 	- 	-     //
//  num_resolve_at_root 	I 	0 	        + 	- 	- 	-     //
////////////////////////////////////////////////////////////////////////////////////////
//  Outer Approximation cuts generation                                               //
////////////////////////////////////////////////////////////////////////////////////////
//  add_only_violated_oa 	S 	no 	        - 	+ 	+ 	+     //
//  cut_strengthening_type 	S 	none 	        - 	+ 	+ 	+     //
//  disjunctive_cut_type 	S 	none 	        - 	+ 	+ 	+     //
//  oa_cuts_log_level 	        I 	0 	        - 	+ 	+ 	+     //
//  oa_cuts_scope 	        S 	global 	        - 	+ 	+ 	+     //
//  tiny_element 	        F 	1e-08 	        - 	+ 	+ 	+     //
//  very_tiny_element 	        F 	1e-17 	        - 	+ 	+ 	+     //
////////////////////////////////////////////////////////////////////////////////////////
//  Output ond log-levels options                                                     //
////////////////////////////////////////////////////////////////////////////////////////
//  bb_log_interval 	        I 	100 	        + 	- 	+ 	+     //
//  bb_log_level 	        I 	1 	        + 	- 	+ 	+     //
//  lp_log_level 	        I 	0 	        - 	- 	+ 	+     //
////////////////////////////////////////////////////////////////////////////////////////
//  Strong branching setup                                                            //
////////////////////////////////////////////////////////////////////////////////////////
//  candidate_sort_criterion 	S 	best-ps-cost 	+ 	+ 	+ 	+     //
//  maxmin_crit_have_sol 	F 	0.1 	        + 	+ 	+ 	+     //
//  maxmin_crit_no_sol 	        F 	0.7 	        + 	+ 	+ 	+     //
//  min_number_strong_branch 	I 	0 	        + 	+ 	+ 	+     //
//  number_before_trust_list 	I 	0 	        + 	+ 	+ 	+     //
//  number_look_ahead 	        I 	0 	        + 	+ 	+ 	+     //
//  number_strong_branch_root 	I 	INT_MAX 	+ 	+ 	+ 	+     //
//  setup_pseudo_frac 	        F 	0.5 	        + 	+ 	+ 	+     //
//  trust_strong_branching_for_pseudo_cost S 	yes 	+ 	+ 	+ 	+     //
////////////////////////////////////////////////////////////////////////////////////////
//  VectorLength diving MIP heuristic                                                 //
////////////////////////////////////////////////////////////////////////////////////////
//  heuristic_dive_MIP_vectorLength S 	no 	        + 	+ 	+ 	+     //
////////////////////////////////////////////////////////////////////////////////////////
//  VectorLength diving heuristic                                                     //
////////////////////////////////////////////////////////////////////////////////////////
//  heuristic_dive_vectorLength S 	no 	        + 	+ 	+ 	+     //
////////////////////////////////////////////////////////////////////////////////////////
//  nlp interface option                                                              //
////////////////////////////////////////////////////////////////////////////////////////
//  file_solution 	        S 	no 	        + 	+ 	+ 	+     //
//  nlp_log_level 	        I 	1 	        + 	+ 	+ 	+     //
//  nlp_solver 	                S 	Ipopt 	        + 	+ 	+ 	+     //
//  warm_start 	                S 	none 	        + 	- 	- 	-     //
////////////////////////////////////////////////////////////////////////////////////////

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <iomanip>
#include <fstream>
#include <string>
#include <new>
#include <exception>
#include <algorithm>

#undef min
#undef max

#include <CoinTime.hpp>
#include <CoinError.hpp>

#include <BonTMINLP.hpp>
#include <BonOsiTMINLPInterface.hpp>
#include <BonIpoptSolver.hpp>
#include <BonCbc.hpp>
#include <BonBonminSetup.hpp>

#include <BonOACutGenerator2.hpp>
#include <BonEcpCuts.hpp>
#include <BonOaNlpOptim.hpp>

#include <scilabjournal.hpp>
#include <call_function.hpp>

extern "C"
{
#include <stack-c.h>
#include <sciprint.h>
#include <Scierror.h>
#include <MALLOC.h>
#include <api_parameters.h>
}

#include <api_scilab.h>

#include <helper.hpp>
#include <def_solvers.hpp>
#include <manage_params.hpp>

//#define DEBUG 1
//#define DBGPRINTF printf
#define DBGPRINTF sciprint

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

static FTAB FTab_bonmin_call_f[] = {{(char *) 0, (voidf) 0}};

int sci_bonmin_objective(double * x, double * f, int n_size_x, double x_new, void * param);
int sci_bonmin_objective_grad(double * x, double * f, int n_size_x, double x_new, void * param);
int sci_bonmin_constraints(double * x, int n_size_x, double * g, int n_size_g, double x_new, void * param);
int sci_bonmin_constraints_jac(double * x, int n_size_x, double new_x, int n_size_g, int nnz_jac, 
			       int * iRow, int * jCol, double * values, void * param);
int sci_bonmin_hessian(double * x, int n_size_x, double new_x, double obj_factor, int n_size_g, double * lambda,
		       double new_lambda, int nele_hess, int * iRow, int * jCol, double * values, void * param);

#define X_POS                   1
#define FOBJ_POS                2
#define DFOBJ_POS               3
#define GOBJ_POS                4
#define DGOBJ_POS               5
#define SPARSE_DGOBJ_POS        6
#define DHOBJ_POS               7
#define SPARSE_DHOBJ_POS        8
#define VAR_TYPE_POS            9
#define VAR_LIN_TYPE_POS        10
#define CONSTR_LIN_TYPE_POS     11
#define CONSTR_RHS_POS          12
#define CONSTR_LHS_POS          13
#define X_LOWER_POS             14
#define X_UPPER_POS             15
#define PARAMS_POS              16
#define PARAM_IN                PARAMS_POS
#define LAST_PARAMS             PARAMS_POS
#define X_SOL_OUT_POS           Rhs+1
#define F_SOL_OUT_POS           Rhs+2
#define EXTRA_OUT_POS           Rhs+3

using namespace Ipopt;
using namespace Bonmin;

// A structure which handles informations related to objective functions
struct sci_bonmin_info
{
  int fobj_lhs,  fobj_rhs,  l_fobj,  fobj_is_list;
  int dfobj_lhs, dfobj_rhs, l_dfobj, dfobj_is_list;
  int gobj_lhs,  gobj_rhs,  l_gobj,  gobj_is_list;
  int dgobj_lhs, dgobj_rhs, l_dgobj, dgobj_is_list;
  int dhobj_lhs, dhobj_rhs, l_dhobj, dhobj_is_list;
  int m_sparse_dgobj, n_sparse_dgobj; // the sparsity structure of the jacobian
  int m_sparse_dhobj, n_sparse_dhobj; // the sparsity structure of the Hessian
  double * sparse_dgobj;
  double * sparse_dhobj;
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
  int error;             // ==1 if an error occured
  int ibegin;            // the position of the top of the stack
  int nbvar;             // number of variables
  int nbconstr;          // number of constraints
  int nnz_jac_g;         // number of non zeroes in Jacobian
  int nnz_h_lag;         // number of non zeroes in Hessian of Lagrangean
  int index_style;       // if the index starts from 0 (C) or 1 (FORTRAN)
  double * var_type;        // position on the stack of the type of variables (int)
  double * var_lin_type;    // position on the stack of the type of variables linearity (int)
  double * constr_lin_type; // position on the stack of the type of constraints linearity (int)
  double * constr_rhs;   // position on the stack of the upper bound for constraints (double)
  double * constr_lhs;      // position on the stack of the lower bound for constraints (double)
  double * x_lower;         // position on the stack of the lower bound for the variables (double)
  double * x_upper;         // position on the stack of the upper bound for the variables (double)
  double * x_0;             // position on the stack of the starting point (double)
  BonminSetup * bonminSetup;
  Bab         * bab;
  bool          HessianPresent;
  double * pdbl_x_sol;
  double f_sol;
  double status;
  double mip_status;
  double bestbound;
  double numnodes;
  double itercount;
  double bonminstat;
};

class MyTMINLP : public TMINLP
{
public:
  MyTMINLP() : printSol_(false), sci_parameters(NULL) {}
  virtual ~MyTMINLP() {}
  MyTMINLP(const MyTMINLP &other) : printSol_(other.printSol_), sci_parameters(other.sci_parameters) {}
  virtual bool get_variables_types(Index n, VariableType* var_types);
  virtual bool get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types);
  virtual bool get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types);
  virtual bool get_nlp_info(Index& n, Index&m, Index& nnz_jac_g, Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style);
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u);
  virtual bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda);
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index *jCol, Number* values);
  virtual bool eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda, bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values);
  virtual void finalize_solution(TMINLP::SolverReturn status, Index n, const Number* x, Number obj_value);
  virtual const SosInfo * sosConstraints() const{return NULL;}
  virtual const BranchingInfo* branchingInfo() const{return NULL;}
  void printSolutionAtEndOfAlgorithm() {printSol_ = true;}
  void set_scilab_parameters(struct sci_bonmin_info * param) {sci_parameters = param;}
  struct sci_bonmin_info * get_scilab_parameters() {return sci_parameters;}
private:
  bool printSol_;
  struct sci_bonmin_info * sci_parameters;
};

bool MyTMINLP::get_variables_types(Index n, VariableType* var_types)
{
#ifdef DEBUG
  DBGPRINTF("scibonmin: get_variables_type\n");
#endif

  // CONTINUOUS 0
  // BINARY     1
  // INTEGER    2

  int i;
  
  if (sci_parameters->var_type)
    {
      for(i=0;i<n;i++) 
        {
          switch((int)sci_parameters->var_type[i])
            {
            case 0:
              var_types[i] = CONTINUOUS;
              break;
            case 1:
              var_types[i] = BINARY;
              break;
            default:
              var_types[i] = INTEGER;
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

bool MyTMINLP::get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_lin_types)
{
#ifdef DEBUG
  DBGPRINTF("scibonmin: get_variables_linearity\n");
#endif

  int i;

  if (sci_parameters->var_lin_type)
    {
      for(i=0;i<n; i++) 
        {
          switch((int)sci_parameters->var_lin_type[i])
            {
            case 0:
              var_lin_types[i] = Ipopt::TNLP::LINEAR;
              break;
            default:
              var_lin_types[i] = Ipopt::TNLP::NON_LINEAR;
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

bool MyTMINLP::get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_lin_types)
{
#ifdef DEBUG
  DBGPRINTF("scibonmin: get_constraints_linearity\n");
#endif

  int i;

  if (sci_parameters->constr_lin_type)
    {
      for(i=0;i<m;i++) 
        {
          switch((int)sci_parameters->constr_lin_type[i])
            {
            case 0:
              const_lin_types[i] = Ipopt::TNLP::LINEAR;
              break;
            default:
              const_lin_types[i] = Ipopt::TNLP::NON_LINEAR;
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

bool MyTMINLP::get_nlp_info(Index& n, Index&m, Index& nnz_jac_g,Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
{
#ifdef DEBUG
  DBGPRINTF("scibonmin: get_nlp_info\n");
#endif

  n = sci_parameters->nbvar;    // number of variable
  m = sci_parameters->nbconstr; // number of constraints
  nnz_jac_g = sci_parameters->nnz_jac_g; // number of non zeroes in Jacobian
  nnz_h_lag = sci_parameters->nnz_h_lag; // number of non zeroes in Hessian of Lagrangean
  index_style = (sci_parameters->index_style==0) ? TNLP::C_STYLE : TNLP::FORTRAN_STYLE;

  return true;
}

bool MyTMINLP::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
{
#ifdef DEBUG
  DBGPRINTF("scibonmin: get_bounds_info\n");
#endif

  int i;
  for(i=0;i<n;i++)
    {
      x_l[i] = sci_parameters->x_lower[i];
      x_u[i] = sci_parameters->x_upper[i];
    }

  for(i=0;i<m;i++)
    {
      g_l[i] = sci_parameters->constr_lhs[i];
      g_u[i] = sci_parameters->constr_rhs[i];
    }

  return true;
}

bool MyTMINLP::get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda)
{
#ifdef DEBUG
  DBGPRINTF("scibonmin: get_starting_point\n");
#endif

  int i;

  if (sci_parameters->x_0)
    {
      for(i=0;i<n;i++) x[i] = sci_parameters->x_0[i];
      
      return true;
    }
  else
    {
      return false;
    }
}

bool MyTMINLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
#ifdef DEBUG
  DBGPRINTF("scibonmin: eval_f\n");
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

bool MyTMINLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
#ifdef DEBUG
  DBGPRINTF("scibonmin: eval_grad_f\n");
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

bool MyTMINLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
#ifdef DEBUG
  DBGPRINTF("scibonmin: eval_g\n");
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

bool MyTMINLP::eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nnz_jac, Index* iRow, Index *jCol, Number* values)
{
#ifdef DEBUG
  DBGPRINTF("scibonmin: eval_jac_g\n");
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

bool MyTMINLP::eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda,
		      bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values)
{
#ifdef DEBUG
  DBGPRINTF("scibonmin: eval_h\n");
#endif

  if (sci_parameters->hessian)
    {
      return (int)(*sci_parameters->hessian)((double *)x, (int)n, (double)new_x, (double)obj_factor, (int)m, 
					     (double *)lambda, (double)new_lambda, 
					     (int)nele_hess, (int *)iRow, (int *)jCol, (double *)values,
					     (void *)sci_parameters);
    }
  else
    {
      return false;
    }
}

void MyTMINLP::finalize_solution(TMINLP::SolverReturn status, Index n, const Number* x, Number obj_value)
{
#ifdef DEBUG
  DBGPRINTF("scibonmin: finalize_solution\n");
#endif

  if (sci_parameters->error==1) status = (Bonmin::TMINLP::SolverReturn)INFEASIBLE;
  
  // SUCCESS 	    0
  // INFEASIBLE     1
  // LIMIT_EXCEEDED 2
  // MINLP_ERROR    3

  if (status == (Bonmin::TMINLP::SolverReturn)SUCCESS)
    {
      memcpy(sci_parameters->pdbl_x_sol,x,n*sizeof(double));
      
      sci_parameters->f_sol = obj_value;
    }
  else
    {
      memcpy(sci_parameters->pdbl_x_sol,sci_parameters->x_0,n*sizeof(double));
      
      sci_parameters->f_sol = 0;
    }
  
  // Status values:
  // SUCCESS 	    0
  // INFEASIBLE     1
  // LIMIT_EXCEEDED 2
  // MINLP_ERROR    3

  if (sci_parameters->error==1)
    {
      sci_parameters->status     = -100;
      sci_parameters->mip_status = -100;
      sci_parameters->bestbound  = 0;
      sci_parameters->numnodes   = 0;
      sci_parameters->itercount  = 0;
      sci_parameters->bonminstat = -100;
    }
  else
    {
      sci_parameters->status = (double)status;
      
      // We get some informations via Bab:
      // MipStatuses mipStatus() const       : return Mip Status 
      // double bestBound()                  : return the best known lower bound on the objective value
      // int numNodes() const                : return the total number of nodes explored
      // int iterationCount()                : return the total number of iterations in the last mip solved. 
      // const double * bestSolution() const : get the best solution known to the problem (is optimal if MipStatus is FeasibleOptimal).
      // double bestObj() const              : return objective value of the bestSolution 
      
      // MipStatuses:
      // 0 - FeasibleOptimal  Optimum solution has been found and its optimality proved.
      // 1 - ProvenInfeasible Problem has been proven to be infeasible.
      // 2 - Feasible         An integer solution to the problem has been found.
      // 3 - NoSolutionKnown  No feasible solution to the problem is known. 
      
      sci_parameters->mip_status = (double)sci_parameters->bab->mipStatus();
      sci_parameters->bestbound  = (double)sci_parameters->bab->bestBound();
      sci_parameters->numnodes   = (double)sci_parameters->bab->numNodes();
      sci_parameters->itercount  = (double)sci_parameters->bab->iterationCount();
      
      // Via bonmin.nonlinearSolver()->is....
      // bool isAbandoned() const                    : Are there a numerical difficulties?
      // bool isProvenOptimal () const               : Is optimality proven?
      // bool isProvenPrimalInfeasible () const      : Is primal infeasiblity proven?
      // bool isProvenDualInfeasible () const        : Is dual infeasiblity proven?
      // bool isPrimalObjectiveLimitReached () const : Is the given primal objective limit reached?
      // bool isDualObjectiveLimitReached () const   : Is the given dual objective limit reached?
      // bool isIterationLimitReached () const       : Iteration limit reached?
      
      sci_parameters->bonminstat  = (double)pow(2.0,1)*sci_parameters->bonminSetup->nonlinearSolver()->isAbandoned();
      sci_parameters->bonminstat += (double)pow(2.0,2)*sci_parameters->bonminSetup->nonlinearSolver()->isProvenOptimal();
      sci_parameters->bonminstat += (double)pow(2.0,3)*sci_parameters->bonminSetup->nonlinearSolver()->isProvenPrimalInfeasible();
      sci_parameters->bonminstat += (double)pow(2.0,4)*sci_parameters->bonminSetup->nonlinearSolver()->isProvenDualInfeasible();
      //sci_parameters->bonminstat += (double)pow(2.0,5)*sci_parameters->bonminSetup->nonlinearSolver()->isPrimalObjectiveLimitReached(); 
      sci_parameters->bonminstat += (double)pow(2.0,6)*sci_parameters->bonminSetup->nonlinearSolver()->isDualObjectiveLimitReached();
      sci_parameters->bonminstat += (double)pow(2.0,7)*sci_parameters->bonminSetup->nonlinearSolver()->isIterationLimitReached();
    }
}
  
extern "C" int scibonmin(char * fname)
{
  // x0 fobj, dfobj, gobj, dgobj, params ...
  int m_x,               n_x,               * x_addr;
  int m_fobj,            n_fobj,            l_fobj,  fobj_is_list;
  int m_dfobj,           n_dfobj,           l_dfobj, dfobj_is_list;
  int m_gobj,            n_gobj,            l_gobj,  gobj_is_list;
  int m_dgobj,           n_dgobj,           l_dgobj, dgobj_is_list;
  int m_dhobj,           n_dhobj,           l_dhobj, dhobj_is_list;
  int m_sparse_dgobj,    n_sparse_dgobj,    * sparse_dgobj_addr = NULL; // sparsity structure of dgobj
  int m_sparse_dhobj,    n_sparse_dhobj,    * sparse_dhobj_addr = NULL; // sparsity structure of dhobj
  int m_var_type,        n_var_type,        * var_type_addr = NULL;
  int m_var_lin_type,    n_var_lin_type,    * var_lin_type_addr = NULL;
  int m_constr_lin_type, n_constr_lin_type, * constr_lin_type_addr = NULL;
  int m_constr_rhs,      n_constr_rhs,      * constr_rhs_addr = NULL;
  int m_constr_lhs,      n_constr_lhs,      * constr_lhs_addr = NULL;
  int m_x_lower,         n_x_lower,         * x_lower_addr = NULL;
  int m_x_upper,         n_x_upper,         * x_upper_addr = NULL;
  double * x = NULL, * x_lower = NULL, * x_upper = NULL;
  double * constr_rhs = NULL, * constr_lhs = NULL;
  double * var_type = NULL, * var_lin_type = NULL, * constr_lin_type = NULL;
  double * sparse_dgobj = NULL, * sparse_dhobj = NULL;
  int Log = 0;
  bool nnz_jac_is_needed = false, nnz_hess_is_needed = false;
  int * param_in_addr = NULL;
  struct sci_bonmin_info * sci_parameters = new struct sci_bonmin_info;
  DerivedHandler * printer = NULL;
  SciErr _SciErr;

  CheckRhs(LAST_PARAMS,LAST_PARAMS);
  CheckLhs(3,3);

  if (Rhs<LAST_PARAMS) 
    {
      Scierror(999,"%s: %d inputs required in call to %s. Bug in bonmin.sci ?...\n",fname, fname, LAST_PARAMS);
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
    (call_f_objective)GetFunctionPtr("bonmin: objective", FOBJ_POS, FTab_bonmin_call_f, 
				     (voidf)sci_bonmin_objective, 
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
	(call_f_objective_grad)GetFunctionPtr("bonmin: objective_grad", DFOBJ_POS, FTab_bonmin_call_f, 
					      (voidf)sci_bonmin_objective_grad, 
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
	(call_f_constraints)GetFunctionPtr("bonmin: constraints", GOBJ_POS, FTab_bonmin_call_f, 
					   (voidf)sci_bonmin_constraints, 
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
	(call_f_constraints_jac)GetFunctionPtr("bonmin: constraints_jac", DGOBJ_POS, FTab_bonmin_call_f, 
					       (voidf)sci_bonmin_constraints_jac, 
					       &l_dgobj, &m_dgobj, &n_dgobj, &dgobj_is_list,
					       sci_parameters->dgobj_parameter_list);

      // If we pass a C function here, we need the 'nnz_jac' option which
      // specifies the number of non zeros elements in the Jacobian
      if (sci_parameters->constraints_jac!=sci_bonmin_constraints_jac)
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
	(call_f_hessian)GetFunctionPtr("bonmin: hessian", DHOBJ_POS, FTab_bonmin_call_f, 
				       (voidf)sci_bonmin_hessian, 
				       &l_dhobj, &m_dhobj, &n_dhobj, &dhobj_is_list,
				       sci_parameters->dhobj_parameter_list);

      // If we pass a C function here, we need the 'nnz_hess' option which
      // specifies the number of non zeros elements in the Hessian
      if (sci_parameters->hessian!=sci_bonmin_hessian)
        {
          nnz_hess_is_needed = true;
        }
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, SPARSE_DHOBJ_POS, &sparse_dhobj_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, sparse_dhobj_addr, &m_sparse_dhobj, &n_sparse_dhobj, &sparse_dhobj); SCICOINOR_ERROR;

  _SciErr = getVarAddressFromPosition(pvApiCtx, VAR_TYPE_POS, &var_type_addr); SCICOINOR_ERROR;
  if (~isEmptyMatrix(pvApiCtx, var_type_addr))
    {
      _SciErr = getMatrixOfDouble(pvApiCtx, var_type_addr, &m_var_type, &n_var_type, &var_type); SCICOINOR_ERROR;
    }
  else
    {
      var_type = NULL;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, VAR_LIN_TYPE_POS, &var_lin_type_addr); SCICOINOR_ERROR;
  if (~isEmptyMatrix(pvApiCtx, var_lin_type_addr))
    {
      _SciErr = getMatrixOfDouble(pvApiCtx, var_lin_type_addr, &m_var_lin_type, &n_var_lin_type, &var_lin_type); SCICOINOR_ERROR;
    }
  else
    {
      var_lin_type = NULL;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, CONSTR_LIN_TYPE_POS, &constr_lin_type_addr); SCICOINOR_ERROR;
  if (~isEmptyMatrix(pvApiCtx, constr_lin_type_addr))
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

  if ((m_x_upper!=m_var_type)&&(n_x_upper!=n_var_type))
    {
      Scierror(999,"%s: x_upper and var_lin_type must have the same size.\n",fname);
      return 0;
    }

  sci_parameters->fobj_lhs      = m_fobj;
  sci_parameters->fobj_rhs      = n_fobj;
  sci_parameters->l_fobj        = l_fobj;
  sci_parameters->fobj_is_list  = fobj_is_list;

  sci_parameters->dfobj_lhs     = m_dfobj;
  sci_parameters->dfobj_rhs     = n_dfobj;
  sci_parameters->l_dfobj       = l_dfobj;
  sci_parameters->dfobj_is_list = dfobj_is_list;

  sci_parameters->gobj_lhs     = m_gobj; 
  sci_parameters->gobj_rhs     = n_gobj;
  sci_parameters->l_gobj       = l_gobj;
  sci_parameters->gobj_is_list = gobj_is_list;
  
  sci_parameters->dgobj_lhs     = m_dgobj; 
  sci_parameters->dgobj_rhs     = n_dgobj;
  sci_parameters->l_dgobj       = l_dgobj;
  sci_parameters->dgobj_is_list = dgobj_is_list;

  sci_parameters->dhobj_lhs     = m_dhobj; 
  sci_parameters->dhobj_rhs     = n_dhobj;
  sci_parameters->l_dhobj       = l_dhobj;
  sci_parameters->dhobj_is_list = dhobj_is_list;

  if ((n_sparse_dgobj!=2)&&(n_sparse_dgobj))
    {
      Scierror(999,"%s: sparse_dgobj must be a mx2 matrix\n",fname);
      return 0;
    }

  if (sci_parameters->HessianPresent)
    {
      if ((n_sparse_dhobj!=2)&&(n_sparse_dhobj!=0))
	{
	  Scierror(999,"%s: sparse_dhobj must be a mx2 matrix\n",fname);
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

  sci_parameters->ibegin          = Rhs + 1;
  sci_parameters->nbvar           = m_x*n_x;

  if (sci_parameters->constraints)
    sci_parameters->nbconstr        = m_constr_rhs * n_constr_rhs;
  else
    sci_parameters->nbconstr        = 0;

  sci_parameters->index_style     = 1; // C_STYLE=0, FORTRAN_STYLE=1
  sci_parameters->var_type        = var_type;
  sci_parameters->var_lin_type    = var_lin_type;
  sci_parameters->constr_lin_type = constr_lin_type;
  sci_parameters->constr_rhs      = constr_rhs;
  sci_parameters->constr_lhs      = constr_lhs;
  sci_parameters->x_lower         = x_lower;
  sci_parameters->x_upper         = x_upper;
  sci_parameters->x_0             = x;

  int m_x_sol       = 1,  n_x_sol      = 1;
  int m_f_sol       = 1,  n_f_sol      = 1;
  // Create the 'extra' structure of type plist
  int m_status      = 1, n_status      = 1;
  int m_mip_status  = 1, n_mip_status  = 1;
  int m_bestbound   = 1, n_bestbound   = 1;
  int m_numnodes    = 1, n_numnodes    = 1;
  int m_itercount   = 1, n_itercount   = 1;
  int m_bonminstat  = 1, n_bonminstat  = 1;
  int m_list_labels = 1, n_list_labels = 7;

  m_x_sol = sci_parameters->nbvar;

  sci_parameters->pdbl_x_sol = (double *)MALLOC(m_x_sol*sizeof(double));

  sci_parameters->f_sol      = 0;
  sci_parameters->status     = 0;
  sci_parameters->mip_status = 0;
  sci_parameters->bestbound  = 0;
  sci_parameters->numnodes   = 0;
  sci_parameters->itercount  = 0;
  sci_parameters->bonminstat = 0;

  SmartPtr<MyTMINLP> myminlp = new MyTMINLP;

  BonminSetup * bonmin = NULL;
  bonmin = new BonminSetup;

  sci_parameters->bonminSetup = bonmin;

  myminlp->set_scilab_parameters(sci_parameters);

  // Get the parameters stored in the plist
  int tmp_res, tmp_int;

  _SciErr = initPList(pvApiCtx, PARAM_IN, &param_in_addr); SCICOINOR_ERROR;
  if (!checkPList(pvApiCtx, param_in_addr))
    {
      Scierror(999, "%s: argument n° %d is not a plist\n", fname, PARAM_IN);
      return 0;
    }

#ifdef DEBUG
  DBGPRINTF("scibonmin: processing options\n");
#endif

  //Here we read several option files
  bonmin->initializeOptionsAndJournalist();

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

  // This reads the default file "bonmin.opt"
  //bonmin->readOptionsFile(); // YC: ou mettre cette fonction: en début ou en fin de lecture d'options ??

  // Add a new journal to print all the bonmin informations in the scilab console
  _SciErr = getIntInPList(pvApiCtx, param_in_addr, "journal_level", &tmp_int, &tmp_res, 0, Log, CHECK_NONE); SCICOINOR_ERROR;
  SmartPtr<Journal> console = new ScilabJournal((EJournalLevel)tmp_int);
  // For Bonmin
  bonmin->journalist()->AddJournal(console);

  // Manage the parameters for Bonmin and Ipopt
  manage_bonmin_params(bonmin->options(), param_in_addr, Log);
  manage_ipopt_params(bonmin->options(), param_in_addr, Log);

#ifdef DEBUG
  bonmin->options()->SetStringValue("print_options_documentation", "yes");
  bonmin->mayPrintDoc();
#endif

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

  //Now initialize from myminlp
  bonmin->initialize(GetRawPtr(myminlp));

  bonmin->options()->GetIntegerValue("bb_log_level",tmp_int,"bonmin");
  printer = new DerivedHandler();
  printer->setLogLevel(tmp_int);
  bonmin->nonlinearSolver()->passInMessageHandler(printer);
  bonmin->continuousSolver()->passInMessageHandler(printer);

  Bab * bb = NULL;
  try 
    {
      bb = new Bab;
      sci_parameters->bab = bb;
      (*bb)(bonmin); //process parameter file using Ipopt and do branch and bound using Cbc
    }
  catch(TNLPSolver::UnsolvedError *) 
    {
      //There has been a failure to solve a problem with Ipopt.
      sciprint("scibonmin: Ipopt has failed to solve a problem\n");
      if (sci_parameters) delete sci_parameters;
      if (bonmin) delete bonmin;
      if (bb) delete bb;
      return 0;
    }
  catch(OsiTMINLPInterface::SimpleError &E) 
    {
      sciprint("OsiTMINLPInterface::SimpleError - %s::%s - %s\n", E.className().c_str(), E.methodName().c_str(), E.message().c_str());
      if (sci_parameters) delete sci_parameters;
      if (bonmin) delete bonmin;
      if (bb) delete bb;
      return 0;
    }
  catch(CoinError &E) 
    {
      sciprint("CoinError - %s::%s - %s\n", E.className().c_str(), E.methodName().c_str(), E.message().c_str());
      if (sci_parameters) delete sci_parameters;
      if (bonmin) delete bonmin;
      if (bb) delete bb;
      return 0;
    }

  _SciErr = createMatrixOfDouble(pvApiCtx, X_SOL_OUT_POS, m_x_sol, n_x_sol, sci_parameters->pdbl_x_sol); SCICOINOR_ERROR;
  _SciErr = createMatrixOfDouble(pvApiCtx, F_SOL_OUT_POS, m_f_sol, n_f_sol, &sci_parameters->f_sol); SCICOINOR_ERROR;

  int * extra_addr = NULL;
  char * ListLabels [] = {"status","mip_status","best_bound","num_nodes","iter_count","bonmin_status"};

  _SciErr = createPList(pvApiCtx, EXTRA_OUT_POS, &extra_addr, (char **)ListLabels, 6); SCICOINOR_ERROR;

  _SciErr = createIntInPList(pvApiCtx, EXTRA_OUT_POS,    extra_addr, "status",        sci_parameters->status); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx, EXTRA_OUT_POS,    extra_addr, "mip_status",    sci_parameters->mip_status); SCICOINOR_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, EXTRA_OUT_POS, extra_addr, "best_bound",    sci_parameters->bestbound); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx, EXTRA_OUT_POS,    extra_addr, "num_nodes",     sci_parameters->numnodes); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx, EXTRA_OUT_POS,    extra_addr, "iter_count",    sci_parameters->itercount); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx, EXTRA_OUT_POS,    extra_addr, "bonmin_status", sci_parameters->bonminstat); SCICOINOR_ERROR;

  LhsVar(1) = X_SOL_OUT_POS;
  LhsVar(2) = F_SOL_OUT_POS;
  LhsVar(3) = EXTRA_OUT_POS;

  // Now destroy some pointer and quit
  if (sci_parameters->pdbl_x_sol) FREE(sci_parameters->pdbl_x_sol);
  if (sci_parameters) delete sci_parameters;
  if (bonmin)         delete bonmin;
  if (bb)             delete bb;
  if (printer)        delete printer;

  return 0;
}

//
// Definition of the "fake" objective and constraints.
// These functions allow to be able to deal with scilab scripts or C functions
//

int sci_bonmin_objective(double * x, double * f, int n_size_x, double x_new, void * param)
{
  int n_tmp_1, m_tmp_1, * tmp_1_addr = NULL;
  int n_tmp_2, m_tmp_2;
  double * tmp_1 = NULL, * tmp_2 = NULL;
  int nbvars_old = Nbvars;
  struct sci_bonmin_info * sci_parameters = (struct sci_bonmin_info *)param;
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
      Scierror(999,"scibonmin: error when calling objective function\n");
      sci_parameters->error = 1;
      Nbvars = nbvars_old;
      return 0;
    }

  if (Err>0) 
    {
      Scierror(999,"scibonmin: error when calling objective function\n");
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

int sci_bonmin_objective_grad(double * x, double * df, int n_size_x, double x_new, void * param)
{
  int n_tmp_1, m_tmp_1, * tmp_1_addr = NULL;
  int n_tmp_2, m_tmp_2;
  int nbvars_old = Nbvars;
  double * tmp_1 = NULL, * tmp_2 = NULL;
  int var_type;
  struct sci_bonmin_info * sci_parameters = (struct sci_bonmin_info *)param;
  SciErr _SciErr;

  ////////////////////////////////////////////////
  // Call to the gradient of objective function //
  ////////////////////////////////////////////////

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
      Scierror(999,"scibonmin: error when calling gradient of objective function\n");
      sci_parameters->error = 1;
      Nbvars = nbvars_old;
      return 0;
    }

  if (Err>0) 
    {
      Scierror(999,"scibonmin: error when calling gradient of objective function\n");
      sci_parameters->error = 1;
      Nbvars = nbvars_old;
      return 0;
    } /* End If */

  _SciErr = getVarAddressFromPosition(pvApiCtx, sci_parameters->ibegin, &tmp_1_addr); SCICOINOR_ERROR;
  _SciErr = getVarType(pvApiCtx, tmp_1_addr, &var_type); SCICOINOR_ERROR;

  if (var_type==sci_sparse)
    {
      Scierror(999,"scibonmin: the gradient of the objective function must be not sparse\n");
      sci_parameters->error = 1;
      Nbvars = nbvars_old;
      return 0;
    }

  // Get dfobj
  _SciErr = getMatrixOfDouble(pvApiCtx, tmp_1_addr, &n_tmp_1, &m_tmp_1, &tmp_1); SCICOINOR_ERROR;

  memcpy(df, tmp_1, n_size_x*sizeof(double));

#ifdef DEBUG
  DBGPRINTF("scibonmin: leaving eval_grad_f\n");
#endif

  Nbvars = nbvars_old;

  return (int)true;
}

int sci_bonmin_constraints(double * x, int n_size_x, double * g, int n_size_g, double x_new, void * param)
{
#ifdef DEBUG
  DBGPRINTF("scibonmin: in eval_g\n");
#endif

  int n_tmp_1, m_tmp_1, * tmp_1_addr = NULL;
  int n_tmp_2, m_tmp_2;
  double * tmp_1 = NULL, * tmp_2 = NULL;
  int nbvars_old = Nbvars;
  struct sci_bonmin_info * sci_parameters = (struct sci_bonmin_info *)param;
  SciErr _SciErr;

  if (sci_parameters->nbconstr==0) return (int)false;

  /////////////////////////////////////
  // Call to the constraint function //
  /////////////////////////////////////

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
      Scierror(999,"scibonmin: error when calling constraint function\n");
      sci_parameters->error = 1;
      Nbvars = nbvars_old;
      return 0;
    }

  if (Err>0) 
    {
      Scierror(999,"scibonmin: error when calling constraint function\n");
      sci_parameters->error = 1;
      Nbvars = nbvars_old;
      return 0;
    } /* End If */

  // Get gobj
  _SciErr = getVarAddressFromPosition(pvApiCtx, sci_parameters->ibegin, &tmp_1_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, tmp_1_addr, &n_tmp_1, &m_tmp_1, &tmp_1); SCICOINOR_ERROR;

  memcpy(g, tmp_1, n_size_g*sizeof(double));

  Nbvars = nbvars_old;

  return (int)true;
}

int sci_bonmin_constraints_jac(double * x, int n_size_x, double new_x, int n_size_g, int nnz_jac, 
			       int * iRow, int * jCol, double * values, void * param)
{
#ifdef DEBUG
  DBGPRINTF("scibonmin: in eval_jac_g\n");
#endif

  int n_tmp_1, m_tmp_1, * tmp_1_addr = NULL;
  int n_tmp_2, m_tmp_2;
  double * tmp_1 = NULL, * tmp_2 = NULL;
  int i, Index, var_type;
  int nbvars_old = Nbvars;
  struct sci_bonmin_info * sci_parameters = (struct sci_bonmin_info *)param;
  SciErr _SciErr;

  if (sci_parameters->nbconstr==0) return (int)false;

#ifdef _MSC_VER
  Nbvars = sci_parameters->ibegin + max(sci_parameters->fobj_rhs,sci_parameters->fobj_lhs);
#else
  Nbvars = sci_parameters->ibegin + std::max(sci_parameters->fobj_rhs,sci_parameters->fobj_lhs);
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
      DBGPRINTF("scibonmin: call to the jacobian\n");
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
      DBGPRINTF("scibonmin: call to the function\n");
#endif

      // Call to the scilab function (constraint function)
      try
	{
	  SciFunction(&sci_parameters->ibegin,&sci_parameters->l_dgobj,&sci_parameters->dgobj_lhs,&sci_parameters->dgobj_rhs);
	}
      catch(...)
	{
	  Scierror(999,"scibonmin: error when calling Jacobian of the constraints function\n");
	  sci_parameters->error = 1;
	  Nbvars = nbvars_old;
	  return 0;
	}

      if (Err>0) 
	{
	  Scierror(999,"scibonmin: error when calling Jacobian of the constraints function\n");
	  sci_parameters->error = 1;
	  Nbvars = nbvars_old;
	  return 0;
	} /* End If */
      
      _SciErr = getVarAddressFromPosition(pvApiCtx, sci_parameters->ibegin, &tmp_1_addr); SCICOINOR_ERROR;
      _SciErr = getVarType(pvApiCtx, tmp_1_addr, &var_type); SCICOINOR_ERROR;

      if (var_type==sci_sparse)
	{
	  Scierror(999,"scibonmin: the Jacobian of the constraints must be not sparse\n");
	  sci_parameters->error = 1;
	  Nbvars = nbvars_old;
	  return 0;
	}

      // Get dG
      _SciErr = getMatrixOfDouble(pvApiCtx, tmp_1_addr, &n_tmp_1, &m_tmp_1, &tmp_1); SCICOINOR_ERROR;

      if (n_tmp_1*m_tmp_1!=nnz_jac)
	{
	  Scierror(999,"scibonmin: the Jacobian of the constraint function must return %d values\n",nnz_jac);
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

int sci_bonmin_hessian(double * x, int n_size_x, double new_x, double obj_factor, int n_size_g, double * lambda,
		       double new_lambda, int nele_hess, int * iRow, int * jCol, double * values, void * param)
{
#ifdef DEBUG
  DBGPRINTF("scibonmin: in eval_h\n");
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
  struct sci_bonmin_info * sci_parameters = (struct sci_bonmin_info *)param;
  SciErr _SciErr;

#ifdef _MSC_VER
  Nbvars = sci_parameters->ibegin + max(sci_parameters->fobj_rhs,sci_parameters->fobj_lhs);
#else
  Nbvars = sci_parameters->ibegin + std::max(sci_parameters->fobj_rhs,sci_parameters->fobj_lhs);
#endif

  if (sci_parameters->HessianPresent)
    {
      // Fill iRow and jCol for the sparsity structure of the hessian
      if ((iRow!=NULL)&&(jCol!=NULL)) 
	{
#ifdef DEBUG
	  DBGPRINTF("scibonmin: in eval_h - sparsity\n");
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
	  DBGPRINTF("scibonmin: in eval_h - hessian\n");
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
	  
	  memcpy(tmp_2, lambda, sizeof(double)*n_size_g);
	  
	  // Store the objective function weight into tmp_3
	  n_tmp_3 = 1;
	  m_tmp_3 = 1;

	  _SciErr = allocMatrixOfDouble(pvApiCtx, sci_parameters->ibegin+2, n_tmp_3, m_tmp_3, &tmp_3); SCICOINOR_ERROR;
	  
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
	  DBGPRINTF("scibonmin: in eval_h - call to the function\n");
	  DBGPRINTF("scibonmin: in eval_h - Hessian = %d\n",sci_parameters->HessianPresent);
#endif

	  // Call to the scilab function (Hessian function)
	  try
	    {
	      SciFunction(&sci_parameters->ibegin,&sci_parameters->l_dhobj,&sci_parameters->dhobj_lhs,&sci_parameters->dhobj_rhs);
	    }
	  catch(...)
	    {
	      Scierror(999,"scibonmin: error when calling Hessian function\n");
	      sci_parameters->error = 1;
	      Nbvars = nbvars_old;
	      return 0;
	    }

	  if (Err>0) 
	    {
	      Scierror(999,"scibonmin: error when calling Hessian function\n");
	      sci_parameters->error = 1;
	      Nbvars = nbvars_old;
	      return 0;
	    } /* End If */
	  
	  _SciErr = getVarAddressFromPosition(pvApiCtx, sci_parameters->ibegin, &tmp_1_addr); SCICOINOR_ERROR;
	  _SciErr = getVarType(pvApiCtx, tmp_1_addr, &var_type); SCICOINOR_ERROR;

	  if (var_type==sci_sparse)
	    {
	      Scierror(999,"scibonmin: the Hessian function must be not sparse\n");
	      sci_parameters->error = 1;
	      Nbvars = nbvars_old;
	      return 0;
	    }
	  
#ifdef DEBUG
	  DBGPRINTF("scibonmin: in eval_h - function called\n");
#endif
	  
	  // We get the resulting vector
	  _SciErr = getMatrixOfDouble(pvApiCtx, tmp_1_addr, &n_tmp_1, &m_tmp_1, &tmp_1); SCICOINOR_ERROR;
	  
	  if (n_tmp_1*m_tmp_1!=nele_hess)
	    {
	      Scierror(999,"scibonmin: the Hessian function must return %d values\n",nele_hess);
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
