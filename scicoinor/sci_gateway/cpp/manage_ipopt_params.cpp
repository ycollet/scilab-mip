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
#include <scilabjournal.hpp>

#include <string>
#include <new>
#include <string>

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

int manage_ipopt_params(Ipopt::SmartPtr<Ipopt::OptionsList> options, int * param_in_addr, int Log)
{
  // Get the parameters stored in the plist
  int     tmp_res, tmp_int;
  char *  tmp_char;
  double  tmp_double;

#ifdef DEBUG
  DBGPRINTF("sciipopt: processing options\n");
#endif

  ////////////
  // Output //
  ////////////
  // verbosity level
  getIntInPList(pvApiCtx, param_in_addr, "print_level", &tmp_int, &tmp_res, 5, Log, CHECK_BOTH, 0, 12);
  if (tmp_res!=-1) options->SetIntegerValue("print_level", tmp_int);  
  // File name of options file (default: ipopt.opt)
  getStringInPList(pvApiCtx, param_in_addr, "option_file_name", &tmp_char, &tmp_res, "ipopt.opt", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("option_file_name", tmp_char);
  FREE(tmp_char);

  /////////////////
  // Termination //
  /////////////////
  // Desired convergence tolerance (relative)
  getDoubleInPList(pvApiCtx, param_in_addr, "tol", &tmp_double, &tmp_res, 1e-8, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("tol", tmp_double);
  // Maximum number of iterations
  getIntInPList(pvApiCtx, param_in_addr, "max_iter", &tmp_int, &tmp_res, 3000, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetIntegerValue("max_iter", tmp_int);
  // Desired threshold for the dual infeasibility
  // Absolute tolerance on the dual infeasibility. Successful termination requires that the max-norm of the (unscaled) dual infeasibility is less than this threshold.
  getDoubleInPList(pvApiCtx, param_in_addr, "dual_inf_tol", &tmp_double, &tmp_res, 1, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("dual_inf_tol", tmp_double);
  // Desired threshold for the constraint violation
  // Absolute tolerance on the constraint violation. Successful termination requires that the max-norm of the (unscaled) constraint violation is less than this threshold.
  getDoubleInPList(pvApiCtx, param_in_addr, "constr_viol_tol", &tmp_double, &tmp_res, 1e-4, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("constr_viol_tol", tmp_double);
  // Acceptance threshold for the complementarity conditions
  // Absolute tolerance on the complementarity. Successful termination requires that the max-norm of the (unscaled) complementarity is less than this threshold.
  getDoubleInPList(pvApiCtx, param_in_addr, "compl_inf_tol", &tmp_double, &tmp_res, 1e-4, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("compl_inf_tol", tmp_double);
  // Acceptable convergence tolerance (relative)
  // Determines which (scaled) overall optimality error is considered to be "acceptable." 
  // There are two levels of termination criteria. If the usual "desired" tolerances (see tol, dual_inf_tol etc) are satisfied at an iteration, 
  // the algorithm immediately terminates with a success message. On the other hand, if the algorithm encounters "acceptable_iter" many iterations 
  // in a row that are considered "acceptable", it will terminate before the desired convergence tolerance is met. This is useful in cases where
  // the algorithm might not be able to achieve the "desired" level of accuracy.
  getDoubleInPList(pvApiCtx, param_in_addr, "acceptable_tol", &tmp_double, &tmp_res, 1e-6, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("acceptable_tol", tmp_double);
  // Acceptance threshold for the constraint violation
  // Absolute tolerance on the constraint violation. "Acceptable" termination requires that the max-norm of the (unscaled) constraint violation is
  // less than this threshold; see also acceptable_tol.
  getDoubleInPList(pvApiCtx, param_in_addr, "acceptable_constr_viol_tol", &tmp_double, &tmp_res, 1e-2, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("acceptable_constr_viol_tol", tmp_double);
  // Acceptance threshold for the dual infeasibility
  // Absolute tolerance on the dual infeasibility. "Acceptable" termination requires that the (max-norm of the unscaled) dual infeasibility is less than this 
  // threshold; see also acceptable_tol.
  getDoubleInPList(pvApiCtx, param_in_addr, "acceptable_dual_inf_tol", &tmp_double, &tmp_res, 1e10, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("acceptable_dual_inf_tol", tmp_double);
  // Acceptance threshold for the complementarity conditions
  // Absolute tolerance on the complementarity. "Acceptable" termination requires that the max-norm of the (unscaled) complementarity is less than this 
  // threshold; see also acceptable_tol.
  getDoubleInPList(pvApiCtx, param_in_addr, "acceptable_compl_inf_tol", &tmp_double, &tmp_res, 1e-2, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("acceptable_compl_inf_tol", tmp_double);
  // "Acceptance" stopping criterion based on objective function change.
  // If the relative change of the objective function (scaled by Max(1,|f(x)|)) is less than this value, this part of the acceptable tolerance termination 
  // is satisfied; see also acceptable_tol. This is useful for the quasi-Newton option, which has trouble to bring down the dual infeasibility.
  getDoubleInPList(pvApiCtx, param_in_addr, "acceptable_obj_change_tol", &tmp_double, &tmp_res, 1e20, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("acceptable_obj_change_tol", tmp_double);
  // Threshold for maximal value of primal iterates
  // If any component of the primal iterates exceeded this value (in absolute terms), the optimization is aborted with the exit message that the iterates
  // seem to be diverging.
  getDoubleInPList(pvApiCtx, param_in_addr, "diverging_iterates_tol", &tmp_double, &tmp_res, 1e20, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("diverging_iterates_tol", tmp_double);
  // Number of "acceptable" iterates before triggering termination.
  // If the algorithm encounters this many successive "acceptable" iterates (see "acceptable_tol"), it terminates, assuming that the problem 
  // has been solved to best possible accuracy given round-off. If it is set to zero, this heuristic is disabled.
  getIntInPList(pvApiCtx, param_in_addr, "acceptable_iter", &tmp_int, &tmp_res, 15, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetIntegerValue("acceptable_iter", tmp_int);
  // Indicates if all variable bounds should be replaced by inequality
  // constraints
  //   This option must be set for the inexact algorithm
  // Possible values:
  //  - no                      [leave bounds on variables]
  //  - yes                     [replace variable bounds by inequality constraints]
  getStringInPList(pvApiCtx, param_in_addr, "replace_bounds", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=1) options->SetStringValue("replace_bounds", tmp_char);
  FREE(tmp_char);
  // Scaling threshold for the NLP error.
  //   (See paragraph after Eqn. (6) in the implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "s_max", &tmp_double, &tmp_res, 100, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("s_max", tmp_double);
  // Maximum number of CPU seconds.
  //   A limit on CPU seconds that Ipopt can use to solve one problem.  If
  //   during the convergence check this limit is exceeded, Ipopt will terminate
  //   with a corresponding error message.
  getDoubleInPList(pvApiCtx, param_in_addr, "max_cpu_time", &tmp_double, &tmp_res, 1e6, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("max_cpu_time", tmp_double);

  /////////////////
  // NLP Scaling //
  /////////////////
  // Scaling factor for the objective function
  getDoubleInPList(pvApiCtx, param_in_addr, "obj_scaling_factor", &tmp_double, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("obj_scaling_factor", tmp_double);
  // Select the technique used for scaling the NLP
  //  none: no problem scaling will be performed
  //  user-scaling: scaling parameters will come from the user
  //  gradient-based: scale the problem so the maximum gradient at the starting point is scaling_max_gradient
  //  equilibration-based: scale the problem so that first derivatives are of order 1 at random points (only available with MC19)
  // Selects the technique used for scaling the problem internally before it is solved. For user-scaling, the parameters come from the NLP.
  //  If you are using AMPL, they can be specified through suffixes ("scaling_factor")
  getStringInPList(pvApiCtx, param_in_addr, "nlp_scaling_method", &tmp_char, &tmp_res, "gradient-based", Log, 
		       CHECK_VALUES, 4, 
		       "none",
		       "user_scaling",
		       "gradient-based",
		       "equilibration-based");
  if (tmp_res!=-1) options->SetStringValue("nlp_scaling_method", tmp_char);
  FREE(tmp_char);
  // Maximum gradient after scaling
  // This is the gradient scaling cut-off. If the maximum gradient is above this value, then gradient based scaling will be performed.
  // Scaling parameters are calculated to scale the maximum gradient back to this value. (This is g_max in Section 3.8 of the implementation paper.)
  //  Note: This option is only used if "nlp_scaling_method" is chosen as "gradient-based".
  getDoubleInPList(pvApiCtx, param_in_addr, "nlp_scaling_max_gradient", &tmp_double, &tmp_res, 100, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("nlp_scaling_max_gradient", tmp_double);
  // Target value for objective function gradient size.
  //   If a positive number is chosen, the scaling factor the objective function
  //   is computed so that the gradient has the max norm of the given size at
  //   the starting point.  This overrides nlp_scaling_max_gradient for the
  //   objective function.
  getDoubleInPList(pvApiCtx, param_in_addr, "nlp_scaling_obj_target_gradient", &tmp_double, &tmp_res, 0, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("nlp_scaling_obj_target_gradient", tmp_double);
  // Target value for constraint function gradient size.
  //   If a positive number is chosen, the scaling factor the constraint
  //   functions is computed so that the gradient has the max norm of the given
  //   size at the starting point.  This overrides nlp_scaling_max_gradient for
  //   the constraint functions.
  getDoubleInPList(pvApiCtx, param_in_addr, "nlp_scaling_constr_target_gradient", &tmp_double, &tmp_res, 0, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("nlp_scaling_constr_target_gradient", tmp_double);

  /////////
  // NLP //
  /////////
  // Factor for initial relaxation of the bounds
  // Before start of the optimization, the bounds given by the user are relaxed.  This option sets the factor for this relaxation. If it 
  // is set to zero, then then bounds relaxation is disabled. (See Eqn.(35) in implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "bound_relax_factor", &tmp_double, &tmp_res, 1e-8, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("bound_relax_factor", tmp_double);
  // Indicates whether final points should be projected into original bounds.
  // * no: Leave final point unchanged
  // * yes: Project final point back into original bounds
  // Ipopt might relax the bounds during the optimization (see, e.g., option "bound_relax_factor"). This option determines whether the final 
  // point should be projected back into the user-provide original bounds after the optimization.
  getStringInPList(pvApiCtx, param_in_addr, "honor_original_bounds", &tmp_char, &tmp_res, "yes", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("honor_original_bounds", tmp_char);
  FREE(tmp_char);
  // Indicates whether it is desired to check for Nan/Inf in derivative matrices
  // * no: Don't check (faster).
  // * yes: Check Jacobians and Hessian for Nan and Inf.
  // Activating this option will cause an error if an invalid number is detected in the constraint Jacobians or the Lagrangian Hessian. 
  // If this is not activated, the test is skipped, and the algorithm might proceed with invalid numbers and fail.
  getStringInPList(pvApiCtx, param_in_addr, "check_derivatives_for_naninf", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("check_derivatives_for_naninf", tmp_char);
  FREE(tmp_char);
  // any bound less or equal this value will be considered -inf (i.e. not lower bounded).
  getDoubleInPList(pvApiCtx, param_in_addr, "nlp_lower_bound_inf", &tmp_double, &tmp_res, -1e19, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("nlp_lower_bound_inf", tmp_double);
  // any bound greater or this value will be considered +inf (i.e. not upper bounded).
  getDoubleInPList(pvApiCtx, param_in_addr, "nlp_upper_bound_inf", &tmp_double, &tmp_res, 1e19, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("nlp_upper_bound_inf", tmp_double);
  // Determines how fixed variables should be handled.
  // * make_parameter: Remove fixed variable from optimization variables
  // * make_constraint: Add equality constraints fixing variables
  // * relax_bounds: Relax fixing bound constraints
  // The main difference between those options is that the starting point in the "make_constraint" case still has the fixed variables at 
  // their given values, whereas in the case "make_parameter" the functions are always evaluated with the fixed values for those variables.  
  // Also, for "relax_bounds", the fixing bound constraints are relaxed (according to "bound_relax_factor"). For both "make_constraints"
  //  and "relax_bounds", bound multipliers are computed for the fixed variables.
  getStringInPList(pvApiCtx, param_in_addr, "fixed_variable_treatment", &tmp_char, &tmp_res, "make_parameter", Log, 
		       CHECK_VALUES, 3, 
		       "make_parameter",
		       "make_constraint",
		       "relax_bounds");
  if (tmp_res!=-1) options->SetStringValue("fixed_variable_treatment", tmp_char);
  FREE(tmp_char);
  //Indicates which linear solver should be used to detect linearly dependent equality constraints.
  // none: don't check; no extra work at beginning
  // mumps: use MUMPS
  // wsmp: use WSMP
  // ma28: use MA28
  // The default and available choices depend on how Ipopt has been compiled. This is experimental and does not work well.
  getStringInPList(pvApiCtx, param_in_addr, "dependency_detector", &tmp_char, &tmp_res, "none", Log, 
		       CHECK_VALUES, 4, 
		       "none",
		       "mumps",
		       "wsmp",
		       "ma28");
  if (tmp_res!=-1) options->SetStringValue("dependency_detector", tmp_char);
  FREE(tmp_char);
  // Indicates if the right hand sides of the constraints should be considered during dependency detection
  // no: only look at gradients
  // yes: also consider right hand side
  getStringInPList(pvApiCtx, param_in_addr, "dependency_detection_with_rhs", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("dependency_detection_with_rhs", tmp_char);
  FREE(tmp_char);
  // Number of linear variables
  // When the Hessian is approximated, it is assumed that the first num_linear_variables variables are linear. The Hessian is then not 
  // approximated in this space. If the get_number_of_nonlinear_variables method in the TNLP is implemented, this option is ignored.
  getIntInPList(pvApiCtx, param_in_addr, "num_linear_variables", &tmp_int, &tmp_res, 0, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetIntegerValue("num_linear_variables", tmp_int);
  // Weight for linear damping term (to handle one-sided bounds).
  //   (see Section 3.7 in implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "kappa_d", &tmp_double, &tmp_res, 1e-5, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("kappa_d", tmp_double);
  // Indicates whether all equality constraints are linear
  // * no: Don't assume that all equality constraints are linear
  // * yes: Assume that equality constraints Jacobian are constant
  // Activating this option will cause Ipopt to ask for the Jacobian of the equality constraints only once from the NLP and reuse this information later.
  getStringInPList(pvApiCtx, param_in_addr, "jac_c_constant", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("jac_c_constant", tmp_char);
  FREE(tmp_char);
  // Indicates whether all inequality constraints are linear
  // * no: Don't assume that all inequality constraints are linear
  // * yes: Assume that equality constraints Jacobian are constant
  // Activating this option will cause Ipopt to ask for the Jacobian of the inequality constraints only once from the NLP and reuse this information later.
  getStringInPList(pvApiCtx, param_in_addr, "jac_d_constant", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("jac_d_constant", tmp_char);
  FREE(tmp_char);
  // Indicates whether the problem is a quadratic problem
  // * no: Assume that Hessian changes
  // * yes: Assume that Hessian is constant
  // Activating this option will cause Ipopt to ask for the Hessian of the Lagrangian function only once from the NLP and reuse this information later.
  getStringInPList(pvApiCtx, param_in_addr, "hessian_constant", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("hessian_constant", tmp_char);
  FREE(tmp_char);

  ////////////////////
  // Initialization //
  ////////////////////
  // Desired minimal relative distance of initial point to bound
  // Determines how much the initial point might have to be modified in order to be sufficiently inside the bounds (together with "bound_push").
  // (This is kappa_2 in Section 3.6 of implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "bound_frac", &tmp_double, &tmp_res, 1e-2, Log, CHECK_BOTH, 0, 0.5);
  if (tmp_res!=-1) options->SetNumericValue("bound_frac", tmp_double);
  // Desired minimal absolute distance of initial point to bound
  // Determines how much the initial point might have to be modified in order to be sufficiently inside the bounds (together with "bound_frac").
  // (This is kappa_1 in "Section 3.6 of implementation paper.
  getDoubleInPList(pvApiCtx, param_in_addr, "bound_push", &tmp_double, &tmp_res, 1e-2, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("bound_push", tmp_double);
  // Desired minimal relative distance of initial slack to bound
  // Determines how much the initial slack variables might have to be modified in order to be sufficiently inside the inequality bounds
  // (together with "slack_bound_push").  (This is kappa_2 in Section 3.6 of implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "slack_bound_frac", &tmp_double, &tmp_res, 1e-2, Log, CHECK_BOTH, 0, 0.5);
  if (tmp_res!=-1) options->SetNumericValue("slack_bound_frac", tmp_double);
  // Desired minimal absolute distance of initial slack to bound
  // Determines how much the initial slack variables might have to be modified in order to be sufficiently inside the inequality bounds
  //  (together with "slack_bound_frac").  (This is kappa_1 in Section 3.6 of implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "slack_bound_push", &tmp_double, &tmp_res, 1e-2, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("slack_bound_push", tmp_double);
  // Initial value for the bound multipliers
  // All dual variables corresponding to bound constraints are initialized to this value.
  getDoubleInPList(pvApiCtx, param_in_addr, "bound_mult_init_val", &tmp_double, &tmp_res, 1, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("bound_mult_init_val", tmp_double);
  // Maximal allowed least-square guess of constraint multipliers
  // Determines how large the initial least-square guesses of the constraint multipliers are allowed to be (in max-norm). If the guess is larger
  // than this value, it is discarded and all constraint multipliers are set to zero. This options is also used when initializing the 
  // restoration phase. By default, "resto.constr_mult_init_max" (the one used in RestoIterateInitializer) is set to zero.
  getDoubleInPList(pvApiCtx, param_in_addr, "constr_mult_init_max", &tmp_double, &tmp_res, 1000, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("constr_mult_init_max", tmp_double);
  // Initialization method for bound multipliers
  // * constant: set all bound multipliers to the value of bound_mult_init_val
  // * mu-based: initialize to mu_init/x_slack
  // This option defines how the iterates for the bound multipliers are initialized.  If "constant" is chosen, then all bound multipliers 
  // are initialized to the value of "bound_mult_init_val".  If "mu-based" is chosen, the each value is initialized to the the value 
  // of "mu_init" divided by the corresponding slack variable. This latter option might be useful if the starting point is close to the
  // optimal solution.
  getStringInPList(pvApiCtx, param_in_addr, "bound_mult_init_method", &tmp_char, &tmp_res, "constant", Log, 
		       CHECK_VALUES, 2, 
		       "constant",
		       "mu-based");
  if (tmp_res!=-1) options->SetStringValue("bound_mult_init_method", tmp_char);
  FREE(tmp_char);
  // Least square initialization of the primal variables
  // no: take user-provided point
  // yes: overwrite user-provided point with least-square estimates
  // If set to yes, Ipopt ignores the user provided point and solves a least square problem for the primal variables (x and s), to fit the 
  // linearized equality and inequality constraints.  This might be useful if the user doesn't know anything about the starting point, or for 
  // solving an LP or QP.
  getStringInPList(pvApiCtx, param_in_addr, "least_square_init_primal", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("least_square_init_primal", tmp_char);
  FREE(tmp_char);
  // Least square initialization of all dual variables
  // no: use bound_mult_init_val and least-square equality constraint multipliers
  // yes: overwrite user-provided point with least-square estimates
  // If set to yes, Ipopt tries to compute least-square multipliers (considering ALL dual variables). If successful, the bound 
  // multipliers are possibly corrected to be at least bound_mult_init_val. This might be useful if the user doesn't know anything 
  // about the starting point, or for solving an LP or QP. This overwrites option "bound_mult_init_method".
  getStringInPList(pvApiCtx, param_in_addr, "least_square_init_duals", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("least_square_init_duals", tmp_char);
  FREE(tmp_char);

  ///////////////////////
  // Barrier Parameter //
  ///////////////////////
  // Indicates if we want to do Mehrotra's algorithm.
  // * no: Do the usual Ipopt algorithm.
  // * yes: Do Mehrotra's predictor-corrector algorithm.
  // If set to yes, Ipopt runs as Mehrotra's predictor-corrector algorithm. This works usually very well for LPs and convex QPs. This
  // automatically disables the line search, and chooses the (unglobalized) adaptive mu strategy with the "probing" oracle, and uses 
  // "corrector_type=affine" without any safeguards; you should not set any of those options explicitly in addition. Also, unless 
  // otherwise specified, the values of "bound_push", "bound_frac", and "bound_mult_init_val" are set more aggressive, and sets 
  // "alpha_for_y=bound_mult".
  getStringInPList(pvApiCtx, param_in_addr, "mehrotra_algorithm", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("mehrotra_algorithm", tmp_char);
  FREE(tmp_char);
  // Indicates if the linear system should be solved quickly.
  //   If set to yes, the algorithm assumes that the linear system that is
  //   solved to obtain the search direction, is solved sufficiently well. In
  //   that case, no residuals are computed, and the computation of the search
  //   direction is a little faster.
  // Possible values:
  //  - no                      [Verify solution of linear system by computing residuals.]
  //  - yes                     [Trust that linear systems are solved well.]
  getStringInPList(pvApiCtx, param_in_addr, "fast_step_computation", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("fast_step_computation", tmp_char);
  FREE(tmp_char);
  // Update strategy for barrier parameter
  // * monotone: use the monotone (Fiacco-McCormick) strategy
  // * adaptive: use the adaptive update strategy
  // Determines which barrier parameter update strategy is to be used.
  getStringInPList(pvApiCtx, param_in_addr, "mu_strategy", &tmp_char, &tmp_res, "monotone", Log, 
		       CHECK_VALUES, 2,
		       "monotone",
		       "adaptive");
  if (tmp_res!=-1) options->SetStringValue("mu_strategy", tmp_char);
  FREE(tmp_char);
  // Oracle for a new barrier parameter in the adaptive strategy
  // * probing: Mehrotra's probing heuristic
  // * loqo: LOQO's centrality rule
  // * quality-function: minimize a quality function
  // Determines how a new barrier parameter is computed in each "free-mode" iteration of the adaptive barrier parameter strategy.
  // (Only considered if "adaptive" is selected for option "mu_strategy").
  getStringInPList(pvApiCtx, param_in_addr, "mu_oracle", &tmp_char, &tmp_res, "quality-function", Log, 
		       CHECK_VALUES, 3,
		       "probing",
		       "loqo",
		       "quality-function");
  if (tmp_res!=-1) options->SetStringValue("mu_oracle", tmp_char);
  FREE(tmp_char);
  getIntInPList(pvApiCtx, param_in_addr, "quality_function_max_section_steps", &tmp_int, &tmp_res, 8, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetIntegerValue("quality_function_max_section_steps", tmp_int);
  // Oracle for the barrier parameter when switching to fixed mode.
  // * probing: Mehrotra's probing heuristic
  // * loqo: LOQO's centrality rule
  // * quality-function: minimize a quality function
  // * average_compl: base on current average complementarity
  // Determines how the first value of the barrier parameter should be computed when switching to the "monotone mode" in the adaptive strategy. 
  // (Only considered if "adaptive" is selected for option "mu_strategy")
  getStringInPList(pvApiCtx, param_in_addr, "fixed_mu_oracle", &tmp_char, &tmp_res, "average_compl", Log, 
		       CHECK_VALUES, 4,
		       "probing",
		       "loqo",
		       "quality-function",
		       "average_compl");
  if (tmp_res!=-1) options->SetStringValue("fixed_mu_oracle", tmp_char);
  FREE(tmp_char);
  // Initial value for the barrier parameter
  // This option determines the initial value for the barrier parameter (mu).It is only relevant in the monotone, Fiacco-McCormick 
  // version of the algorithm. (i.e., if "mu_strategy" is chosen as "monotone")
  getDoubleInPList(pvApiCtx, param_in_addr, "mu_init", &tmp_double, &tmp_res, 0.1, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("mu_init", tmp_double);
  // Factor for initialization of maximum value for barrier parameter.
  // This option determines the upper bound on the barrier parameter. This upper bound is computed as the average complementarity at the initial 
  // point times the value of this option. (Only used if option "mu_strategy" is chosen as "adaptive".)
  getDoubleInPList(pvApiCtx, param_in_addr, "mu_max_fact", &tmp_double, &tmp_res, 1000, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("mu_max_fact", tmp_double);
  // Maximal value for barrier parameter for adaptive strategy
  // This option specifies an upper bound on the barrier parameter in the adaptive mu selection mode. If this option is set, it overwrites the 
  // effect of mu_max_fact. (Only used if option "mu_strategy" is chosen as "adaptive".)
  getDoubleInPList(pvApiCtx, param_in_addr, "mu_max", &tmp_double, &tmp_res, 100000, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("mu_max", tmp_double);
  // Minimum value for barrier parameter.
  // This option specifies the lower bound on the barrier parameter in the adaptive mu selection mode. By default, it is set to the minimum of 1e-11 and 
  // min("tol","compl_inf_tol")/("barrier_tol_factor"+1), which should be a reasonable value. (Only used if option "mu_strategy" is chosen as "adaptive".)
  getDoubleInPList(pvApiCtx, param_in_addr, "mu_min", &tmp_double, &tmp_res, 1e-11, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("mu_min", tmp_double);
  // Factor for mu in barrier stop test.
  // The convergence tolerance for each barrier problem in the monotone mode is the value of the barrier parameter times "barrier_tol_factor".
  // This option is also used in the adaptive mu strategy during the monotone mode. (This is kappa_epsilon in implementation paper).
  getDoubleInPList(pvApiCtx, param_in_addr, "barrier_tol_factor", &tmp_double, &tmp_res, 10, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("barrier_tol_factor", tmp_double);
  // Determines linear decrease rate of barrier parameter.
  // For the Fiacco-McCormick update procedure the new barrier parameter mu is obtained by taking the minimum of mu*"mu_linear_decrease_factor"
  // and mu^"superlinear_decrease_power". (This is kappa_mu in implementation paper.) This option is also used in the adaptive mu 
  // strategy during the monotone mode.
  getDoubleInPList(pvApiCtx, param_in_addr, "mu_linear_decrease_factor", &tmp_double, &tmp_res, 0.2, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("mu_linear_decrease_factor", tmp_double);
  // Determines superlinear decrease rate of barrier parameter.
  // For the Fiacco-McCormick update procedure the new barrier parameter mu is obtained by taking the minimum of mu*"mu_linear_decrease_factor"
  // and mu^"superlinear_decrease_power". (This is theta_mu in implementation paper.) This option is also used in the adaptive mu 
  // strategy during the monotone mode.
  getDoubleInPList(pvApiCtx, param_in_addr, "mu_superlinear_decrease_power", &tmp_double, &tmp_res, 1.5, Log, CHECK_BOTH, 1, 2);
  if (tmp_res!=-1) options->SetNumericValue("mu_superlinear_decrease_power", tmp_double);
  // Allow skipping of barrier problem if barrier test is already met.
  // no: Take at least one iteration per barrier problem
  // yes: Allow fast decrease of mu if barrier test it met
  // If set to "no", the algorithm enforces at least one iteration per barrier problem, even if the barrier test is already met for the updated barrier parameter.
  getStringInPList(pvApiCtx, param_in_addr, "mu_allow_fast_monotone_decrease", &tmp_char, &tmp_res, "yes", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("mu_allow_fast_monotone_decrease", tmp_char);
  FREE(tmp_char);
  // Lower bound on fraction-to-the-boundary parameter tau
  // (This is tau_min in the implementation paper.)  This option is also used in the adaptive mu strategy during the monotone mode.
  getDoubleInPList(pvApiCtx, param_in_addr, "tau_min", &tmp_double, &tmp_res, 0.99, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("tau_min", tmp_double);
  // Maximum value of the centering parameter.
  //   This is the upper bound for the centering parameter chosen by the quality
  //   function based barrier parameter update. (Only used if option "mu_oracle"
  //   is set to "quality-function".)
  getDoubleInPList(pvApiCtx, param_in_addr, "sigma_max", &tmp_double, &tmp_res, 100, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("sigma_max", tmp_double);
  // Minimum value of the centering parameter.
  //   This is the lower bound for the centering parameter chosen by the quality
  //   function based barrier parameter update. (Only used if option "mu_oracle"
  //   is set to "quality-function".)
  getDoubleInPList(pvApiCtx, param_in_addr, "sigma_min", &tmp_double, &tmp_res, 100, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("sigma_min", tmp_double);
  // Norm used for components of the quality function.
  //   (Only used if option "mu_oracle" is set to "quality-function".)
  // Possible values:
  //  - 1-norm                  [use the 1-norm (abs sum)]
  //  - 2-norm-squared          [use the 2-norm squared (sum of squares)]
  //  - max-norm                [use the infinity norm (max)]
  //  - 2-norm                  [use 2-norm]
  getStringInPList(pvApiCtx, param_in_addr, "quality_function_norm_type", &tmp_char, &tmp_res, "2-norm-squared", Log, 
		       CHECK_VALUES, 4, 
		       "1-norm", 
		       "2-norm-squared",
		       "max-norm",
		       "2-norm");
  if (tmp_res!=-1) options->SetStringValue("quality_function_norm_type", tmp_char);
  FREE(tmp_char);
  // The penalty term for centrality that is included in quality function.
  //   This determines whether a term is added to the quality function to
  //   penalize deviation from centrality with respect to complementarity.  The
  //   complementarity measure here is the xi in the Loqo update rule. (Only
  //   used if option "mu_oracle" is set to "quality-function".)
  // Possible values:
  //  - none                    [no penalty term is added]
  //  - log                     [complementarity * the log of the centrality measure]
  //  - reciprocal              [complementarity * the reciprocal of the centrality measure]
  //  - cubed-reciprocal        [complementarity * the reciprocal of the centrality measure cubed]
  getStringInPList(pvApiCtx, param_in_addr, "quality_function_centrality", &tmp_char, &tmp_res, "none", Log, 
		       CHECK_VALUES, 4, 
		       "none", 
		       "log",
		       "reciprocal",
		       "cubed-reciprocal");
  if (tmp_res!=-1) options->SetStringValue("quality_function_centrality", tmp_char);
  FREE(tmp_char);
  // The balancing term included in the quality function for centrality.
  //   This determines whether a term is added to the quality function that
  //   penalizes situations where the complementarity is much smaller than dual
  //   and primal infeasibilities. (Only used if option "mu_oracle" is set to
  //   "quality-function".)
  // Possible values:
  //  - none                    [no balancing term is added]
  //  - cubic                   [Max(0,Max(dual_inf,primal_inf)-compl)^3]
  getStringInPList(pvApiCtx, param_in_addr, "quality_function_balancing_term", &tmp_char, &tmp_res, "none", Log, 
		       CHECK_VALUES, 2, 
		       "none", 
		       "cubic");
  if (tmp_res!=-1) options->SetStringValue("quality_function_balancing_term", tmp_char);
  FREE(tmp_char);
  // Tolerance for the section search procedure determining the optimal
  // centering parameter (in sigma space).
  //   The golden section search is performed for the quality function based mu
  //   oracle. (Only used if option "mu_oracle" is set to "quality-function".)
  getDoubleInPList(pvApiCtx, param_in_addr, "quality_function_section_sigma_tol", &tmp_double, &tmp_res, 0.01, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("quality_function_section_sigma_tol", tmp_double);
  // Tolerance for the golden section search procedure determining the optimal
  // centering parameter (in the function value space).
  //   The golden section search is performed for the quality function based mu
  //   oracle. (Only used if option "mu_oracle" is set to "quality-function".)
  getDoubleInPList(pvApiCtx, param_in_addr, "quality_function_section_qf_tol", &tmp_double, &tmp_res, 0, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("quality_function_section_qf_tol", tmp_double);
  // Factor limiting the deviation of dual variables from primal estimates.
  // If the dual variables deviate from their primal estimates, a correction is performed. (See Eqn. (16) in the implementation paper.) 
  // Setting the value to less than 1 disables the correction.
  getDoubleInPList(pvApiCtx, param_in_addr, "kappa_sigma", &tmp_double, &tmp_res, 1e10, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("kappa_sigma", tmp_double);
  // Globalization method used in backtracking line search
  // filter: Filter method
  // cg-penalty: Chen-Goldfarb penalty function
  // penalty: Standard penalty function
  getStringInPList(pvApiCtx, param_in_addr, "line_search_method", &tmp_char, &tmp_res, "cg-penalty", Log, 
		       CHECK_VALUES, 3, 
		       "filter",
		       "cg-penalty",
		       "penalty");
  if (tmp_res!=-1) options->SetStringValue("line_search_method", tmp_char);
  FREE(tmp_char);
  // Globalization strategy for the adaptive mu selection mode.
  // kkt-error: nonmonotone decrease of kkt-error
  // obj-constr-filter: 2-dim filter for objective and constraint violation
  // never-monotone-mode: disables globalization
  // To achieve global convergence of the adaptive version, the algorithm has to switch to the monotone mode (Fiacco-McCormick approach) when 
  // convergence does not seem to appear.  This option sets the criterion used to decide when to do this switch. (Only used if option 
  // "mu_strategy" is chosen as "adaptive".)
  getStringInPList(pvApiCtx, param_in_addr, "adaptive_mu_globalization", &tmp_char, &tmp_res, "obj-constr-filter", Log, 
		       CHECK_VALUES, 3, 
		       "kkt-error",
		       "obj_constr_filter",
		       "never_monotone_mode");
  if (tmp_res!=-1) options->SetStringValue("adaptive_mu_globalization", tmp_char);
  FREE(tmp_char);
  // Maximum number of iterations requiring sufficient progress.
  // For the "kkt-error" based globalization strategy, sufficient progress must be made for "adaptive_mu_kkterror_red_iters
  // iterations. If this number of iterations is exceeded, the globalization strategy switches to the monotone mode.
  getIntInPList(pvApiCtx, param_in_addr, "adaptive_mu_kkterror_red_iters", &tmp_int, &tmp_res, 4, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetIntegerValue("adaptive_mu_kkterror_red_iters", tmp_int);
  // Sufficient decrease factor for "kkt-error" globalization strategy.
  // For the "kkt-error" based globalization strategy, the error must decrease by this factor to be deemed sufficient decrease.
  getDoubleInPList(pvApiCtx, param_in_addr, "adaptive_mu_kkterror_red_fact", &tmp_double, &tmp_res, 0.9999, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("adaptive_mu_kkterror_red_fact", tmp_double);
  // Factor determining width of margin for obj-constr-filter adaptive globalization strategy.
  // When using the adaptive globalization strategy, "obj-constr-filter" sufficient progress for a filter entry is defined as 
  // follows: (new obj) < (filter obj) - filter_margin_fact*(new constr-viol) OR (new constr-viol) < (filter constr-viol) - 
  // filter_margin_fact*(new constr-viol).  For the description of the "kkt-error-filter" option see "filter_max_margin.
  getDoubleInPList(pvApiCtx, param_in_addr, "filter_margin_fact", &tmp_double, &tmp_res, 1e-5, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("filter_margin_fact", tmp_double);
  // Maximum width of margin in obj-constr-filter adaptive globalization strategy.
  getDoubleInPList(pvApiCtx, param_in_addr, "filter_max_margin", &tmp_double, &tmp_res, 1.0, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("filter_max_margin", tmp_double);
  // Indicates if the previous iterate should be restored if the monotone mode is entered.
  // no: don't restore accepted iterate
  // yes: restore accepted iterate
  // When the globalization strategy for the adaptive barrier algorithm switches to the monotone mode, it can either start 
  // from the most recent iterate (no), or from the last iterate that was accepted (yes).
  getStringInPList(pvApiCtx, param_in_addr, "adaptive_mu_restore_previous_iterate", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("adaptive_mu_restore_previous_iterate", tmp_char);
  FREE(tmp_char);
  // Determines the initial value of the barrier parameter when switching to the monotone mode.
  // When the globalization strategy for the adaptive barrier algorithm switches to the monotone mode and fixed_mu_oracle is chosen as 
  // "average_compl", the barrier parameter is set to the current average complementarity times the value of "adaptive_mu_monotone_init_factor"
  getDoubleInPList(pvApiCtx, param_in_addr, "adaptive_mu_monotone_init_factor", &tmp_double, &tmp_res, 0.8, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("adaptive_mu_monotone_init_factor", tmp_double);
  // Norm used for the KKT error in the adaptive mu globalization strategies.
  // 1-norm: use the 1-norm (abs sum)
  // 2-norm-squared: use the 2-norm squared (sum of squares)
  // max-norm: use the infinity norm (max)
  // 2-norm: use 2-norm
  // When computing the KKT error for the globalization strategies, the norm to be used is specified with this option. Note, this options is also used 
  // in the QualityFunctionMuOracle.
  getStringInPList(pvApiCtx, param_in_addr, "adaptive_mu_kkt_norm_type", &tmp_char, &tmp_res, "2-norm-squared", Log, 
		       CHECK_VALUES, 4, 
		       "1-norm",
		       "2-norm-squared",
		       "max-norm",
		       "2-norm");
  if (tmp_res!=-1) options->SetStringValue("adaptive_mu_kkt_norm_type", tmp_char);
  FREE(tmp_char);

  ////////////////////////
  // Multiplier Updates //
  ////////////////////////

  // Fractional reduction of the trial step size in the backtracking line search.
  //   At every step of the backtracking line search, the trial step size is
  //   reduced by this factor.
  getDoubleInPList(pvApiCtx, param_in_addr, "alpha_red_factor", &tmp_double, &tmp_res, 0.5, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("alpha_red_factor", tmp_double);
  // Always accept the first trial step.
  //   Setting this option to "yes" essentially disables the line search and
  //   makes the algorithm take aggressive steps, without global convergence guarantees.
  // Possible values:
  //  - no                      [don't arbitrarily accept the full step]
  //  - yes                     [always accept the full step]
  getStringInPList(pvApiCtx, param_in_addr, "accept_every_trial_step", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("accept_every_trial_step", tmp_char);
  FREE(tmp_char);
  // Accept a trial point after maximal this number of steps.
  //   Even if it does not satisfy line search conditions.
  getDoubleInPList(pvApiCtx, param_in_addr, "accept_after_max_steps", &tmp_double, &tmp_res, -1, Log, CHECK_MIN, -1);
  if (tmp_res!=-1) options->SetNumericValue("accept_after_max_steps", tmp_double);

  // Method to determine the step size for constraint multipliers.
  // * primal: use primal step size
  // * bound_mult: use step size for the bound multipliers (good for LPs)
  // * min: use the min of primal and bound multipliers
  // * max: use the max of primal and bound multipliers
  // * full: take a full step of size one
  // * min_dual_infeas: choose step size minimizing new dual infeasibility
  // * safe_min_dual_infeas: like "min_dual_infeas", but safeguarded by "min" and "max"
  // * primal-and-full: use the primal step size, and full step if delta_x <= alpha_for_y_tol
  // * dual-and-full: use the dual step size, and full step if delta_x <= alpha_for_y_tol
  // * acceptor: Call LSAcceptor to get step size for y
  // This option determines how the step size (alpha_y) will be calculated when updating the constraint multipliers.
  getStringInPList(pvApiCtx, param_in_addr, "alpha_for_y", &tmp_char, &tmp_res, "primal", Log, 
		       CHECK_VALUES, 10, 
		       "primal",
		       "bound-mult",
		       "min",
		       "max",
		       "full",
		       "min-dual-infeas",
		       "safer-min-dual-infeas",
		       "primal-and-full",
		       "dual-and-full",
		       "acceptor");
  if (tmp_res!=-1) options->SetStringValue("alpha_for_y", tmp_char);
  FREE(tmp_char);
  // Tolerance for switching to full equality multiplier steps
  // This is only relevant if "alpha_for_y" is chosen "primal-and-full" or "dual-and-full". The step size for the equality constraint 
  // multipliers is taken to be one if the max-norm of the primal step is less than this tolerance
  getDoubleInPList(pvApiCtx, param_in_addr, "alpha_for_y_tol", &tmp_double, &tmp_res, 10.0, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("alpha_for_y_tol", tmp_double);
  // Tolerance for detecting numerically insignificant steps
  // If the search direction in the primal variables (x and s) is, in relative terms for each component, less than this value, the 
  // algorithm accepts the full step without line search. If this happens repeatedly, the algorithm will terminate with a corresponding exit 
  // message. The default value is 10 times machine precision.  
  getDoubleInPList(pvApiCtx, param_in_addr, "tiny_step_tol", &tmp_double, &tmp_res, 10.0*std::numeric_limits<double>::epsilon(), Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("tiny_step_tol", tmp_double);
  // Tolerance for quitting because of numerically insignificant steps.
  // If the search direction in the primal variables (x and s) is, in relative terms for each component, repeatedly less than tiny_step_tol,
  // and the step in the y variables is smaller than this threshold, the algorithm will terminate.
  getDoubleInPList(pvApiCtx, param_in_addr, "tiny_step_y_tol", &tmp_double, &tmp_res, 1e-2, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("tiny_step_y_tol", tmp_double);
  // Tells the algorithm to recalculate the equality and inequality multipliers as least square estimates.
  // no: use the Newton step to update the multipliers
  // yes: use least-square multiplier estimates
  // This asks the algorithm to recompute the multipliers, whenever the current infeasibility is less than recalc_y_feas_tol. 
  // Choosing yes might be helpful in the quasi-Newton option.  However, each recalculation requires an extra factorization of the linear 
  // system.  If a limited memory quasi-Newton option is chosen, this is used by default.
  getStringInPList(pvApiCtx, param_in_addr, "recalc_y", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("recalc_y", tmp_char);
  FREE(tmp_char);
  // Feasibility threshold for recomputation of multipliers.
  // If recalc_y is chosen and the current infeasibility is less than this value, then the multipliers are recomputed.
  getDoubleInPList(pvApiCtx, param_in_addr, "recalc_y_feas_tol", &tmp_double, &tmp_res, 1e-6, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("recalc_y_feas_tol", tmp_double);
  // Correction size for very small slacks.
  //   Due to numerical issues or the lack of an interior, the slack variables
  //   might become very small.  If a slack becomes very small compared to
  //   machine precision, the corresponding bound is moved slightly.  This
  //   parameter determines how large the move should be.  Its default value is
  //   mach_eps^{3/4}.  (See also end of Section 3.5 in implementation paper -
  //   but actual implementation might be somewhat different.)  
  getDoubleInPList(pvApiCtx, param_in_addr, "slack_move", &tmp_double, &tmp_res, 1.81899e-12, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("slack_move", tmp_double);

  /////////////////
  // Line Search //
  /////////////////
  // Maximal number of second order correction trial steps
  getIntInPList(pvApiCtx, param_in_addr, "max_soc", &tmp_int, &tmp_res, 4, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetIntegerValue("max_soc", tmp_int);
  // Factor in the sufficient reduction rule for second order correction.
  //   This option determines how much a second order correction step must
  //   reduce the constraint violation so that further correction steps are
  //   attempted.  (See Step A-5.9 of Algorithm A in the implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "kappa_soc", &tmp_double, &tmp_res, 0.99, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("kappa_soc", tmp_double);
  // Determines the upper bound on the acceptable increase of barrier objective function.
  //   Trial points are rejected if they lead to an increase in the barrier
  //   objective function by more than obj_max_inc orders of magnitude.
  getDoubleInPList(pvApiCtx, param_in_addr, "obj_max_inc", &tmp_double, &tmp_res, 5, Log, CHECK_MIN, 1);
  if (tmp_res!=-1) options->SetNumericValue("obj_max_inc", tmp_double);
  // Maximal allowed number of filter resets
  //   A positive number enables a heuristic that resets the filter, whenever in
  //   more than "filter_reset_trigger" successive iterations the last rejected
  //   trial steps size was rejected because of the filter.  This option
  //   determine the maximal number of resets that are allowed to take place.
  getIntInPList(pvApiCtx, param_in_addr, "max_filter_resets", &tmp_int, &tmp_res, 5, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetIntegerValue("max_filter_resets", tmp_int);
  // Number of iterations that trigger the filter reset.
  //   If the filter reset heuristic is active and the number of successive
  //   iterations in which the last rejected trial step size was rejected
  //   because of the filter, the filter is reset.
  getIntInPList(pvApiCtx, param_in_addr, "filter_reset_trigger", &tmp_int, &tmp_res, 5, Log, CHECK_MIN, 1);
  if (tmp_res!=-1) options->SetIntegerValue("filter_reset_trigger", tmp_int);
  // Trigger counter for watchdog procedure
  // If the number of successive iterations in which the backtracking line search did not accept the first trial point exceeds this number, the 
  // watchdog procedure is activated.  Choosing "0" here disables the watchdog procedure.
  getIntInPList(pvApiCtx, param_in_addr, "watchdog_shortened_iter_trigger", &tmp_int, &tmp_res, 10, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetIntegerValue("watchdog_shortened_iter_trigger", tmp_int);
  // Maximum number of watchdog iterations.
  // This option determines the number of trial iterations allowed before the watchdog procedure is aborted and the algorithm returns to the stored point.
  getIntInPList(pvApiCtx, param_in_addr, "watchdog_trial_iter_max", &tmp_int, &tmp_res, 3, Log, CHECK_MIN, 1);
  if (tmp_res!=-1) options->SetIntegerValue("watchdog_trial_iter_max", tmp_int);
  // Determines upper bound for constraint violation in the filter.
  //   The algorithmic parameter theta_max is determined as theta_max_fact times
  //   the maximum of 1 and the constraint violation at initial point.  Any
  //   point with a constraint violation larger than theta_max is unacceptable
  //   to the filter (see Eqn. (21) in the implementation paper).
  getDoubleInPList(pvApiCtx, param_in_addr, "theta_max_fact", &tmp_double, &tmp_res, 10000, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("theta_max_fact", tmp_double);
  // Determines constraint violation threshold in the switching rule.
  //   The algorithmic parameter theta_min is determined as theta_min_fact times
  //   the maximum of 1 and the constraint violation at initial point.  The
  //   switching rules treats an iteration as an h-type iteration whenever the
  //   current constraint violation is larger than theta_min (see paragraph
  //   before Eqn. (19) in the implementation paper).
  getDoubleInPList(pvApiCtx, param_in_addr, "theta_min_fact", &tmp_double, &tmp_res, 0.0001, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("theta_min_fact", tmp_double);
  // Relaxation factor in the Armijo condition.
  //   (See Eqn. (20) in the implementation paper)
  getDoubleInPList(pvApiCtx, param_in_addr, "eta_phi", &tmp_double, &tmp_res, 1e-8, Log, CHECK_BOTH, 0, 0.5);
  if (tmp_res!=-1) options->SetNumericValue("eta_phi", tmp_double);
  // Multiplier for constraint violation in the switching rule.
  //   (See Eqn. (19) in the implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "delta", &tmp_double, &tmp_res, 1, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("delta", tmp_double);
  // Exponent for linear barrier function model in the switching rule.
  //   (See Eqn. (19) in the implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "s_phi", &tmp_double, &tmp_res, 2.3, Log, CHECK_MIN, 1);
  if (tmp_res!=-1) options->SetNumericValue("s_phi", tmp_double);
  // Exponent for current constraint violation in the switching rule.
  //   (See Eqn. (19) in the implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "s_theta", &tmp_double, &tmp_res, 1.1, Log, CHECK_MIN, 1);
  if (tmp_res!=-1) options->SetNumericValue("s_theta", tmp_double);
  // Relaxation factor in the filter margin for the barrier function.
  //   (See Eqn. (18a) in the implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "gamma_phi", &tmp_double, &tmp_res, 1e-8, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("gamma_phi", tmp_double);
  // Relaxation factor in the filter margin for the constraint violation.
  //   (See Eqn. (18b) in the implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "gamma_theta", &tmp_double, &tmp_res, 1e-5, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("gamma_theta", tmp_double);
  // Safety factor for the minimal step size (before switching to restoration phase).
  //   (This is gamma_alpha in Eqn. (20) in the implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "alpha_min_frac", &tmp_double, &tmp_res, 0.05, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("alpha_min_frac", tmp_double);
  // * none: no corrector
  // * affine: corrector step towards mu=0
  // * primal-dual: corrector step towards current mu
  getStringInPList(pvApiCtx, param_in_addr, "corrector_type", &tmp_char, &tmp_res, "none", Log, 
		       CHECK_VALUES, 3,
		       "none",
		       "affine",
		       "primal-dual");
  if (tmp_res!=-1) options->SetStringValue("corrector_type", tmp_char);
  FREE(tmp_char);
  // Skip the corrector step in negative curvature iteration (unsupported!).
  //   The corrector step is not tried if negative curvature has been
  //   encountered during the computation of the search direction in the current
  //   iteration. This option is only used if "mu_strategy" is "adaptive".
  // Possible values:
  //  - no                      [don't skip]
  //  - yes                     [skip]
  getStringInPList(pvApiCtx, param_in_addr, "skip_corr_if_neg_curv", &tmp_char, &tmp_res, "yes", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("skip_corr_if_neg_curv", tmp_char);
  FREE(tmp_char);
  // Skip the corrector step during monotone barrier parameter mode (unsupported!).
  //   The corrector step is not tried if the algorithm is currently in the
  //   monotone mode (see also option "barrier_strategy").This option is only
  //   used if "mu_strategy" is "adaptive".
  // Possible values:
  //  - no                      [don't skip]
  //  - yes                     [skip]
  getStringInPList(pvApiCtx, param_in_addr, "skip_corr_in_monotone_mode", &tmp_char, &tmp_res, "yes", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("skip_corr_in_monotone_mode", tmp_char);
  FREE(tmp_char);
  // Complementarity tolerance factor for accepting corrector step (unsupported!).
  //   This option determines the factor by which complementarity is allowed to
  //   increase for a corrector step to be accepted.
  getDoubleInPList(pvApiCtx, param_in_addr, "corrector_compl_avrg_red_fact", &tmp_double, &tmp_res, 1, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("corrector_compl_avrg_red_fact", tmp_double);
  // Initial value of the penalty parameter.
  getDoubleInPList(pvApiCtx, param_in_addr, "nu_init", &tmp_double, &tmp_res, 1e-6, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("nu_init", tmp_double);
  // Increment of the penalty parameter.
  getDoubleInPList(pvApiCtx, param_in_addr, "nu_inc", &tmp_double, &tmp_res, 0.0001, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("nu_inc", tmp_double);
  // Value in penalty parameter update formula.
  getDoubleInPList(pvApiCtx, param_in_addr, "rho", &tmp_double, &tmp_res, 0.1, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("rho", tmp_double);

  ////////////////
  // Warm Start //
  ////////////////
  // Enables to specify bound multiplier values
  // * no: do not use the warm start initialization
  // * yes: use the warm start initialization
  getStringInPList(pvApiCtx, param_in_addr, "warm_start_init_point", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("warm_start_init_point", tmp_char);
  FREE(tmp_char);
  // Indicates whether a problem with a structure identical to the previous one is to be solved.
  //   If "yes" is chosen, then the algorithm assumes that an NLP is now to be
  //   solved, whose structure is identical to one that already was considered
  //   (with the same NLP object).
  // Possible values:
  //  - no                      [Assume this is a new problem.]
  //  - yes                     [Assume this is problem has known structure]
  getStringInPList(pvApiCtx, param_in_addr, "warm_start_same_structure", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("warm_start_same_structure", tmp_char);
  FREE(tmp_char);
  // Enables to specify how much should variables should be pushed inside the feasible region
  getDoubleInPList(pvApiCtx, param_in_addr, "warm_start_bound_push", &tmp_double, &tmp_res, 1e-3, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("warm_start_bound_push", tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "warm_start_bound_frac", &tmp_double, &tmp_res, 1e-3, Log, CHECK_BOTH, 0, 0.5);
  if (tmp_res!=-1) options->SetNumericValue("warm_start_bound_frac", tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "warm_start_slack_bound_frac", &tmp_double, &tmp_res, 1e-3, Log, CHECK_BOTH, 0, 0.5);
  if (tmp_res!=-1) options->SetNumericValue("warm_start_slack_bound_frac", tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "warm_start_slack_bound_push", &tmp_double, &tmp_res, 1e-3, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("warm_start_slack_bound_push", tmp_double);
  // Enables to specify how much should bound multipliers should be pushed inside the feasible region
  getDoubleInPList(pvApiCtx, param_in_addr, "warm_start_mult_bound_push", &tmp_double, &tmp_res, 1e-3, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("warm_start_mult_bound_push", tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "warm_start_mult_init_max", &tmp_double, &tmp_res, 1e6, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("warm_start_mult_init_max", tmp_double);
  // Tells algorithm whether to use the GetWarmStartIterate method in the NLP.
  // Possible values:
  //  - no                      [call GetStartingPoint in the NLP]
  //  - yes                     [call GetWarmStartIterate in the NLP]
  getStringInPList(pvApiCtx, param_in_addr, "warm_start_entire_iterate", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("warm_start_entire_iterate", tmp_char);
  FREE(tmp_char);

  getDoubleInPList(pvApiCtx, param_in_addr, "warm_start_target_mu", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("warm_start_target_mu", tmp_double);

  ///////////////////////
  // Restoration Phase //
  ///////////////////////
  // Enable heuristics to quickly detect an infeasible problem
  // * no: the problem probably be feasible
  // * yes: the problem has a good chance to be infeasible
  // This options is meant to activate heuristics that may speed up the infeasibility determination if you expect that there is a good chance for the problem to be
  // infeasible. In the filter line search procedure, the restoration phase is called more quickly than usually, and more reduction in 
  // the constraint violation is enforced before the restoration phase is left. If the problem is square, this option is enabled automatically.
  getStringInPList(pvApiCtx, param_in_addr, "expect_infeasible_problem", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("expect_infeasible_problem", tmp_char);
  FREE(tmp_char);
  // Threshold for disabling "expect_infeasible_problem" option.
  // If the constraint violation becomes smaller than this threshold, the "expect_infeasible_problem" heuristics in the filter line 
  // search are disabled. If the problem is square, this options is set to 0.
  getDoubleInPList(pvApiCtx, param_in_addr, "expect_infeasible_problem_ctol", &tmp_double, &tmp_res, 1e3, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("expect_infeasible_problem_ctol", tmp_double);
  // Multiplier threshold for activating "expect_infeasible_problem" option.
  //   If the max norm of the constraint multipliers becomes larger than this
  //   value and "expect_infeasible_problem" is chosen, then the restoration phase is entered.
  getDoubleInPList(pvApiCtx, param_in_addr, "expect_infeasible_problem_ytol", &tmp_double, &tmp_res, 1e8, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("expect_infeasible_problem_ytol", tmp_double);
  // Tells algorithm to switch to restoration phase in first iteration.
  // * no: don't force start in restoration phase
  // * yes: force start in restoration phase
  // Setting this option to "yes" forces the algorithm to switch to the feasibility restoration phase in the first iteration. If the initial 
  // point is feasible, the algorithm will abort with a failure.
  getStringInPList(pvApiCtx, param_in_addr, "start_with_resto", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("start_with_resto", tmp_char);
  FREE(tmp_char);
  // Required reduction in primal-dual error in the soft restoration phase.
  // The soft restoration phase attempts to reduce the primal-dual error with regular steps. If the damped 
  // primal-dual step (damped only to satisfy the fraction-to-the-boundary rule) is not decreasing the primal-dual error 
  // by at least this factor, then the regular restoration phase is called. Choosing "0" here disables the soft restoration phase.
  getDoubleInPList(pvApiCtx, param_in_addr, "soft_resto_pderror_reduction_factor", &tmp_double, &tmp_res, 0.9999, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("soft_resto_pderror_reduction_factor", tmp_double);
  // Maximum number of iterations performed successively in soft restoration phase.
  // If the soft restoration phase is performed for more than so many iterations in a row, the regular restoration phase is called.
  getIntInPList(pvApiCtx, param_in_addr, "max_soft_resto_iters", &tmp_int, &tmp_res, 10, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetIntegerValue("max_soft_resto_iters", tmp_int);
  // Required infeasibility reduction in restoration phase
  // The restoration phase algorithm is performed, until a point is found that is acceptable to the filter and the infeasibility has been
  // reduced by at least the fraction given by this option.
  getDoubleInPList(pvApiCtx, param_in_addr, "required_infeasibility_reduction", &tmp_double, &tmp_res, 0.9, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("required_infeasibility_reduction", tmp_double);
  // Maximum number of successive iterations in restoration phase.
  // The algorithm terminates with an error message if the number of iterations successively taken in the restoration phase exceeds this number.
  getIntInPList(pvApiCtx, param_in_addr, "max_resto_iter", &tmp_int, &tmp_res, 3000000, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetIntegerValue("max_resto_iter", tmp_int);
  // Threshold for resetting bound multipliers after the restoration phase.
  // After returning from the restoration phase, the bound multipliers are updated with a Newton step for complementarity. Here, the
  // change in the primal variables during the entire restoration phase is taken to be the corresponding primal Newton step. 
  // However, if after the update the largest bound multiplier exceeds the threshold specified by this option, the multipliers are all reset to 1.
  getDoubleInPList(pvApiCtx, param_in_addr, "bound_mult_reset_threshold", &tmp_double, &tmp_res, 1000, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("bound_mult_reset_threshold", tmp_double);
  // Threshold for resetting equality and inequality multipliers after restoration phase.
  // After returning from the restoration phase, the constraint multipliers are recomputed by a least square estimate.  This option triggers when
  // those least-square estimates should be ignored.
  getDoubleInPList(pvApiCtx, param_in_addr, "constr_mult_reset_threshold", &tmp_double, &tmp_res, 0, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("constr_mult_reset_threshold", tmp_double);
  // Determines if the original objective function should be evaluated at restoration phase trial points.
  // Setting this option to "yes" makes the restoration phase algorithm evaluate the objective function of the original problem at every trial 
  // point encountered during the restoration phase, even if this value is not required. In this way, it is guaranteed that the original 
  // objective function can be evaluated without error at all accepted iterates; otherwise the algorithm might fail at a point where the 
  // restoration phase accepts an iterate that is good for the restoration phase problem, but not the original problem. On the other hand, if 
  // the evaluation of the original objective is expensive, this might be costly.
  // * no: skip evaluation
  // * yes: evaluate at every trial point
  getStringInPList(pvApiCtx, param_in_addr, "evaluate_orig_obj_at_resto_trial", &tmp_char, &tmp_res, "yes", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("evaluate_orig_obj_at_resto_trial", tmp_char);
  FREE(tmp_char);
  // Penalty parameter in the restoration phase objective function.
  // This is the parameter rho in equation (31a) in the Ipopt implementation paper.
  getDoubleInPList(pvApiCtx, param_in_addr, "resto_penalty_parameter", &tmp_double, &tmp_res, 1000.0, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("resto_penalty_parameter", tmp_double);

  ///////////////////
  // Linear Solver //
  ///////////////////
  // Linear solver to be used for step calculation
  // * ma27: use the Harwell routine MA27
  // * ma57: use the Harwell routine MA57
  // * pardiso: use the Pardiso package
  // * wsmp: use WSMP package
  // * mumps: use MUMPS package
  // * custom: use custom linear solver
  // Determines which linear algebra package is to be used for the solution of the augmented linear system (for obtaining the search directions).
  // Note, the code must have been compiled with the linear solver you want to choose. Depending on your Ipopt installation, not all options are available.
  getStringInPList(pvApiCtx, param_in_addr, "linear_solver", &tmp_char, &tmp_res, "mumps", Log, 
		       CHECK_VALUES, 6,
		       "ma27",
		       "ma57",
		       "pardiso",
		       "wsmp",
		       "mumps",
		       "custom");
  if (tmp_res!=-1) options->SetStringValue("linear_solver", tmp_char);
  FREE(tmp_char);
  // Method for scaling the linear systems
  // * none: no scaling will be performed
  // * mc19: use the Harwell routine MC19
  // "Determines the method used to compute symmetric scaling factors for the augmented system (see also the "linear_scaling_on_demand" option).
  // This scaling is independent of the NLP problem scaling.  By default, MC19 is only used if MA27 or MA57 are selected as linear solvers. 
  // This option is only available if Ipopt has been compiled with MC19.
  getStringInPList(pvApiCtx, param_in_addr, "linear_system_scaling", &tmp_char, &tmp_res, "none", Log, 
		       CHECK_VALUES, 2, 
		       "none",
		       "mc19");
  if (tmp_res!=-1) options->SetStringValue("linear_system_scaling", tmp_char);
  FREE(tmp_char);
  // Enables heuristic for scaling only when seems required
  // * no: Always scale the linear system.
  // * yes: Start using linear system scaling if solutions seem not good.
  getStringInPList(pvApiCtx, param_in_addr, "linear_scaling_on_demand", &tmp_char, &tmp_res, "yes", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("linear_scaling_on_demand", tmp_char);
  FREE(tmp_char);
  // Maximal number of iterative refinement steps per linear system solve
  getIntInPList(pvApiCtx, param_in_addr, "max_refinement_steps", &tmp_int, &tmp_res, 10, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetIntegerValue("max_refinement_steps", tmp_int);
  // Minimum number of iterative refinement steps per linear system solve
  getIntInPList(pvApiCtx, param_in_addr, "min_refinement_steps", &tmp_int, &tmp_res, 1, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetIntegerValue("min_refinement_steps", tmp_int);
  // Iterative refinement tolerance
  //   Iterative refinement is performed until the residual test ratio is less
  //   than this tolerance (or until "max_refinement_steps" refinement steps are performed).
  getDoubleInPList(pvApiCtx, param_in_addr, "residual_ratio_max", &tmp_double, &tmp_res, 1e-10, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("residual_ratio_max", tmp_double);
  // Threshold for declaring linear system singular after failed iterative refinement.
  //   If the residual test ratio is larger than this value after failed
  //   iterative refinement, the algorithm pretends that the linear system is singular.
  getDoubleInPList(pvApiCtx, param_in_addr, "residual_ratio_singular", &tmp_double, &tmp_res, 1e-5, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("residual_ratio_singular", tmp_double);
  // Minimal required reduction of residual test ratio in iterative refinement.
  //   If the improvement of the residual test ratio made by one iterative
  //   refinement step is not better than this factor, iterative refinement is aborted.
  getDoubleInPList(pvApiCtx, param_in_addr, "residual_improvement_factor", &tmp_double, &tmp_res, 1, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("residual_improvement_factor", tmp_double);
  // Tolerance for heuristic to ignore wrong inertia.
  //   If positive, incorrect inertia in the augmented system is ignored, and we
  //   test if the direction is a direction of positive curvature.  This
  //   tolerance determines when the direction is considered to be sufficiently positive.
  getDoubleInPList(pvApiCtx, param_in_addr, "neg_curv_test_tol", &tmp_double, &tmp_res, 0, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("neg_curv_test_tol", tmp_double);

  //////////////////////////
  // Hessian Perturbation //
  //////////////////////////
  // Maximum value of regularization parameter for handling negative curvature.
  // In order to guarantee that the search directions are indeed proper descent directions, Ipopt requires that the inertia of the 
  // (augmented) linear system for the step computation has the correct number of negative and positive eigenvalues. The idea 
  // is that this guides the algorithm away from maximizers and makes Ipopt more likely converge to first order optimal points that 
  // are minimizers. If the inertia is not correct, a multiple of the identity matrix is added to the Hessian of the Lagrangian in the 
  // augmented system. This parameter gives the maximum value of the regularization parameter. If a regularization of that size is 
  // not enough, the algorithm skips this iteration and goes to the restoration phase. (This is delta_w^max in the implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "max_hessian_perturbation", &tmp_double, &tmp_res, 1e20, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("max_hessian_perturbation", tmp_double);
  // Smallest perturbation of the Hessian block.
  // The size of the perturbation of the Hessian block is never selected smaller than this value, unless no perturbation is necessary. (This 
  // is delta_w^min in implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "min_hessian_perturbation", &tmp_double, &tmp_res, 1e-20, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("min_hessian_perturbation", tmp_double);
  // Size of first x-s perturbation tried.
  // The first value tried for the x-s perturbation in the inertia correction scheme.
  // (This is delta_0 in the implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "first_hessian_perturbation", &tmp_double, &tmp_res, 1e-4, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("first_hessian_perturbation", tmp_double);
  // Increase factor for x-s perturbation for very first perturbation.
  // The factor by which the perturbation is increased when a trial value was not sufficient - this value is used for the computation of the 
  // very first perturbation and allows a different value for for the first perturbation than that used for the remaining perturbations. 
  // (This is bar_kappa_w^+ in the implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "perturb_inc_fact_first", &tmp_double, &tmp_res, 100, Log, CHECK_MIN, 1);
  if (tmp_res!=-1) options->SetNumericValue("perturb_inc_fact_first", tmp_double);
  // Increase factor for x-s perturbation.
  // The factor by which the perturbation is increased when a trial value was not sufficient - this value is used for the computation of 
  // all perturbations except for the first. (This is kappa_w^+ in the implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "perturb_inc_fact", &tmp_double, &tmp_res, 8, Log, CHECK_MIN, 1);
  if (tmp_res!=-1) options->SetNumericValue("perturb_inc_fact", tmp_double);
  // Decrease factor for x-s perturbation.
  // The factor by which the perturbation is decreased when a trial value is deduced from the size of the most recent successful perturbation. 
  // (This is kappa_w^- in the implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "perturb_dec_fact", &tmp_double, &tmp_res, 0.333333, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("perturb_dec_fact", tmp_double);
  // Size of the regularization for rank-deficient constraint Jacobians.
  // (This is bar delta_c in the implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "jacobian_regularization_value", &tmp_double, &tmp_res, 1e-8, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("jacobian_regularization_value", tmp_double);
  // Exponent for mu in the regularization for rank-deficient constraint Jacobians.
  // (This is kappa_c in the implementation paper.)
  getDoubleInPList(pvApiCtx, param_in_addr, "jacobian_regularization_exponent", &tmp_double, &tmp_res, 0.25, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("jacobian_regularization_exponent", tmp_double);
  // Active permanent perturbation of constraint linearization.
  // no: perturbation only used when required
  // yes: always use perturbation
  // This options makes the delta_c and delta_d perturbation be used for the computation of every search direction.  Usually, it is only used 
  // when the iteration matrix is singular.
  getStringInPList(pvApiCtx, param_in_addr, "perturb_always_cd", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("perturb_always_cd", tmp_char);
  FREE(tmp_char);

  //////////////////
  // Quasi-Newton //
  //////////////////
  // Can enable Quasi-Newton approximation of hessian
  // * exact: Use second derivatives provided by the NLP.
  // * limited-memory: Perform a limited-memory quasi-Newton approximation
  // This determines which kind of information for the Hessian of the Lagrangian function is used by the algorithm.
  getStringInPList(pvApiCtx, param_in_addr, "hessian_approximation", &tmp_char, &tmp_res, "limited-memory", Log, 
		       CHECK_VALUES, 2,
		       "exact",
		       "limited-memory");
  if (tmp_res!=-1) options->SetStringValue("hessian_approximation", tmp_char);
  FREE(tmp_char);
  // Indicates in which subspace the Hessian information is to be approximated.
  // nonlinear-variables: only in space of nonlinear variables.
  // all-variables: in space of all variables (without slacks)
  getStringInPList(pvApiCtx, param_in_addr, "hessian_approximation_space", &tmp_char, &tmp_res, "nonlinear-variables", Log, 
		       CHECK_VALUES, 2,
		       "nonlinear-variables",
		       "all-variables");
  if (tmp_res!=-1) options->SetStringValue("hessian_approximation_space", tmp_char);
  FREE(tmp_char);
  // Maximum size of the history for the limited quasi-Newton Hessian approximation.
  // This option determines the number of most recent iterations that are taken into account for the limited-memory quasi-Newton approximation.
  getIntInPList(pvApiCtx, param_in_addr, "limited_memory_max_history", &tmp_int, &tmp_res, 6, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetIntegerValue("limited_memory_max_history", tmp_int);
  // Threshold for successive iterations where update is skipped.
  // If the update is skipped more than this number of successive iterations, we quasi-Newton approximation is reset.
  getIntInPList(pvApiCtx, param_in_addr, "limited_memory_max_skipping", &tmp_int, &tmp_res, 2, Log, CHECK_MIN, 1);
  if (tmp_res!=-1) options->SetIntegerValue("limited_memory_max_skipping", tmp_int);
  // Quasi-Newton update formula for the limited memory approximation.
  // bfgs: BFGS update (with skipping)
  // sr1:  SR1 (not working well)
  // Determines which update formula is to be used for the limited-memory quasi-Newton approximation.
  getStringInPList(pvApiCtx, param_in_addr, "limited_memory_update_type", &tmp_char, &tmp_res, "bfgs", Log, 
		       CHECK_VALUES, 2,
		       "bfgs",
		       "sr1");
  if (tmp_res!=-1) options->SetStringValue("limited_memory_update_type", tmp_char);
  FREE(tmp_char);
  // Initialization strategy for the limited memory quasi-Newton approximation.
  // scalar1: sigma = s^Ty/s^Ts
  // scalar2: sigma = y^Ty/s^Ty
  // constant: sigma = limited_memory_init_val
  // Determines how the diagonal Matrix B_0 as the first term in the limited memory approximation should be computed.
  getStringInPList(pvApiCtx, param_in_addr, "limited_memory_initialization", &tmp_char, &tmp_res, "scalar1", Log, 
		       CHECK_VALUES, 3, 
		       "scalar1",
		       "scalar2",
		       "constant");
  if (tmp_res!=-1) options->SetStringValue("limited_memory_initialization", tmp_char);
  FREE(tmp_char);
  // Upper bound on value for B0 in low-rank update.
  // The starting matrix in the low rank update, B0, is chosen to be this multiple of the identity in the first iteration (when no updates have 
  // been performed yet), and is constantly chosen as this value, if "limited_memory_initialization" is "constant".
  getDoubleInPList(pvApiCtx, param_in_addr, "limited_memory_init_val_max", &tmp_double, &tmp_res, 1e8, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("limited_memory_init_val_max", tmp_double);
  // Lower bound on value for B0 in low-rank update.
  // The starting matrix in the low rank update, B0, is chosen to be this multiple of the identity in the first iteration (when no updates have 
  // been performed yet), and is constantly chosen as this value, if "limited_memory_initialization" is "constant".
  getDoubleInPList(pvApiCtx, param_in_addr, "limited_memory_init_val_min", &tmp_double, &tmp_res, 1e-8, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("limited_memory_init_val_min", tmp_double);
  // Value for B0 in low-rank update
  // The starting matrix in the low rank update, B0, is chosen to be this multiple of the identity in the first iteration (when no updates have
  // "been performed yet), and is constantly chosen as this value, if "limited_memory_initialization" is "constant".
  getDoubleInPList(pvApiCtx, param_in_addr, "limited_memory_init_val", &tmp_double, &tmp_res, 1, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("limited_memory_init_val", tmp_double);

  /////////////////////
  // Derivative Test //
  /////////////////////
  // Enable derivative checker
  // * none: do not perform derivative test
  // * first-order: perform test of first derivatives at starting point
  // * second-order: perform test of first and second derivatives at starting point
  // If this option is enabled, a (slow) derivative test will be performed before the optimization. The test is performed at the user provided 
  // starting point and marks derivative values that seem suspicious
  getStringInPList(pvApiCtx, param_in_addr, "derivative_test", &tmp_char, &tmp_res, "none", Log, 
		       CHECK_VALUES, 4,
		       "none",
		       "first-order",
		       "second-order",
		       "only-second-order");
  if (tmp_res!=-1) options->SetStringValue("derivative_test", tmp_char);
  FREE(tmp_char);
  // Index of first quantity to be checked by derivative checker
  //   If this is set to -2, then all derivatives are checked.  Otherwise, for
  //   the first derivative test it specifies the first variable for which the
  //   test is done (counting starts at 0).  For second derivatives, it
  //   specifies the first constraint for which the test is done; counting of
  //   constraint indices starts at 0, and -1 refers to the objective function
  //   Hessian.
  getIntInPList(pvApiCtx, param_in_addr, "derivative_test_first_index", &tmp_int, &tmp_res, -2, Log, CHECK_MIN, -2);
  if (tmp_res!=-1) options->SetIntegerValue("derivative_test_first_index", tmp_int);
  // Size of the finite difference perturbation in derivative test.
  // This determines the relative perturbation of the variable entries.
  getDoubleInPList(pvApiCtx, param_in_addr, "derivative_test_perturbation", &tmp_double, &tmp_res, 1e-8, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("derivative_test_perturbation", tmp_double);
  // Threshold for indicating wrong derivative.
  // If the relative deviation of the estimated derivative from the given one is larger than this value, the corresponding derivative is marked as wrong.
  getDoubleInPList(pvApiCtx, param_in_addr, "derivative_test_tol", &tmp_double, &tmp_res, 1e-4, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("derivative_test_tol", tmp_double);
  // Indicates whether information for all estimated derivatives should be printed.
  // * no: Print only suspect derivatives
  // * yes: Print all derivatives
  getStringInPList(pvApiCtx, param_in_addr, "derivative_test_print_all", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("derivative_test_print_all", tmp_char);
  FREE(tmp_char);
  // Maximal perturbation of an evaluation point.
  // If a random perturbation of a points is required, this number indicates the maximal perturbation. This is for example used when 
  // determining the center point at which the finite difference derivative test is executed.
  getDoubleInPList(pvApiCtx, param_in_addr, "point_perturbation_radius", &tmp_double, &tmp_res, 10, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("point_perturbation_radius", tmp_double);
  // Specifies technique to compute constraint Jacobian
  // exact: user-provided derivatives
  // finite-difference-values: user-provided structure, values by finite differences
  getStringInPList(pvApiCtx, param_in_addr, "jacobian_approximation", &tmp_char, &tmp_res, "exact", Log, 
		       CHECK_VALUES, 2, 
		       "exact",
		       "finite-difference-values");
  if (tmp_res!=-1) options->SetStringValue("jacobian_approximation", tmp_char);
  FREE(tmp_char);
  // Size of the finite difference perturbation for derivative approximation.
  // This determines the relative perturbation of the variable entries.
  getDoubleInPList(pvApiCtx, param_in_addr, "findiff_perturbation", &tmp_double, &tmp_res, 1e-7, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("findiff_perturbation", tmp_double);

#ifdef USE_MA27
  ////////////////////////
  // MA27 Linear Solver //
  ////////////////////////
  // Pivot tolerance for the linear solver MA27
  // A smaller number pivots for sparsity, a larger number pivots for stability.  This option is only available if Ipopt has been compiled with MA27.
  getDoubleInPList(pvApiCtx, param_in_addr, "ma27_pivtol", &tmp_double, &tmp_res, 1e-8, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("ma27_pivtol", tmp_double);
  // Maximal pivot tolerance for the linear solver MA27
  // Ipopt may increase pivtol as high as pivtolmax to get a more accurate solution to the linear system.  This option is only available if 
  // Ipopt has been compiled with MA27.
  getDoubleInPList(pvApiCtx, param_in_addr, "ma27_pivtolmax", &tmp_double, &tmp_res, 1e-4, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("ma27_pivtolmax", tmp_double);
  // Integer workspace memory for MA27.
  // The initial integer workspace memory = liw_init_factor * memory required by unfactored system. Ipopt will increase the workspace 
  // size by meminc_factor if required.  This option is only available if Ipopt has been compiled with MA27.
  getDoubleInPList(pvApiCtx, param_in_addr, "ma27_liw_init_factor", &tmp_double, &tmp_res, 5, Log, CHECK_MIN, 1);
  if (tmp_res!=-1) options->SetNumericValue("ma27_liw_init_factor", tmp_double);
  // Real workspace memory for MA27.
  // The initial real workspace memory = la_init_factor * memory required by unfactored system. Ipopt will increase the workspace
  // size by meminc_factor if required.  This option is only available if Ipopt has been compiled with MA27.
  getDoubleInPList(pvApiCtx, param_in_addr, "ma27_la_init_factor", &tmp_double, &tmp_res, 5, Log, CHECK_MIN, 1);
  if (tmp_res!=-1) options->SetNumericValue("ma27_la_init_factor", tmp_double);
  // Increment factor for workspace size for MA27.
  // If the integer or real workspace is not large enough, Ipopt will increase its size by this factor.  This option is only 
  // available if Ipopt has been compiled with MA27.
  getDoubleInPList(pvApiCtx, param_in_addr, "ma27_meminc_factor", &tmp_double, &tmp_res, 5, Log, CHECK_MIN, 1);
  if (tmp_res!=-1) options->SetNumericValue("ma27_meminc_factor", tmp_double);
  // Always pretend inertia is correct.
  // no: check inertia
  // yes: skip inertia check
  // Setting this option to "yes" essentially disables inertia check. This option makes the algorithm non-robust and easily fail, but it 
  // might give some insight into the necessity of inertia control.
  getStringInPList(pvApiCtx, param_in_addr, "ma27_skip_inertia_check", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("ma27_skip_inertia_check", tmp_char);
  FREE(tmp_char);
  // Enables MA27's ability to solve a linear system even if the matrix is singular.
  // no: Don't have MA27 solve singular systems
  // yes: Have MA27 solve singular systems
  // Setting this option to \"yes\" means that Ipopt will call MA27 to compute solutions for right hand sides, even if MA27 has detected that 
  // the matrix is singular (but is still able to solve the linear system). In some cases this might be better than using Ipopt's heuristic of 
  // small perturbation of the lower diagonal of the KKT matrix.
  getStringInPList(pvApiCtx, param_in_addr, "ma27_ignore_singularity", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("ma27_ignore_singularity", tmp_char);
  FREE(tmp_char);
#endif

#ifdef USE_MA28
  ////////////////////////
  // MA28 Linear Solver //
  ////////////////////////
  // Pivot tolerance for linear solver MA28.
  //   This is used when MA28 tries to find the dependent constraints.
  getDoubleInPList(pvApiCtx, param_in_addr, "ma28_pivtol", &tmp_double, &tmp_res, 0.01, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("ma27_pivtol", tmp_double);
#endif

#ifdef USE_MA57
  ////////////////////////
  // MA57 Linear Solver //
  ////////////////////////
  // Pivot tolerance for the linear solver MA57
  // A smaller number pivots for sparsity, a larger number pivots for stability. This option is only available if Ipopt has been compiled with MA57.
  getDoubleInPList(pvApiCtx, param_in_addr, "ma57_pivtol", &tmp_double, &tmp_res, 1e-8, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("ma57_pivtol", tmp_double);
  // Maximal pivot tolerance for the linear solver MA57
  // Ipopt may increase pivtol as high as ma57_pivtolmax to get a more accurate solution to the linear system.  This option is only available 
  // if Ipopt has been compiled with MA57.
  getDoubleInPList(pvApiCtx, param_in_addr, "ma57_pivtolmax", &tmp_double, &tmp_res, 1e-4, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("ma57_pivtolmax", tmp_double);
  // Safety factor for work space memory allocation for the linear solver MA57.
  // If 1 is chosen, the suggested amount of work space is used.  However, choosing a larger number might avoid reallocation if the suggest values 
  // do not suffice.  This option is only available if Ipopt has been compiled with MA57.
  getIntInPList(pvApiCtx, param_in_addr, "ma57_pre_alloc", &tmp_int, &tmp_res, 3, Log, CHECK_MIN, 1);
  if (tmp_res!=-1) options->SetIntegerValue("ma57_pre_alloc", tmp_int);
  // Controls pivot order in MA57
  //   This is INCTL(6) in MA57.
  getIntInPList(pvApiCtx, param_in_addr, "ma57_pivot_order", &tmp_int, &tmp_res, 5, Log, CHECK_BOTH, 0, 5);
  if (tmp_res!=-1) options->SetIntegerValue("ma57_pivot_order", tmp_int);
#endif

#ifdef USE_MUMPS
  /////////////////////////
  // MUMPS Linear Solver //
  /////////////////////////
  // Pivot tolerance for the linear solver MUMPS.
  // A smaller number pivots for sparsity, a larger number pivots for stability.  This option is only available if Ipopt has been compiled with MUMPS.
  getDoubleInPList(pvApiCtx, param_in_addr, "mumps_pivtol", &tmp_double, &tmp_res, 1e-6, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("mumps_pivtol", tmp_double);
  // Maximum pivot tolerance for the linear solver MUMPS.
  // Ipopt may increase pivtol as high as pivtolmax to get a more accurate solution to the linear system.  This option is only available if 
  // Ipopt has been compiled with MUMPS.
  getDoubleInPList(pvApiCtx, param_in_addr, "mumps_pivtolmax", &tmp_double, &tmp_res, 0.1, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("mumps_pivtolmax", tmp_double);
  // Percentage increase in the estimated working space for MUMPS.
  // In MUMPS when significant extra fill-in is caused by numerical pivoting, larger values of mumps_mem_percent may help use the workspace more efficiently.
  getIntInPList(pvApiCtx, param_in_addr, "mumps_mem_percent", &tmp_int, &tmp_res, 1000, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetIntegerValue("mumps_mem_percent", tmp_int);
  // Controls permuting and scaling in MUMPS
  // This is ICTL(6) in MUMPS.
  getIntInPList(pvApiCtx, param_in_addr, "mumps_permuting_scaling", &tmp_int, &tmp_res, 7, Log, CHECK_BOTH, 0, 7);
  if (tmp_res!=-1) options->SetIntegerValue("mumps_permuting_scaling", tmp_int);
  // Controls pivot order in MUMPS
  // This is ICTL(8) in MUMPS.
  getIntInPList(pvApiCtx, param_in_addr, "mumps_pivot_order", &tmp_int, &tmp_res, 7, Log, CHECK_BOTH, 0, 7);
  if (tmp_res!=-1) options->SetIntegerValue("mumps_pivot_order", tmp_int);
  // Controls scaling in MUMPS
  // This is ICTL(8) in MUMPS.
  getIntInPList(pvApiCtx, param_in_addr, "mumps_scaling", &tmp_int, &tmp_res, 77, Log, CHECK_BOTH, -2, 77);
  if (tmp_res!=-1) options->SetIntegerValue("mumps_scaling", tmp_int);
  // Pivot threshold for detection of linearly dependent constraints in MUMPS.
  // When MUMPS is used to determine linearly dependent constraints, this is determines the threshold for a pivot to be considered zero. This is CNTL(3) in MUMPS.
  getDoubleInPList(pvApiCtx, param_in_addr, "mumps_dep_tol", &tmp_double, &tmp_res, -1.0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("mumps_dep_tol", tmp_double);
#endif

#ifdef USE_PARDISO
  ///////////////////////////
  // Pardiso Linear Solver //
  ///////////////////////////
  // Matching strategy for linear solver Pardiso
  // * complete: Match complete (IPAR(13)=1)
  // * complete+2x2: Match complete+2x2 (IPAR(13)=2)
  // * constraints: Match constraints (IPAR(13)=3)
  // This is IPAR(13) in Pardiso manual.  This option is only available if Ipopt has been compiled with Pardiso.
  getStringInPList(pvApiCtx, param_in_addr, "pardiso_matching_strategy", &tmp_char, &tmp_res, "complete+2x2", Log, 
		       CHECK_VALUES, 3,
		       "complete",
		       "complete+2x2",
		       "constraints");
  if (tmp_res!=-1) options->SetStringValue("pardiso_matching_strategy", tmp_char);
  FREE(tmp_char);
  // Enables out-of-core version of linear solver Pardiso
  // Setting this option to a positive integer k makes Pardiso work in the out-of-core variant where the factor is split in 2^k subdomains. This 
  // is IPARM(50) in the Pardiso manual.  This option is only available if Ipopt has been compiled with Pardiso.
  getIntInPList(pvApiCtx, param_in_addr, "pardiso_out_of_core_power", &tmp_int, &tmp_res, 0, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetIntegerValue("pardiso_out_of_core_power", tmp_int);
  // Pardiso message level
  //   This determines the amount of analysis output from the Pardiso solver.
  //   This is MSGLVL in the Pardiso manual.
  getIntInPList(pvApiCtx, param_in_addr, "pardiso_msglvl", &tmp_int, &tmp_res, 0, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetIntegerValue("pardiso_msglvl", tmp_int);
  // Maximum number of Krylov-Subspace Iteration
  //   DPARM(1)
  getIntInPList(pvApiCtx, param_in_addr, "pardiso_max_iter", &tmp_int, &tmp_res, 500, Log, CHECK_MIN, 1);
  if (tmp_res!=-1) options->SetIntegerValue("pardiso_max_iter", tmp_int);
  // Relative Residual Convergence
  //   DPARM(2)
  getDoubleInPList(pvApiCtx, param_in_addr, "pardiso_iter_relative_tol", &tmp_double, &tmp_res, 1e-6, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("pardiso_iter_relative_tol", tmp_double);
  // Maximum Size of Coarse Grid Matrix
  //   DPARM(3)
  getDoubleInPList(pvApiCtx, param_in_addr, "pardiso_iter_coarse_size", &tmp_double, &tmp_res, 5000, Log, CHECK_MIN, 1);
  if (tmp_res!=-1) options->SetNumericValue("pardiso_iter_coarse_size", tmp_double);
  // Maximum Size of Grid Levels
  //   DPARM(4)
  getDoubleInPList(pvApiCtx, param_in_addr, "pardiso_iter_max_levels", &tmp_double, &tmp_res, 10000, Log, CHECK_MIN, 1);
  if (tmp_res!=-1) options->SetNumericValue("pardiso_iter_max_levels", tmp_double);
  // dropping value for incomplete factor
  //   DPARM(5)
  getDoubleInPList(pvApiCtx, param_in_addr, "pardiso_iter_dropping_factor", &tmp_double, &tmp_res, 0.5, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("pardiso_iter_dropping_factor", tmp_double);
  // dropping value for sparsify schur complement factor
  //   DPARM(6)
  getDoubleInPList(pvApiCtx, param_in_addr, "pardiso_iter_dropping_schur", &tmp_double, &tmp_res, 0.1, Log, CHECK_BOTH, 0, 1);
  if (tmp_res!=-1) options->SetNumericValue("pardiso_iter_dropping_schur", tmp_double);
  // max fill for each row
  //   DPARM(7)
  getDoubleInPList(pvApiCtx, param_in_addr, "pardiso_iter_max_row_fill", &tmp_double, &tmp_res, 10000000, Log, CHECK_MIN, 1);
  if (tmp_res!=-1) options->SetNumericValue("pardiso_iter_max_row_fill", tmp_double);
  //   DPARM(8)
  getDoubleInPList(pvApiCtx, param_in_addr, "pardiso_iter_inverse_norm_factor", &tmp_double, &tmp_res, 5e6, Log, CHECK_MIN, 1);
  if (tmp_res!=-1) options->SetNumericValue("pardiso_iter_inverse_norm_factor", tmp_double);
  // Switch on iterative solver in Pardiso library
  // Possible values:
  //  - no                      []
  //  - yes                     []
  getStringInPList(pvApiCtx, param_in_addr, "pardiso_iterative", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("pardiso_iterative", tmp_char);
  FREE(tmp_char);
  // Maximal number of decreases of drop tolerance during one solve.
  //   This is relevant only for iterative Pardiso options.
  getIntInPList(pvApiCtx, param_in_addr, "pardiso_max_droptol_corrections", &tmp_int, &tmp_res, 4, Log, CHECK_MIN, 1);
  if (tmp_res!=-1) options->SetIntegerValue("pardiso_max_droptol_corrections", tmp_int);

  // Toggle for handling case when elements were perturbed by Pardiso.
  // no: Always redo symbolic factorization when elements were perturbed
  // yes: Only redo symbolic factorization when elements were perturbed if also the inertia was wrong
  // This option is only available if Ipopt has been compiled with Pardiso.
  getStringInPList(pvApiCtx, param_in_addr, "pardiso_redo_symbolic_fact_only_if_inertia_wrong", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("pardiso_redo_symbolic_fact_only_if_inertia_wrong", tmp_char);
  FREE(tmp_char);
  // Interpretation of perturbed elements
  // no: Don't assume that matrix is singular if elements were perturbed after recent symbolic factorization
  // yes: Assume that matrix is singular if elements were perturbed after recent symbolic factorization
  // This option is only available if Ipopt has been compiled with Pardiso.
  getStringInPList(pvApiCtx, param_in_addr, "pardiso_repeated_perturbation_means_singular", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("pardiso_repeated_perturbation_means_singular", tmp_char);
  FREE(tmp_char);
  // Always pretent inertia is correct.
  // no: check inertia
  // yes: skip inertia check
  // Setting this option to "yes" essentially disables inertia check. This option makes the algorithm non-robust and easily fail, but it 
  // might give some insight into the necessity of inertia control.
  getStringInPList(pvApiCtx, param_in_addr, "pardiso_skip_inertia_check", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=-1) options->SetStringValue("pardiso_skip_inertia_check", tmp_char);
  FREE(tmp_char);
#endif

#ifdef USE_WSMP
  ////////////////////////
  // WSMP Linear Solver //
  ////////////////////////
  // Number of threads to be used in WSMP
  // This determines on how many processors WSMP is running on. This option is only available if Ipopt has been compiled with WSMP.
  getIntInPList(pvApiCtx, param_in_addr, "wsmp_num_threads", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("wsmp_num_threads", tmp_int);
  // Determines how ordering is done in WSMP
  // This corresponds to the value of WSSMP's IPARM(16). This option is only available if Ipopt has been compiled with WSMP.
  getIntInPList(pvApiCtx, param_in_addr, "wsmp_ordering_option", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("wsmp_ordering_option", tmp_int);
  // Pivot tolerance for the linear solver WSMP.
  // A smaller number pivots for sparsity, a larger number pivots for stability.  This option is only available if Ipopt has been compiled with WSMP.
  getDoubleInPList(pvApiCtx, param_in_addr, "wsmp_pivtol", &tmp_int, &tmp_res, 1e-4, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("wsmp_pivtol", tmp_double);
  // Maximum pivot tolerance for the linear solver WSMP.
  // Ipopt may increase pivtol as high as pivtolmax to get a more accurate solution to the linear system. This option is only available if Ipopt 
  // has been compiled with WSMP.
  getDoubleInPList(pvApiCtx, param_in_addr, "wsmp_pivtolmax", &tmp_int, &tmp_res, 0.1, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("wsmp_pivtolmax", tmp_double);
  // Determines how the matrix is scaled by WSMP.
  // This corresponds to the value of WSSMP's IPARM(10). This option is only available if Ipopt has been compiled with WSMP.
  getIntInPList(pvApiCtx, param_in_addr, "wsmp_scaling", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("wsmp_scaling", tmp_int);
#endif

  // wantsol: solution report without -AMPL: sum of
  // 1 ==> write .sol file
  // 2 ==> print primal variable values
  // 4 ==> print dual variable values
  // 8 ==> do not print solution message
  getIntInPList(pvApiCtx, param_in_addr, "wantsol", &tmp_int, &tmp_res, 8, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("wantsol", tmp_int);

  // print_user_options: If selected, the algorithm will print the list of all options set by
  //                     the user including their values and whether they have been used.  In
  //                     some cases this information might be incorrect, due to the internal program flow.
  // - no, don't print options (default value)
  // - yes, print options
  getStringInPList(pvApiCtx, param_in_addr, "print_user_options", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");

  if (tmp_res!=-1) options->SetStringValue("print_user_options", tmp_char);
  FREE(tmp_char);

  // print_options_documentation: If selected, the algorithm will print the list of all available
  //                              algorithmic options with some documentation before solving the
  //                              optimization problem.
  // - no, don't print list (default value)
  // - yes, print list

  getStringInPList(pvApiCtx, param_in_addr, "print_options_documentation", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");

  if (tmp_res!=-1) options->SetStringValue("print_options_documentation", tmp_char);
  FREE(tmp_char);

  return 0;
}
