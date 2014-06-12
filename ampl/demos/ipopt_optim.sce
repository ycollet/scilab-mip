lines(0);

stacksize('max');

path = get_absolute_file_path('ipopt_optim.sce');

exec(path + 'nl_data.sce');

Solve_macminlp = %T; // 43  files
Solve_coinor   = %F; // 24  files
Solve_asl      = %F; // 4   files
Solve_modnl    = %F; // 898 files

UseDenseDg  = %F;
UseSparseDg = %T;

AddDeltaToX0 = %T;
///////////////////////////
// Set global paremeters //
///////////////////////////

nl_index = 36; // Bug avec bonmin_optim.sce en premier + pb 18 de macminlp; // CoinOR + Index 13  = hs100 // modnl: 52 - 62 - 462 - 489 ??
              // ModNL  + Index 440 = rosenbr // macminlp: 8 - 24 - 26 - [37 - 43]??
                
if Solve_macminlp then nl_filename = MacMINLP(nl_index); end
if Solve_coinor   then nl_filename = CoinOR(nl_index);   end
if Solve_asl      then nl_filename = ASL(nl_index);      end
if Solve_modnl    then nl_filename = ModNL(nl_index);    end

///////////////////////////////
// Load and test the problem //
///////////////////////////////

printf('\nOptimization of the %s problem.\n\n',basename(nl_filename));

[asl, x0, lower, upper, v, constr_lhs, constr_rhs] = ampl_init(nl_filename);

if AddDeltaToX0 then
  x0 = x0 + 0.1*ones(x0);
end
if length(constr_rhs)==0 then 
  printf('\nno constraints. Stops\n');
  return
end

if UseDenseDg then
  ///////////////////////////////////////////////
  // Dense objective and constraints functions //
  ///////////////////////////////////////////////

  deff('y=f(x)','[y,tmp] = ampl_evalf(asl,x);');
  deff('y=df(x)','[y,tmp] = ampl_evalg(asl,x);');

  deff('y=g(x)','[tmp,y] = ampl_evalf(asl,x);');
  deff('y=dg(x)','[tmp,y] = ampl_evalg(asl,x); ...
                  y = matrix(y,1,length(y));');

  // Defined the sparsity structure of the problem
  [tmp,dg0] = ampl_evalg(asl,x0);

  // For inequality constraints
  // Be careful, Scilab stores matrix by column.
  // A = [1 2 3; 4 5 6; 7 8 9; 10 11 12];
  // matrix(A,1,length(1)) =
  // ans  =
  // 
  //    1.    4.    7.    10.    2.    5.    8.    11.    3.    6.    9.    12.

  sparse_dg = [];
  Index = 1;
  for j=1:length(x0)
    for i=1:length(constr_rhs)
      sparse_dg(Index,1) = i;
      sparse_dg(Index,2) = j;
      Index = Index + 1;
    end
  end
end

if UseSparseDg then
  ////////////////////////////////////////////
  // Dense objective and sparse constraints //
  ////////////////////////////////////////////

  deff('y=f(x)','[y,tmp] = ampl_evalf(asl,x);');
  deff('y=df(x)','[y,tmp] = ampl_evalg(asl,x);');

  deff('y=g(x)','[tmp,y] = ampl_evalf(asl,x);');
  deff('y=dg(x)','y = ampl_eval_spst_g_val(asl,x);');

  [irow, jcol] = ampl_eval_spst_g_rc(asl,x0);

  sparse_dg = [];
  for i=1:length(irow)
    sparse_dg(i,1) = irow(i); // We have transposed dg
    sparse_dg(i,2) = jcol(i);
  end
end

// No Hessian informations
h  = [];
dh = [];
sparse_dh = [];

// Define type of variables and constraints
nb_constr = length(constr_rhs);
nb_var    = length(x0);

var_lin_type    = ones(nb_var,1);    // 0 Linear - 1 Non-Linear
constr_lin_type = ones(nb_constr,1); // 0 Linear - 1 Non-Linear
var_type = ampl_get_type(asl);
printf('\n');
printf('number of variables:   %d\n', nb_var);
printf('number of constraints: %d\n', nb_constr);
printf('variables type: %s\n', var_type);

////////////////////////////////////////////////////////////////////////

params = init_param();

/////////////
// Journal //
/////////////

// Journal verbosity level. 
// J_INSUPPRESSIBLE -1
// J_NONE 	          0
// J_ERROR 	         1
// J_STRONGWARNING   2
// J_SUMMARY 	       3
// J_WARNING 	       4
// J_ITERSUMMARY     5
// J_DETAILED        6
// J_MOREDETAILED    7
// J_VECTOR 	        8
// J_MOREVECTOR      9
// J_MATRIX 	        10
// J_MOREMATRIX      11
// J_ALL 	           12
// J_LAST_LEVEL      13
params = add_param(params,"journal_level",5);

////////////
// Output //
////////////

// Output verbosity level. 
// Sets the default verbosity level for console output. The larger this value the more detailed is the output. 
// The valid range for this integer option is -2 <= print_level <= 12 and its default value is 5.
params = add_param(params,"print_level",5);

// File name of options file (to overwrite default).
// By default, the name of the Ipopt options file is "ipopt.opt" - or something else if specified in the IpoptApplication::Initialize call. 
// If this option is set by SetStringValue BEFORE the options file is read, it specifies the name of the options file. 
// It does not make any sense to specify this option within the options file. The default value for this string option is "".
// Possible values:
// * Any acceptable standard file name
params = add_param(params,"option_file_name","ipopt.opt");

params = add_param(params,"jacobian_approximation", "finite-difference-values");
//params = add_param(params,"jacobian_approximation", "exact");

/////////////////
// Termination //
/////////////////

// Desired convergence tolerance (relative).
// Determines the convergence tolerance for the algorithm. The algorithm terminates successfully, if the (scaled) NLP error becomes smaller than this value,
// and if the (absolute) criteria according to "dual_inf_tol", "primal_inf_tol", and "cmpl_inf_tol" are met.
// (This is epsilon_tol in Eqn. (6) in implementation paper). See also "acceptable_tol" as a second termination criterion. 
// Note, some other algorithmic features also use this quantity to determine thresholds etc. 
// The valid range for this real option is 0 < tol < +inf and its default value is  1.10^-08. 
params = add_param(params,"tol",1e-8);

// Maximum number of iterations.
// The algorithm terminates with an error message if the number of iterations exceeded this number. The valid range for this integer option is 
// 0 <= max_iter < +inf and its default value is 3000
params = add_param(params,"max_iter",1000);

// Desired threshold for the dual infeasibility.
// Absolute tolerance on the dual infeasibility. Successful termination requires that the max-norm of the (unscaled) dual infeasibility is 
// less than this threshold. The valid range for this real option is 0 < dual_inf_tol < +inf and its default value is 1
params = add_param(params,"dual_inf_tol",1);

// Desired threshold for the constraint violation. 
// Absolute tolerance on the constraint violation. Successful termination requires that the max-norm of the (unscaled) constraint violation 
// is less than this threshold. The valid range for this real option is 0 < constr_viol_tol < +inf and its default value is 0.0001
params = add_param(params,"constr_viol_tol",1e-4);

// Desired threshold for the complementarity conditions.
// Absolute tolerance on the complementarity. Successful termination requires that the max-norm of the (unscaled) complementarity is less 
// than this threshold. The valid range for this real option is 0 < compl_inf_tol < +inf and its default value is 0.0001
params = add_param(params,"compl_inf_tol",1e-4);

// "Acceptable" convergence tolerance (relative).
// Determines which (scaled) overall optimality error is considered to be "acceptable." There are two levels of termination criteria. 
// If the usual "desired" tolerances (see tol, dual_inf_tol etc) are satisfied at an iteration, the algorithm immediately terminates with a success message. 
// On the other hand, if the algorithm encounters "acceptable_iter" many iterations in a row that are considered "acceptable", it will terminate before the 
// desired convergence tolerance is met. This is useful in cases where the algorithm might not be able to achieve the "desired" level of accuracy. 
// The valid range for this real option is 0 < acceptable_tol < +inf and its default value is 1.10^-6 
params = add_param(params,"acceptable_tol",1e-6);

//"Acceptance" threshold for the constraint violation.
// Absolute tolerance on the constraint violation. "Acceptable" termination requires that the max-norm of the (unscaled) constraint violation 
// is less than this threshold; see also acceptable_tol. 
// The valid range for this real option is 0 < acceptable_constr_viol_tol < +inf and its default value is 0.01
params = add_param(params,"acceptable_constr_viol_tol",1e-2);

// "Acceptance" threshold for the dual infeasibility.
// Absolute tolerance on the dual infeasibility. "Acceptable" termination requires that the (max-norm of the unscaled) dual infeasibility 
// is less than this threshold; see also acceptable_tol. 
// The valid range for this real option is 0 < acceptable_dual_inf_tol < +inf and its default value is 1.10^+10 
params = add_param(params,"acceptable_dual_inf_tol",1e10);

// "Acceptance" threshold for the complementarity conditions.
// Absolute tolerance on the complementarity. "Acceptable" termination requires that the max-norm of the (unscaled) complementarity is less 
// than this threshold; see also acceptable_tol. 
// The valid range for this real option is 0 < acceptable_compl_inf_tol < +inf and its default value is 0.01
params = add_param(params,"acceptable_compl_inf_tol",1e-2);

// Threshold for maximal value of primal iterates.
// If any component of the primal iterates exceeded this value (in absolute terms), the optimization is aborted with the exit message that 
// the iterates seem to be diverging. 
// The valid range for this real option is 0 < diverging_iterates_tol < +inf and its default value is 1.10^+20
params = add_param(params,"diverging_iterates_tol",1e20);

/////////////////
// NLP Scaling //
/////////////////

// Scaling factor for the objective function.
// This option sets a scaling factor for the objective function. The scaling is seen internally by Ipopt but the unscaled objective is reported 
// in the console output. If additional scaling parameters are computed (e.g. user-scaling or gradient-based), both factors are multiplied. 
// If this value is chosen to be negative, Ipopt will maximize the objective function instead of minimizing it. 
// The valid range for this real option is -inf < obj_scaling_factor < +inf and its default value is 1
params = add_param(params,"obj_scaling_factor",1);

// Select the technique used for scaling the NLP.
// Selects the technique used for scaling the problem internally before it is solved. For user-scaling, the parameters come from the NLP. 
// If you are using AMPL, they can be specified through suffixes ("scaling_factor") The default value for this string option is "gradient-based".
// Possible values:
//   * none: no problem scaling will be performed
//   * user-scaling: scaling parameters will come from the user
//   * gradient-based: scale the problem so the maximum gradient at the starting point is scaling_max_gradient
//   * equilibration-based: scale the problem so that first derivatives are of order 1 at random points (only available with MC19)
//params = add_param(params,"nlp_scaling_method","gradient-based");
params = add_param(params,"nlp_scaling_method","none");

// Maximum gradient after NLP scaling.
// This is the gradient scaling cut-off. If the maximum gradient is above this value, then gradient based scaling will be performed. 
// Scaling parameters are calculated to scale the maximum gradient back to this value. (This is g_max in Section 3.8 of the implementation paper.) 
// Note: This option is only used if "nlp_scaling_method" is chosen as "gradient-based". 
// The valid range for this real option is 0 < nlp_scaling_max_gradient < +inf and its default value is 100
params = add_param(params,"nlp_scaling_max_gradient",100);

/////////
// NLP //
/////////

// Factor for initial relaxation of the bounds.
// Before start of the optimization, the bounds given by the user are relaxed. This option sets the factor for this relaxation. 
// If it is set to zero, then then bounds relaxation is disabled. (See Eqn.(35) in implementation paper.) 
// The valid range for this real option is 0 <= bound_relax_factor < +inf and its default value is 1.10^-8 
params = add_param(params,"bound_relax_factor",1e-8);

// Indicates whether final points should be projected into original bounds.
// Ipopt might relax the bounds during the optimization (see, e.g., option "bound_relax_factor"). 
// This option determines whether the final point should be projected back into the user-provide original bounds after the optimization. 
// The default value for this string option is "yes".
// Possible values:
//  * no: Leave final point unchanged
//  * yes: Project final point back into original bounds
params = add_param(params,"honor_original_bounds","yes");

// Indicates whether it is desired to check for Nan/Inf in derivative matrices
// Activating this option will cause an error if an invalid number is detected in the constraint Jacobians or the Lagrangian Hessian. 
// If this is not activated, the test is skipped, and the algorithm might proceed with invalid numbers and fail. 
// The default value for this string option is "no".
// Possible values:
//  * no: Don't check (faster).
//  * yes: Check Jacobians and Hessian for Nan and Inf.
//params = add_param(params,"check_derivatives_for_naninf","no");
params = add_param(params,"check_derivatives_for_naninf","yes");

// any bound less or equal this value will be considered -inf (i.e. not lower bounded).
// The valid range for this real option is -inf < nlp_lower_bound_inf < +inf and its default value is -1.10^+19
params = add_param(params,"nlp_lower_bound_inf", -1e19);

// any bound greater or this value will be considered +inf (i.e. not upper bounded).
// The valid range for this real option is -inf < nlp_upper_bound_inf < +inf and its default value is 1.10^+19
params = add_param(params,"nlp_upper_bound_inf", 1e19);

// Determines how fixed variables should be handled.
// The main difference between those options is that the starting point in the "make_constraint" case still has the fixed variables at their 
// given values, whereas in the case "make_parameter" the functions are always evaluated with the fixed values for those variables. 
// Also, for "relax_bounds", the fixing bound constraints are relaxed (according to" bound_relax_factor"). 
// For both "make_constraints" and "relax_bounds", bound multipliers are computed for the fixed variables. 
// The default value for this string option is "make_parameter".
// Possible values:
//  * make_parameter: Remove fixed variable from optimization variables
//  * make_constraint: Add equality constraints fixing variables
//  * relax_bounds: Relax fixing bound constraints
params = add_param(params,"fixed_variable_treatment","make_parameter");

// Indicates whether all equality constraints are linear
// Activating this option will cause Ipopt to ask for the Jacobian of the equality constraints only once from the NLP and reuse this information later. 
// The default value for this string option is "no".
// Possible values:
//  * no: Don't assume that all equality constraints are linear
//  * yes: Assume that equality constraints Jacobian are constant
params = add_param(params,"jac_c_constant","no");

// Indicates whether all inequality constraints are linear
// Activating this option will cause Ipopt to ask for the Jacobian of the inequality constraints only once from the NLP and reuse this information later. 
// The default value for this string option is "no".
// Possible values:
//  * no: Don't assume that all inequality constraints are linear
//  * yes: Assume that equality constraints Jacobian are constant
params = add_param(params,"jac_d_constant","no");

// Indicates whether the problem is a quadratic problem 
// Activating this option will cause Ipopt to ask for the Hessian of the Lagrangian function only once from the NLP and reuse this information later. 
// The default value for this string option is "no".
// Possible values:
//  * no: Assume that Hessian changes
//  * yes: Assume that Hessian is constant
params = add_param(params,"hessian_constant","no");

////////////////////
// Initialization //
////////////////////

// Desired minimum relative distance from the initial point to bound.
// Determines how much the initial point might have to be modified in order to be sufficiently inside the bounds (together with "bound_push"). 
// (This is kappa_2 in Section 3.6 of implementation paper.) 
// The valid range for this real option is 0 < bound_frac <= 0.5 and its default value is 0.01
params = add_param(params,"bound_frac",1e-2);

// Desired minimum absolute distance from the initial point to bound.
// Determines how much the initial point might have to be modified in order to be sufficiently inside the bounds (together with "bound_frac"). 
// (This is kappa_1 in Section 3.6 of implementation paper.) The valid range for this real option is 0 < bound_push < +inf and its default value is 0.01
params = add_param(params,"bound_push",1e-2);

// Desired minimum relative distance from the initial slack to bound.
// Determines how much the initial slack variables might have to be modified in order to be sufficiently inside the inequality bounds 
// (together with "slack_bound_push"). (This is kappa_2 in Section 3.6 of implementation paper.) 
// The valid range for this real option is 0 < slack_bound_frac <= 0.5 and its default value is 0.01
params = add_param(params,"slack_bound_frac",1e-2);

// Desired minimum absolute distance from the initial slack to bound.
// Determines how much the initial slack variables might have to be modified in order to be sufficiently inside the inequality bounds 
// (together with "slack_bound_frac"). (This is kappa_1 in Section 3.6 of implementation paper.) 
// The valid range for this real option is 0 < slack_bound_push < +inf and its default value is 0.01
params = add_param(params,"slack_bound_push",1e-2);

// Initial value for the bound multipliers.
// All dual variables corresponding to bound constraints are initialized to this value. 
// The valid range for this real option is 0 < bound_mult_init_val < +inf and its default value is 1
params = add_param(params,"bound_mult_init_val",1);

// Maximum allowed least-square guess of constraint multipliers.
// Determines how large the initial least-square guesses of the constraint multipliers are allowed to be (in max-norm). 
// If the guess is larger than this value, it is discarded and all constraint multipliers are set to zero. 
// This options is also used when initializing the restoration phase. 
// By default, "resto.constr_mult_init_max" (the one used in RestoIterateInitializer) is set to zero. 
// The valid range for this real option is 0 <= constr_mult_init_max < +inf and its default value is 1000
params = add_param(params,"constr_mult_init_max",1000);

// Initialization method for bound multipliers
// This option defines how the iterates for the bound multipliers are initialized. If "constant" is chosen, then all bound multipliers 
// are initialized to the value of "bound_mult_init_val". If "mu-based" is chosen, the each value is initialized to the the value of "mu_init" 
// divided by the corresponding slack variable. This latter option might be useful if the starting point is close to the optimal solution. 
// The default value for this string option is "constant".
// Possible values:
//  * constant: set all bound multipliers to the value of bound_mult_init_val
//  * mu-based: initialize to mu_init/x_slack
params = add_param(params,"bound_mult_init_method","constant");

///////////////////////
// Barrier Parameter //
///////////////////////

// Indicates if we want to do Mehrotra's algorithm.
// If set to yes, Ipopt runs as Mehrotra's predictor-corrector algorithm. This works usually very well for LPs and convex QPs. 
// This automatically disables the line search, and chooses the (unglobalized) adaptive mu strategy with the "probing" oracle, 
// and uses "corrector_type=affine" without any safeguards; you should not set any of those options explicitly in addition. 
// Also, unless otherwise specified, the values of "bound_push", "bound_frac", and "bound_mult_init_val" are set more aggressive, 
// and sets "alpha_for_y=bound_mult". 
// The default value for this string option is "no".
// Possible values:
//  * no: Do the usual Ipopt algorithm.
//  * yes: Do Mehrotra's predictor-corrector algorithm.
params = add_param(params,"mehrotra_algorithm","no");

// Update strategy for barrier parameter. $ \;$
// Determines which barrier parameter update strategy is to be used. The default value for this string option is "monotone".
// Possible values:
//  * monotone: use the monotone (Fiacco-McCormick) strategy
//  * adaptive: use the adaptive update strategy
params = add_param(params,"mu_strategy","monotone");

// Oracle for a new barrier parameter in the adaptive strategy. 
// Determines how a new barrier parameter is computed in each "free-mode" iteration of the adaptive barrier parameter strategy. 
// (Only considered if "adaptive" is selected for option "mu_strategy"). 
// The default value for this string option is "quality-function".
// Possible values:
//  * probing: Mehrotra's probing heuristic
//  * loqo: LOQO's centrality rule
//  * quality-function: minimize a quality function
params = add_param(params,"mu_oracle","quality-function");

// Maximum number of search steps during direct search procedure determining the optimal centering parameter.
// The golden section search is performed for the quality function based mu oracle. (Only used if option "mu_oracle" is set to "quality-function".) 
// The valid range for this integer option is 0 <= quality_function_max_section_steps < +inf and its default value is 8
params = add_param(params,"quality_function_max_section_steps",8);

// Oracle for the barrier parameter when switching to fixed mode.
// Determines how the first value of the barrier parameter should be computed when switching to the "monotone mode" in the adaptive strategy. 
// (Only considered if "adaptive" is selected for option "mu_strategy".) 
// The default value for this string option is "average_compl".
// Possible values:
//  * probing: Mehrotra's probing heuristic
//  * loqo: LOQO's centrality rule
//  * quality-function: minimize a quality function
//  * average_compl: base on current average complementarity
params = add_param(params,"fixed_mu_oracle","average_compl");

// Initial value for the barrier parameter.
// This option determines the initial value for the barrier parameter (mu). It is only relevant in the monotone, Fiacco-McCormick version of the algorithm. 
// (i.e., if "mu_strategy" is chosen as "monotone") 
// The valid range for this real option is 0 < mu_init < +inf and its default value is 0.1 
params = add_param(params,"mu_init",0.1);

// Factor for initialization of maximum value for barrier parameter.
// This option determines the upper bound on the barrier parameter. This upper bound is computed as the average complementarity at the initial point 
// times the value of this option. (Only used if option "mu_strategy" is chosen as "adaptive".) 
// The valid range for this real option is 0 < mu_max_fact < +inf and its default value is 1000
params = add_param(params,"mu_max_fact",1000);

// Maximum value for barrier parameter. 
// This option specifies an upper bound on the barrier parameter in the adaptive mu selection mode. If this option is set, it overwrites 
// the effect of mu_max_fact. (Only used if option "mu_strategy" is chosen as "adaptive".) 
// The valid range for this real option is 0 < mu_max < +inf and its default value is 100000
params = add_param(params,"mu_max",100000);

// Minimum value for barrier parameter.
// This option specifies the lower bound on the barrier parameter in the adaptive mu selection mode. By default, it is set to the minimum 
// of 1e-11 and min("tol","compl_inf_tol")/("barrier_tol_fact- or"+1), which should be a reasonable value. (Only used if option "mu_strategy" 
// is chosen as "adaptive".) 
// The valid range for this real option is 0 < mu_min < +inf and its default value is 1.10^-11
params = add_param(params,"mu_min",1e-11);

// Factor for mu in barrier stop test.
// The convergence tolerance for each barrier problem in the monotone mode is the value of the barrier parameter times "barrier_tol_factor". 
// This option is also used in the adaptive mu strategy during the monotone mode. (This is kappa_epsilon in implementation paper). 
// The valid range for this real option is 0 < barrier_tol_factor < +inf and its default value is 10
params = add_param(params,"barrier_tol_factor",10);

// Determines linear decrease rate of barrier parameter. 
// For the Fiacco-McCormick update procedure the new barrier parameter mu is obtained by taking the minimum of mu*"mu_linear_decrease_factor" and 
// mu"superlinear_decrease_power". (This is kappa_mu in implementation paper.) This option is also used in the adaptive mu strategy during the monotone mode. 
// The valid range for this real option is 0 < mu_linear_decrease_factor < 1 and its default value is 0.2
params = add_param(params,"mu_linear_decrease_factor",0.2);

// Determines superlinear decrease rate of barrier parameter.
// For the Fiacco-McCormick update procedure the new barrier parameter mu is obtained by taking the minimum of mu*"mu_linear_decrease_factor" and 
// mu"superlinear_decrease_power". (This is theta_mu in implementation paper.) This option is also used in the adaptive mu strategy during the monotone mode. 
// The valid range for this real option is 1 < mu_superlinear_decrease_power < 2 and its default value is 1.5
params = add_param(params,"mu_superlinear_decrease_power",1.5);

////////////////////////
// Multiplier Updates //
////////////////////////

// Method to determine the step size for constraint multipliers. $ \;$
// This option determines how the step size (alpha_y) will be calculated when updating the constraint multipliers. The default value for this string option is "primal".
// Possible values:
//  * primal: use primal step size
//  * bound_mult: use step size for the bound multipliers (good for LPs)
//  * min: use the min of primal and bound multipliers
//  * max: use the max of primal and bound multipliers
//  * full: take a full step of size one
//  * min_dual_infeas: choose step size minimizing new dual infeasibility
//  * safe_min_dual_infeas: like "min_dual_infeas", but safeguarded by "min" and "max"
//  * primal-and-full: use the primal step size, and full step if delta_x <= alpha_for_y_tol
//  * dual-and-full: use the dual step size, and full step if delta_x <= alpha_for_y_tol
//  * acceptor: Call LSAcceptor to get step size for y
params = add_param(params,"alpha_for_y","primal");

// Tells the algorithm to recalculate the equality and inequality multipliers as least square estimates. $ \;$
// This asks the algorithm to recompute the multipliers, whenever the current infeasibility is less than recalc_y_feas_tol. Choosing yes might be helpful in the quasi-Newton option. However, each recalculation requires an extra factorization of the linear system. If a limited memory quasi-Newton option is chosen, this is used by default. The default value for this string option is "no".
// Possible values:
//  * no: use the Newton step to update the multipliers
//  * yes: use least-square multiplier estimates
params = add_param(params,"recalc_y","no");

// Feasibility threshold for recomputation of multipliers.
// If recalc_y is chosen and the current infeasibility is less than this value, then the multipliers are recomputed. 
// The valid range for this real option is 0 < recalc_y_feas_tol < +inf and its default value is 1.10^-6
params = add_param(params,"recalc_y_feas_tol",1e-6);

/////////////////
// Line Search //
/////////////////

// Maximum number of second order correction trial steps at each iteration.
// Choosing 0 disables the second order corrections. (This is pmax of Step A-5.9 of Algorithm A in the implementation paper.) 
// The valid range for this integer option is 0 <= max_soc < +inf and its default value is 4
params = add_param(params,"max_soc",4);

// Number of shortened iterations that trigger the watchdog.
// If the number of successive iterations in which the backtracking line search did not accept the first trial point exceeds this number, 
// the watchdog procedure is activated. Choosing "0" here disables the watchdog procedure. 
// The valid range for this integer option is 0 <= watchdog_shortened_iter_trigger < +inf and its default value is 10
params = add_param(params,"watchdog_shortened_iter_trigger",10);

// Maximum number of watchdog iterations.
// This option determines the number of trial iterations allowed before the watchdog procedure is aborted and the algorithm returns to the stored point. 
// The valid range for this integer option is 1 <= watchdog_trial_iter_max < +inf and its default value is 3
params = add_param(params,"watchdog_trial_iter_max",3);

// The type of corrector steps that should be taken (unsupported!).
// If "mu_strategy" is "adaptive", this option determines what kind of corrector steps should be tried. The default value for this string option is "none".
// Possible values:
//  * none: no corrector
//  * affine: corrector step towards mu=0
//  * primal-dual: corrector step towards current mu
params = add_param(params,"corrector_type","none");

////////////////
// Warm Start //
////////////////

// Warm-start for initial point
// Indicates whether this optimization should use a warm start initialization, where values of primal and dual variables are given 
// (e.g., from a previous optimization of a related problem.) 
// The default value for this string option is "no".
// Possible values:
//  * no: do not use the warm start initialization
//  * yes: use the warm start initialization
params = add_param(params,"warm_start_init_point","no");

// same as bound_push for the regular initializer.
// The valid range for this real option is 0 < warm_start_bound_push < +inf and its default value is 0.001
params = add_param(params,"warm_start_bound_push",1e-3);

// same as bound_frac for the regular initializer.
// The valid range for this real option is 0 < warm_start_bound_frac <= 0.5 and its default value is 0.001
params = add_param(params,"warm_start_bound_frac",1e-3);

// same as slack_bound_frac for the regular initializer.
// The valid range for this real option is 0 < warm_start_slack_bound_frac <= 0.5 and its default value is 0.001
params = add_param(params,"warm_start_slack_bound_frac",1e-3);

// same as slack_bound_push for the regular initializer.
// The valid range for this real option is 0 < warm_start_slack_bound_push < +inf and its default value is 0.001
params = add_param(params,"warm_start_slack_bound_push",1e-3);

// same as mult_bound_push for the regular initializer.
// The valid range for this real option is 0 < warm_start_mult_bound_push < +inf and its default value is 0.001
params = add_param(params,"warm_start_mult_bound_push",1e-3);

// Maximum initial value for the equality multipliers.
// The valid range for this real option is -inf < warm_start_mult_init_max < +inf and its default value is 1.10^+6 
params = add_param(params,"warm_start_mult_init_max",1e6);

///////////////////////
// Restoration Phase //
///////////////////////

// Enable heuristics to quickly detect an infeasible problem.
// This options is meant to activate heuristics that may speed up the infeasibility determination if you expect that there is a good chance for 
// the problem to be infeasible. In the filter line search procedure, the restoration phase is called more quickly than usually, and more reduction 
// in the constraint violation is enforced before the restoration phase is left. If the problem is square, this option is enabled automatically. 
// The default value for this string option is "no".
// Possible values:
//  * no: the problem probably be feasible
//  * yes: the problem has a good chance to be infeasible
params = add_param(params,"expect_infeasible_problem","no");
//params = add_param(params,"expect_infeasible_problem","yes");

// Threshold for disabling "expect_infeasible_problem" option.
// If the constraint violation becomes smaller than this threshold, the "expect_infeasible_problem" heuristics in the filter line search are disabled. 
// If the problem is square, this options is set to 0. 
// The valid range for this real option is 0 <= expect_infeasible_problem_ctol < +inf and its default value is 0.001
params = add_param(params,"expect_infeasible_problem_ctol",1e-3);

// Tells algorithm to switch to restoration phase in first iteration.
// Setting this option to "yes" forces the algorithm to switch to the feasibility restoration phase in the first iteration. 
// If the initial point is feasible, the algorithm will abort with a failure. 
// The default value for this string option is "no".
// Possible values:
//  * no: don't force start in restoration phase
//  * yes: force start in restoration phase
params = add_param(params,"start_with_resto","no");
//params = add_param(params,"start_with_resto","yes");

// Required reduction in primal-dual error in the soft restoration phase.
// The soft restoration phase attempts to reduce the primal-dual error with regular steps. If the damped primal-dual step (damped only to satisfy 
// the fraction-to-the-boundary rule) is not decreasing the primal-dual error by at least this factor, then the regular restoration phase is called. 
// Choosing "0" here disables the soft restoration phase. 
// The valid range for this real option is 0 <= soft_resto_pderror_reduction_factor < +inf and its default value is 0.9999
params = add_param(params,"soft_resto_pderror_reduction_factor",0.9999);

// Required reduction of infeasibility before leaving restoration phase.
// The restoration phase algorithm is performed, until a point is found that is acceptable to the filter and the infeasibility has been reduced by 
// at least the fraction given by this option. 
// The valid range for this real option is 0 <= required_infeasibility_reduction < 1 and its default value is 0.9
params = add_param(params,"required_infeasibility_reduction",0.9);

// Threshold for resetting bound multipliers after the restoration phase.
// After returning from the restoration phase, the bound multipliers are updated with a Newton step for complementarity. 
// Here, the change in the primal variables during the entire restoration phase is taken to be the corresponding primal Newton step. 
// However, if after the update the largest bound multiplier exceeds the threshold specified by this option, the multipliers are all reset to 1. 
// The valid range for this real option is 0 <= bound_mult_reset_threshold < +inf and its default value is 1000
params = add_param(params,"bound_mult_reset_threshold",1000);

// Threshold for resetting equality and inequality multipliers after restoration phase.
// After returning from the restoration phase, the constraint multipliers are recomputed by a least square estimate. 
// This option triggers when those least-square estimates should be ignored. 
// The valid range for this real option is 0 <= constr_mult_reset_threshold < +inf and its default value is 0. 
params = add_param(params,"constr_mult_reset_threshold",0);

// Determines if the original objective function should be evaluated at restoration phase trial points.
// Setting this option to "yes" makes the restoration phase algorithm evaluate the objective function of the original problem at every trial point 
// encountered during the restoration phase, even if this value is not required. In this way, it is guaranteed that the original objective function 
// can be evaluated without error at all accepted iterates; otherwise the algorithm might fail at a point where the restoration phase accepts an iterate 
// that is good for the restoration phase problem, but not the original problem. On the other hand, if the evaluation of the original objective is expensive, 
// this might be costly.
// The default value for this string option is "yes".
// Possible values:
//  * no: skip evaluation
//  * yes: evaluate at every trial point
params = add_param(params,"evaluate_orig_obj_at_resto_trial","yes");
//params = add_param(params,"evaluate_orig_obj_at_resto_trial","no");

///////////////////
// Linear Solver //
///////////////////

// Linear solver used for step computations. 
// Determines which linear algebra package is to be used for the solution of the augmented linear system (for obtaining the search directions). 
// Note, the code must have been compiled with the linear solver you want to choose. Depending on your Ipopt installation, not all options are available. 
// The default value for this string option is "ma27".
// Possible values:
//  * ma27: use the Harwell routine MA27
//  * ma57: use the Harwell routine MA57
//  * pardiso: use the Pardiso package
//  * wsmp: use WSMP package
//  * mumps: use MUMPS package
//  * custom: use custom linear solver
if MSDOS then
  params = add_param(params,"linear_solver","ma27"); // YC trouver avec quel solveur on a compile - normalement, ma27
else
  params = add_param(params,"linear_solver","mumps"); // YC: linux - pour l'instant 
end

// Method for scaling the linear system.
// Determines the method used to compute symmetric scaling factors for the augmented system (see also the "linear_scaling_on_demand" option). 
// This scaling is independent of the NLP problem scaling. By default, MC19 is only used if MA27 or MA57 are selected as linear solvers. 
// This option is only available if Ipopt has been compiled with MC19. 
// The default value for this string option is "mc19".
// Possible values:
//  * none: no scaling will be performed
//  * mc19: use the Harwell routine MC19
if MSDOS then
  params = add_param(params,"linear_system_scaling","mc19"); // YC: normalement mc19
else
  params = add_param(params,"linear_system_scaling","none"); // YC: linux - pour l'instant
end

// Flag indicating that linear scaling is only done if it seems required.
// This option is only important if a linear scaling method (e.g., mc19) is used. If you choose "no", then the scaling factors are computed for 
// every linear system from the start. This can be quite expensive. Choosing "yes" means that the algorithm will start the scaling method only 
// when the solutions to the linear system seem not good, and then use it until the end. 
// The default value for this string option is "yes".
// Possible values:
//  * no: Always scale the linear system.
//  * yes: Start using linear system scaling if solutions seem not good.
params = add_param(params,"linear_scaling_on_demand","yes");

// Maximum number of iterative refinement steps per linear system solve.
// Iterative refinement (on the full unsymmetric system) is performed for each right hand side. This option determines the maximum number of iterative 
// refinement steps. 
// The valid range for this integer option is 0 <= max_refinement_steps < +inf and its default value is 10 
params = add_param(params,"max_refinement_steps",10);

// Minimum number of iterative refinement steps per linear system solve.
// Iterative refinement (on the full unsymmetric system) is performed for each right hand side. This option determines the minimum number of iterative 
// refinements (i.e. at least "min_refinement_steps" iterative refinement steps are enforced per right hand side.) 
// The valid range for this integer option is 0 <= min_refinement_steps < +inf and its default value is 1
params = add_param(params,"min_refinement_steps",1);

//////////////////////////
// Hessian Perturbation //
//////////////////////////

// Maximum value of regularization parameter for handling negative curvature.
// In order to guarantee that the search directions are indeed proper descent directions, Ipopt requires that the inertia of the (augmented) linear system for 
// the step computation has the correct number of negative and positive eigenvalues. The idea is that this guides the algorithm away from maximizers and makes 
// Ipopt more likely converge to first order optimal points that are minimizers. If the inertia is not correct, a multiple of the identity matrix is added to 
// the Hessian of the Lagrangian in the augmented system. This parameter gives the maximum value of the regularization parameter. 
// If a regularization of that size is not enough, the algorithm skips this iteration and goes to the restoration phase. 
// (This is delta_wmax in the implementation paper.) 
// The valid range for this real option is 0 < max_hessian_perturbation < +inf and its default value is 1.10^+20 
params = add_param(params,"max_hessian_perturbation",1e20);

// Smallest perturbation of the Hessian block.
// The size of the perturbation of the Hessian block is never selected smaller than this value, unless no perturbation is necessary. 
// (This is delta_wmin in implementation paper.) 
// The valid range for this real option is 0 <= min_hessian_perturbation < +inf and its default value is 1.10^-20
params = add_param(params,"min_hessian_perturbation",1e-20);

// Size of first x-s perturbation tried.
// The first value tried for the x-s perturbation in the inertia correction scheme.(This is delta_0 in the implementation paper.) 
// The valid range for this real option is 0 < first_hessian_perturbation < +inf and its default value is 0.0001
params = add_param(params,"first_hessian_perturbation",1e-4);

// Increase factor for x-s perturbation for very first perturbation.
// The factor by which the perturbation is increased when a trial value was not sufficient - this value is used for the computation of the very first 
// perturbation and allows a different value for for the first perturbation than that used for the remaining perturbations. (This is bar_kappa_w+ in 
// the implementation paper.) 
// The valid range for this real option is 1 < perturb_inc_fact_first < +inf and its default value is 100
params = add_param(params,"perturb_inc_fact_first",100);

// Increase factor for x-s perturbation.
// The factor by which the perturbation is increased when a trial value was not sufficient - this value is used for the computation of all perturbations 
// except for the first. (This is kappa_w+ in the implementation paper.) 
// The valid range for this real option is 1 < perturb_inc_fact < +inf and its default value is 8
params = add_param(params,"perturb_inc_fact",8);

// Decrease factor for x-s perturbation.
// The factor by which the perturbation is decreased when a trial value is deduced from the size of the most recent successful perturbation. 
// (This is kappa_w- in the implementation paper.) 
// The valid range for this real option is 0 < perturb_dec_fact < 1 and its default value is 0.333333
params = add_param(params,"perturb_dec_fact",0.333333);

// Size of the regularization for rank-deficient constraint Jacobians.
// (This is bar delta_c in the implementation paper.)
// The valid range for this real option is 0 <= jacobian_regularization_value < +inf and its default value is 1.10^-8 
params = add_param(params,"jacobian_regularization_value",1e-8);

//////////////////
// Quasi-Newton //
//////////////////

// Indicates what Hessian information is to be used.
// This determines which kind of information for the Hessian of the Lagrangian function is used by the algorithm. 
// The default value for this string option is "exact".
// Possible values:
//  * exact: Use second derivatives provided by the NLP.
//  * limited-memory: Perform a limited-memory quasi-Newton approximation
//params = add_param(params,"hessian_approximation","exact");
params = add_param(params,"hessian_approximation","limited-memory");

// Maximum size of the history for the limited quasi-Newton Hessian approximation.
// This option determines the number of most recent iterations that are taken into account for the limited-memory quasi-Newton approximation. 
// The valid range for this integer option is 0 <= limited_memory_max_history < +inf and its default value is 6 
params = add_param(params,"limited_memory_max_history",6);

// Threshold for successive iterations where update is skipped.
// If the update is skipped more than this number of successive iterations, we quasi-Newton approximation is reset. *
// The valid range for this integer option is 1 <= limited_memory_max_skipping < +inf and its default value is 2
params = add_param(params,"limited_memory_max_skipping",2);

/////////////////////
// Derivative Test //
/////////////////////

// Enable derivative checker
// If this option is enabled, a (slow) derivative test will be performed before the optimization. The test is performed at the user provided 
// starting point and marks derivative values that seem suspicious 
// The default value for this string option is "none".
// Possible values:
//  * none: do not perform derivative test
//  * first-order: perform test of first derivatives at starting point
//  * second-order: perform test of first and second derivatives at starting point
//params = add_param(params,"derivative_test","none");
params = add_param(params,"derivative_test","first-order");

// Size of the finite difference perturbation in derivative test.
// This determines the relative perturbation of the variable entries. 
// The valid range for this real option is 0 < derivative_test_perturbation < +inf and its default value is 1.10^-8. 
params = add_param(params,"derivative_test_perturbation",1e-8);

// Threshold for indicating wrong derivative.
// If the relative deviation of the estimated derivative from the given one is larger than this value, the corresponding derivative is marked as wrong. 
// The valid range for this real option is 0 < derivative_test_tol < +inf and its default value is 0.0001
params = add_param(params,"derivative_test_tol",1e-4);

// Indicates whether information for all estimated derivatives should be printed.
// Determines verbosity of derivative checker. 
// The default value for this string option is "no".
// Possible values:
//  * no: Print only suspect derivatives
//  * yes: Print all derivatives
params = add_param(params,"derivative_test_print_all","no");

// Maximal perturbation of an evaluation point.
// If a random perturbation of a points is required, this number indicates the maximal perturbation. 
// This is for example used when determining the center point at which the finite difference derivative test is executed. 
// The valid range for this real option is 0 <= point_perturbation_radius < +inf and its default value is 10
params = add_param(params,"point_perturbation_radius",10);

////////////////////////
// MA27 Linear Solver //
////////////////////////
// Pivot tolerance for the linear solver MA27.
// A smaller number pivots for sparsity, a larger number pivots for stability. This option is only available if Ipopt has been compiled with MA27. 
// The valid range for this real option is 0 < ma27_pivtol < 1 and its default value is 1.10^-8 
params = add_param(params,"ma27_pivtol",1e-8);

// Maximum pivot tolerance for the linear solver MA27.
// Ipopt may increase pivtol as high as pivtolmax to get a more accurate solution to the linear system. 
// This option is only available if Ipopt has been compiled with MA27. 
// The valid range for this real option is 0 < ma27_pivtolmax < 1 and its default value is 0.0001
params = add_param(params,"ma27_pivtolmax",1e-4);

// Integer workspace memory for MA27.
// The initial integer workspace memory = liw_init_factor * memory required by unfactored system. Ipopt will increase the workspace size by meminc_factor 
// if required. This option is only available if Ipopt has been compiled with MA27. 
// The valid range for this real option is 1 <= ma27_liw_init_factor < +inf and its default value is 5
params = add_param(params,"ma27_liw_init_factor",5);

// Real workspace memory for MA27.
// The initial real workspace memory = la_init_factor * memory required by unfactored system. Ipopt will increase the workspace size by meminc_factor 
// if required. This option is only available if Ipopt has been compiled with MA27. 
// The valid range for this real option is 1 <= ma27_la_init_factor < +inf and its default value is 5
params = add_param(params,"ma27_la_init_factor",5);

// Increment factor for workspace size for MA27.
// If the integer or real workspace is not large enough, Ipopt will increase its size by this factor. This option is only available if 
// Ipopt has been compiled with MA27. 
// The valid range for this real option is 1 <= ma27_meminc_factor < +inf and its default value is 10
params = add_param(params,"ma27_meminc_factor",5);

////////////////////////
// MA57 Linear Solver //
////////////////////////

// Pivot tolerance for the linear solver MA57.
// A smaller number pivots for sparsity, a larger number pivots for stability. This option is only available if Ipopt has been compiled with MA57. 
// The valid range for this real option is 0 < ma57_pivtol < 1 and its default value is 1.10^-8 
//params = add_param(params,"ma57_pivtol",1e-8);

// Maximum pivot tolerance for the linear solver MA57.
// Ipopt may increase pivtol as high as ma57_pivtolmax to get a more accurate solution to the linear system. 
// This option is only available if Ipopt has been compiled with MA57. 
// The valid range for this real option is 0 < ma57_pivtolmax < 1 and its default value is 0.0001
//params = add_param(params,"ma57_pivtolmax",1e-4);

// Safety factor for work space memory allocation for the linear solver MA57.
// If 1 is chosen, the suggested amount of work space is used. However, choosing a larger number might avoid reallocation if the suggest values 
// do not suffice. This option is only available if Ipopt has been compiled with MA57. 
// The valid range for this real option is 1 <= ma57_pre_alloc < +inf and its default value is 3
//params = add_param(params,"ma57_pre_alloc",3);

/////////////////////////
// MUMPS Linear Solver //
/////////////////////////

// Pivot tolerance for the linear solver MUMPS.
// A smaller number pivots for sparsity, a larger number pivots for stability. This option is only available if Ipopt has been compiled with MUMPS. 
// The valid range for this real option is 0 <= mumps_pivtol <= 1 and its default value is 1.10^-6 
//params = add_param(params,"mumps_pivtol",1e-6);

// Maximum pivot tolerance for the linear solver MUMPS.
// Ipopt may increase pivtol as high as pivtolmax to get a more accurate solution to the linear system. This option is only available if Ipopt 
// has been compiled with MUMPS. 
// The valid range for this real option is 0 <= mumps_pivtolmax <= 1 and its default value is 0.1
//params = add_param(params,"mumps_pivtolmax",0.1);

// Percentage increase in the estimated working space for MUMPS.
// In MUMPS when significant extra fill-in is caused by numerical pivoting, larger values of mumps_mem_percent may help use the workspace more efficiently. 
// The valid range for this integer option is 0 <= mumps_mem_percent < +inf and its default value is 1000
//params = add_param(params,"mumps_mem_percent",1000);

// Controls permuting and scaling in MUMPS
// This is ICTL(6) in MUMPS. 
// The valid range for this integer option is 0 <= mumps_permuting_scaling <= 7 and its default value is 7
//params = add_param(params,"mumps_permuting_scaling",7);

// Controls pivot order in MUMPS
// This is ICTL(7) in MUMPS.
// The valid range for this integer option is 0 <= mumps_pivot_order <= 7 and its default value is 7
//params = add_param(params,"mumps_pivot_order",7);

// Controls scaling in MUMPS
// This is ICTL(8) in MUMPS.
// The valid range for this integer option is -2 <= mumps_scaling <= 7 and its default value is 7
//params = add_param(params,"mumps_scaling",7);

///////////////////////////
// Pardiso Linear Solver //
///////////////////////////

// Matching strategy to be used by Pardiso
// This is IPAR(13) in Pardiso manual. This option is only available if Ipopt has been compiled with Pardiso. 
// The default value for this string option is "complete+2x2".
// Possible values:
//  * complete: Match complete (IPAR(13)=1)
//  * complete+2x2: Match complete+2x2 (IPAR(13)=2)
//  * constraints: Match constraints (IPAR(13)=3)
//params = add_param(params,"pardiso_matching_strategy","complete+2x2");

// Enables out-of-core variant of Pardiso
// Setting this option to a positive integer k makes Pardiso work in the out-of-core variant where the factor is split in 2k subdomains. 
// This is IPARM(50) in the Pardiso manual. This option is only available if Ipopt has been compiled with Pardiso. 
// The valid range for this integer option is 0 <= pardiso_out_of_core_power < +inf and its default value is 0
//params = add_param(params,"pardiso_out_of_core_power",0);

////////////////////////
// WSMP Linear Solver //
////////////////////////

// Number of threads to be used in WSMP
// This determines on how many processors WSMP is running on. This option is only available if Ipopt has been compiled with WSMP. 
// The valid range for this integer option is 1 <= wsmp_num_threads < +inf and its default value is 1
//params = add_param(params,"wsmp_num_threads",1);

// Determines how ordering is done in WSMP
// This corresponds to the value of WSSMP's IPARM(16). This option is only available if Ipopt has been compiled with WSMP. 
// The valid range for this integer option is -2 <= wsmp_ordering_option <= 3 and its default value is 1
//params = add_param(params,"wsmp_ordering_option",1);

// Pivot tolerance for the linear solver WSMP.
// A smaller number pivots for sparsity, a larger number pivots for stability. This option is only available if Ipopt has been compiled with WSMP. 
// The valid range for this real option is 0 < wsmp_pivtol < 1 and its default value is 0.0001
//params = add_param(params,"wsmp_pivtol",1e-4);

// Maximum pivot tolerance for the linear solver WSMP.
// Ipopt may increase pivtol as high as pivtolmax to get a more accurate solution to the linear system. 
// This option is only available if Ipopt has been compiled with WSMP. 
// The valid range for this real option is 0 < wsmp_pivtolmax < 1 and its default value is 0.1
//params = add_param(params,"wsmp_pivtolmax",0.1);

// Determines how the matrix is scaled by WSMP.
// This corresponds to the value of WSSMP's IPARM(10). This option is only available if Ipopt has been compiled with WSMP. 
// The valid range for this integer option is 0 <= wsmp_scaling <= 3 and its default value is 0
//params = add_param(params,"wsmp_scaling",0);

tic();
[x_sol, f_sol, extra] = ipopt(x0, f, df, nb_constr, g, dg, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);
t = toc();

// Value of extra('status'):
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
printf('status = %d\n',extra('status'));

printf('iteration count = %d\n', extra('it_count'));
printf('cpu time        = %f\n', extra('cpu_time'));
printf('number of objective function evaluation              = %d\n', extra('fobj_eval'));
printf('number of gradient of objective function evaluation  = %d\n', extra('fobj_grad_eval'));
printf('number of constraint function evaluation             = %d\n',extra('constr_eval'));
printf('number of gradient of constraint function evaluation = %d\n',extra('constr_jac_eval'));
printf('number of hessian function evaluation                = %d\n',extra('hess_eval'));
printf('dual infeasibility   = %f\n', extra('dual_inf'));
printf('constraint violation = %f\n', extra('constr_viol'));
printf('complementarity      = %f\n', extra('complementarity'));
printf('kkt error            = %f\n', extra('kkt_error'));
printf('time needed for resolution = %f\n',t);

printf('\nconstraints result:\n');
printf('lhs <= g(x_sol) <= rhs - lambda \n');

disp(string(constr_lhs) + ' <= ' + string(g(x_sol)) + ' <= ' + string(constr_rhs) + ' - lambda = ' + string(extra('lambda')'));

// Release the memory allocated for the selected AMPL problem
ampl_free(asl);
asl = [];

