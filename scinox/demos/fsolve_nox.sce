// the pipe test problem on the scilab fsolve_nox function

lines(0);
funcprot_old = funcprot();
funcprot(0);

global jac_ev;
global fun_ev;

jac_ev = 0;
fun_ev = 0;

/////////////////////
// Test from Petsc //
/////////////////////

function [F,J] = petsc_pb_1(x)
  global jac_ev;
  jac_ev = jac_ev + 1;
  global fun_ev;
  fun_ev = fun_ev + 1;
  
  F(1,1) = sin(3.0*x(1)) + x(1);
  F(2,1) = x(2);
  
  J(1,1) = 3*cos(3*x(1)) + 1;
  J(1,2) = 0;
  J(2,1) = 0;
  J(2,2) = 1;
endfunction

////////////////////////
// Rosenbrock problem //
////////////////////////

function [F,J] = genrosenbrock(x)
  global jac_ev;
  jac_ev = jac_ev + 1;
  global fun_ev;
  fun_ev = fun_ev + 1;
  // Derived from the generalized n-dimensional 
  // Rosenbrock function, returns a dense Jacobian 
  n = length(x);         
  if n == 0, error('Input vector, x, is empty.'); end
  if modulo(n,2) ~= 0, 
    error('Input vector, x ,must have an even number of components.');
  end
  // Evaluate the vector function
  odds  = 1:2:n;
  evens = 2:2:n;
  F = zeros(n,1);
  F(odds,1)  = 1-x(odds);
  F(evens,1) = 10.*(x(evens)-x(odds).^2); 
  // Evaluate the Jacobian matrix if nargout > 1
  [nargout,nargin] = argn();
  if nargout > 1
    c = -ones(n/2,1);    C = sparse([odds; odds]',c,[n,n]);
    d = 10*ones(n/2,1);  D = sparse([evens; evens]',d,[n,n]);
    e = -20.*x(odds);    E = sparse([evens; odds]',e,[n,n]);
    J = C + D + E;
  end
endfunction


function [f,J] = genrosenbrockdense(x)
  [nargout,nargin] = argn();
   
  if nargout == 1 
    f = genrosenbrock(x);
  end
  
  if nargout > 1 
    [f,J] = genrosenbrock(x);
    //J = derivative(genrosenbrock,x);
    J = full(J);
  end
endfunction

////////////////////////////
// Set the solver options //
////////////////////////////

params_nox = init_param();


/////////////////////////////////////
// Set the nonlinear solver method //
/////////////////////////////////////

// 'non_linear_solver':
// - 1 'Line Search Based'
// - 2 'Trust Region Based'
// - 3 'Inexact Trust Region Based'
// - 4 'Tensor Based'
params_nox = add_param(params_nox,'non_linear_solver',1);

//////////////////////////////////////////////////
// Parameters for the Trust Region Based Method //
//////////////////////////////////////////////////

// 'ns_trb_minimum_trust_region_radius':
// ($\Delta_{\min}$) - Minimum allowable trust region radius. Defaults to 1.0e-6.
params_nox = add_param(params_nox,'ns_trb_minimum_trust_region_radius',1.0e-6);

// 'ns_trb_maximum_trust_region_radius':
// ($\Delta_{\max}$) - Maximum allowable trust region radius. Defaults to 1.0e+10.
params_nox = add_param(params_nox,'ns_trb_maximum_trust_region_radius',1.0e10);

// 'ns_trb_minimum_improvement_ratio':
// ($\rho_{\min}$) - Minimum improvement ratio to accept the step. Defaults to 1.0e-4.
params_nox = add_param(params_nox,'ns_trb_minimum_improvement_ratio',1.0e-4);

// 'ns_trb_contraction_trigger_ratio':
// ($\rho_{\rm s}$) - If the improvement ratio is less than this value, then the trust region
// is contracted by the amount specified by the 'ns_trb_contraction_factor'. Must be larger than 'ns_trb_minimum_improvement_ratio'. Defaults to 0.1.
params_nox = add_param(params_nox,'ns_trb_contraction_trigger_ratio',0.1);

// 'ns_trb_contraction_factor':
// ($\beta_{\rm s}$) - See above. Defaults to 0.25.
params_nox = add_param(params_nox,'ns_trb_contraction_factor',0.25);

// 'ns_trb_expansion_trigger_ratio':
// ($\rho_{\rm e}$) - If the improvement ratio is greater than this value, then the trust region is 
// contracted by the amount specified by the 'Expansion Factor'. Defaults to 0.75.
params_nox = add_param(params_nox,'ns_trb_expansion_trigger_ratio',0.75);

// 'ns_trb_expansion_factor':
// ($\beta_{\rm e}$) - See above. Defaults to 4.0.
params_nox = add_param(params_nox,'ns_trb_expansion_factor',4.0);

// 'ns_trb_recovery_step':
// Defaults to 1.0.
params_nox = add_param(params_nox,'ns_trb_recovery_step',1.0);

// 'ns_trb_use_ared_pred_ratio_calculation':
// Defaults to false. If set to true, this option replaces the algorithm used to compute 
// the improvement ratio, $ \rho $, as described above. The improvement ratio is replaced by an 'Ared/Pred' sufficient decrease
// criteria similar to that used in line search algorithms 
// (see Eisenstat and Walker, SIAM Journal on Optimization V4 no. 2 (1994) pp 393-422):
params_nox = add_param(params_nox,'ns_trb_use_ared_pred_ratio_calculation',0);

//////////////////////////////////////////////////////////
// Parameters for the Inexact Trust Region Based Method //
//////////////////////////////////////////////////////////

// 'ns_itrb_inner_iteration_method':
// Choice of trust region algorithm to use. Choices are:
// o 1 - 'Standard Trust Region'
// o 2 - 'Inexact Trust Region'
params_nox = add_param(params_nox,'ns_itrb_inner_iteration_method',1);

// 'ns_itrb_minimum_trust_region_radius':
// ($\Delta_{\min}$) - Minimum allowable trust region radius. Defaults to 1.0e-6.
params_nox = add_param(params_nox,'ns_itrb_minimum_trust_region_radius',1.0e-6);

// 'ns_itrb_maximum_trust_region_radius':
// ($\Delta_{\max}$) - Minimum allowable trust region radius. Defaults to 1.0e+10.
params_nox = add_param(params_nox,'ns_itrb_maximum_trust_region_radius',1.0e10);

// 'ns_itrb_minimum_improvement_ratio':
// ($\rho_{\min}$) - Minimum improvement ratio to accept the step. Defaults to 1.0e-4.
params_nox = add_param(params_nox,'ns_itrb_minimum_improvement_ratio',1.0e-4);

// 'ns_itrb_contraction_trigger_ratio':
// ($\rho_{\rm s}$) - If the improvement ratio is less than this value, then the trust region is
// contracted by the amount specified by the 'ns_itrb_contraction_factor'. Must be larger than 'ns_itrb_minimum_improvement_ratio'. Defaults to 0.1.
params_nox = add_param(params_nox,'ns_itrb_contraction_trigger_ratio',0.1);

// 'ns_itrb_contraction_factor':
// ($\beta_{\rm s}$) - See above. Defaults to 0.25.
params_nox = add_param(params_nox,'ns_itrb_contraction_factor',0.25);

// 'ns_itrb_expansion_trigger_ratio':
// ($\rho_{\rm e}$) - If the improvement ratio is greater than this value, then the trust region is 
// contracted by the amount specified by the 'ns_itrb_expansion_factor'. Defaults to 0.75.
params_nox = add_param(params_nox,'ns_itrb_expansion_trigger_ratio',0.75);

// 'ns_itrb_expansion_factor':
// ($\beta_{\rm e}$) - See above. Defaults to 4.0.
params_nox = add_param(params_nox,'ns_itrb_expansion_factor',4.0);

// 'ns_itrb_recovery_step':
// Defaults to 1.0.
params_nox = add_param(params_nox,'ns_itrb_recovery_step',1.0);

// 'ns_itrb_use_ared_pred_ratio_calculation':
// Defaults to false. If set to true, this option replaces the algorithm used to
// compute the improvement ratio, $ \rho $, as described above. The improvement ratio is replaced by an 'Ared/Pred' sufficient decrease
// criteria similar to that used in line search algorithms (see Eisenstat and Walker, SIAM Journal on Optimization V4 no. 2 (1994) pp 393-422):
params_nox = add_param(params_nox,'ns_itrb_use_ared_pred_ratio_calculation',0);

// 'ns_itrb_use_cauchy_in_newton_direction':
// Used only by the 'Inexact Trust Region' algorithm. If set to true, the initial guess
// for the Newton direction computation will use the Cauchy direction as the initial guess. Defaults to false.
params_nox = add_param(params_nox,'ns_itrb_use_cauchy_in_newton_direction',0);

// 'ns_itrb_use_dogleg_segment_minimization':
// Used only by the 'Inexact Trust Region' algorithm. If set to true, 
// the $ \tau $ parameter is minimized over the dogleg line segments instead of being computed at the trust regioin radius. 
// Used only by the 'Inexact Trust Region' algorithm. Defaults to false.
params_nox = add_param(params_nox,'ns_itrb_use_dogleg_segment_minimization',0);

////////////////////////////////////////////
// Parameters for the Tensor Based Method //
////////////////////////////////////////////

// 'ns_tb_direction_method':
// Name of the direction to be computed in this solver. 'Tensor' and 'Newton' are the only two valid choices. 
//   A sublist by this name specifies all of the parameters to be passed to the linear solver. See below under 'Linear Solver'.
// - 1 - 'Tensor'
// - 2 - 'Newton'
params_nox = add_param(params_nox,'ns_tb_direction_method',2);

// 'ns_tb_rescue_bad_newton_solve':
// If the linear solve does not meet the tolerance specified by the forcing term, 
// then use the step anyway. Defaults to true.
params_nox = add_param(params_nox,'ns_tb_rescue_bad_newton_solve',1);

// Below is a partial list of standard parameters usually available in common linear solvers. 
// Check with the specific linear solver being used for other parameters.
//  o 'ns_tb_linear_solver_max_iterations' - Maximum number of Arnoldi iterations (also max Krylov space dimension)
//  o 'ns_tb_linear_solver_tolerance' - Relative tolerance for solving local model [default = 1e-4]
//  o 'ns_tb_linear_solver_output_frequency' - Print output at every number of iterations [default = 20]
params_nox = add_param(params_nox,'ns_tb_linear_solver_max_iterations',300);
params_nox = add_param(params_nox,'ns_tb_linear_solver_tolerance',1e-4);
params_nox = add_param(params_nox,'ns_tb_linear_solver_output_frequency',20);

// 'ns_tb_line_search_method':
// Because the tensor step is not guaranteed to be a descent direction on 
// the function, not all 'basic' line search approaches would be appropriate. Thus, the LineSearch classes available to Newton's
// method (e.g., Polynomial, More-Thuente) are not used here. Instead, this solver class approriately handles technical considerations
// for tensor methods with its own set of global strategies. The following parameters specify the specific options for this line search:
// Valid choices are:
//  o 1 - 'Curvilinear' - Backtrack along the 'curvilinear' path that spans the tensor direction and the Newton direction and
//         that maintains monotonicity on the tensor model. Recommended because it tends to be more robust and efficient than the
//         other choices. [Default]
//  o 2 - 'Standard' - Backtrack along tensor direction unless it is not a descent direction, in which case backtrack along Newton direction.
//  o 3 - 'Dual' - Backtrack along both the Newton and tensor directions and choose the better of the two.
//  o 4 - 'Full Step' - Only use the full step and do not backtrack along both the Newton and tensor directions and choose the better of the two.
params_nox = add_param(params_nox,'ns_tb_line_search_method',1);

// 'ns_tb_line_search_lambda_selection':
// Flag for how to calculate the next linesearch parameter lambda. Valid choices are 'Quadratic' and 
// 'Halving' (default). Quadratic constructs a quadratic interpolating polynomial from the last trial point and uses the minimum
// of this function as the next trial lambda (bounded by 0.1). Halving divides the linesearch parameter by 2 before each trial,
// which is simpler but tends to generate longer steps than quadratic.
params_nox = add_param(params_nox,'ns_tb_line_search_lambda_selection',2);

// 'ns_tb_default_step':
// Starting value of the linesearch parameter (defaults to 1.0)
params_nox = add_param(params_nox,'ns_tb_default_step',1.0);

// 'ns_tb_minimum_step':
// Minimum acceptable linesearch parameter before the linesearch terminates (defaults to 1.0e-12). 
// If there are many linesearch failures, then lowering this value is one thing to try.
params_nox = add_param(params_nox,'ns_tb_minimum_step',1.0e-12);

// 'ns_tb_recovery_step_type':
// Determines the step size to take when the line search fails. Choices are:
//  o 1 - 'Constant' [default] - Uses a constant value set in 'Recovery Step'.
//  o 2 - 'Last Computed Step' - Uses the last value computed by the line search algorithm.
params_nox = add_param(params_nox,'ns_tb_recovery_step_type',1);

// 'ns_tb_recovery_step':
// Step parameter to take when the line search fails (defaults to value for 'Default Step')
// - 1 - Default Step
// - 2 - Minimum Step
// - 3 - Recovery Step
params_nox = add_param(params_nox,'ns_tb_recovery_step',1);
   
// 'ns_tb_max_iters':
// Maximum number of iterations (i.e., backtracks)
params_nox = add_param(params_nox,'ns_tb_max_iters',300);

/////////////////////////
// Printing Parameters //
/////////////////////////

// Error 	Errors are always printed.
// Warning                  2^0
// OuterIteration           2^1
// InnerIteration           2^2
// Parameters               2^3
// Details                  2^4
// OuterIterationStatusTest 2^5
// LinearSolverDetails      2^6
// TestDetails              2^7
// Debug                    2^12

params_nox = add_param(params_nox,'printing_warning',0);
params_nox = add_param(params_nox,'printing_outeriterations',0);
params_nox = add_param(params_nox,'printing_inneriterations',0);
params_nox = add_param(params_nox,'printing_parameters',0);
params_nox = add_param(params_nox,'printing_details',0);
params_nox = add_param(params_nox,'printing_outeriterationstatustest',0);
params_nox = add_param(params_nox,'printing_linearsolverdetails',0);
params_nox = add_param(params_nox,'printing_testdetails',0);
params_nox = add_param(params_nox,'printing_debug',0);

/////////////////////////////////////////
// Selection of the line search Method //
/////////////////////////////////////////

// 'ls_method':
// - 1 'BackTrack'
// - 2 'Full Step'
// - 3 'Polynomial'
// - 4 'NonLinearCG'
// - 5 'More''-Thuente'

params_nox = add_param(params_nox,'ls_method',5);

////////////////////////////////////////////
// Parameters for the Backtracking Method //
////////////////////////////////////////////

// 'ls_bt_default_step':
// starting step length (defaults to 1.0) 
params_nox = add_param(params_nox,'ls_bt_default_step',1.0);

// 'ls_bt_minimum_step':
// minimum acceptable step length (defaults to 1.0e-12) 
params_nox = add_param(params_nox,'ls_bt_minimum_step',1.0e-12);

// 'ls_bt_recovery_step':
// step to take when the line search fails (defaults to value for 'Default Step') 
params_nox = add_param(params_nox,'ls_bt_recovery_step',1.0);

// 'ls_bt_max_iters':
// maximum number of iterations (i.e., RHS computations) 
params_nox = add_param(params_nox,'ls_bt_max_iters',20);

// 'ls_bt_decrease_condition': 
params_nox = add_param(params_nox,'ls_bt_decrease_condition',1);

// 'ls_bt_reduction_factor':
// A multiplier between zero and one that reduces the step size between line search iterations
params_nox = add_param(params_nox,'ls_bt_reduction_factor',0.5);

/////////////////////////////////////////
// Parameters for the Full Step Method //
/////////////////////////////////////////

// 'ls_fs_full_step':
// length of a full step (defaults to 1.0) 
params_nox = add_param(params_nox,'ls_fs_full_step',1.0);

//////////////////////////////////////////
// Parameters for the Polynomial Method //
//////////////////////////////////////////

// 'ls_pol_default_step':
// Starting step length, i.e., $\lambda_0$. Defaults to 1.0.
params_nox = add_param(params_nox,'ls_pol_default_step',1.0);

// 'ls_pol_max_iters':
// Maximum number of line search iterations. The search fails if the number of iterations exceeds this value. Defaults to 100.
params_nox = add_param(params_nox,'ls_pol_max_iters',100);

// 'ls_pol_minimum_step':
// Minimum acceptable step length. The search fails if the computed $\lambda_k$ is less than this value. Defaults to 1.0e-12.
params_nox = add_param(params_nox,'ls_pol_minimum_step',1.0e-12);

// 'ls_pol_recovery_step_type':
// Determines the step size to take when the line search fails. Choices are:
// - 1 'Constant' - Uses a constant value set in 'Recovery Step'. 
// - 2 'Last Computed Step' - Uses the last value computed by the line search algorithm. 
params_nox = add_param(params_nox,'ls_pol_recovery_step_type',1);

// 'ls_pol_recovery_step':
// The value of the step to take when the line search fails. Only used if the 'Recovery Step Type' is set to 'Constant'. 
// Defaults to value for 'Default Step'.
params_nox = add_param(params_nox,'ls_pol_recovery_step',1.0);

// 'ls_pol_interpolation_type':
// Type of interpolation that should be used.
params_nox = add_param(params_nox,'ls_pol_interpolation_type',1);

// 'ls_pol_min_bounds_factor':
// Choice for $ \gamma_{min} $, i.e., the factor that limits the minimum size of the new step based on the previous step. Defaults to 0.1.
params_nox = add_param(params_nox,'ls_pol_min_bounds_factor',0.1);

// 'ls_pol_max_bounds_factor':
// Choice for $ \gamma_{max} $, i.e., the factor that limits the maximum size of the new step based on the previous step. Defaults to 0.5.
params_nox = add_param(params_nox,'ls_pol_max_bounds_factor',0.5);

// 'ls_pol_sufficient_decrease_condition':
// The decrease condition
params_nox = add_param(params_nox,'ls_pol_sufficient_decrease_condition',2);

// 'ls_pol_alpha_factor':
// Parameter choice for sufficient decrease condition. See checkConvergence() for details. Defaults to 1.0e-4.
params_nox = add_param(params_nox,'ls_pol_alpha_factor',1.0e-4);

// 'ls_pol_force_interpolation':
// Set to true if at least one interpolation step should be used. The default is false which means that the line search will stop 
// if the default step length satisfies the convergence criteria. Defaults to false.
params_nox = add_param(params_nox,'ls_pol_force_interpolation',0);

// 'ls_pol_maximum_iteration_for_increase':
// Maximum index of the nonlinear iteration for which we allow a relative increase. See checkConvergence() 
// for further details. Defaults to 0 (zero).
params_nox = add_param(params_nox,'ls_pol_maximum_iteration_for_increase',0);

// 'ls_pol_allowed_relative_increase':
// See checkConvergence() for details. Defaults to 100.
params_nox = add_param(params_nox,'ls_pol_allowed_relative_increase',100);

/////////////////////////////////////////////
// Parameters for the More'-Thuente Method //
/////////////////////////////////////////////

// 'ls_mt_sufficient_decrease_condition':
// Choice to use for the sufficient decrease condition. Options are 'Ared/Pred' or 'Armijo-Goldstein' (defaults to 'Armijo-Goldstein').
// 1. 'Armijo-Goldstein' conditions: $ f(x_{n-1}+ \lambda s) \le f(x_{n-1}) +\alpha \lambda f'(x_{n-1}) $
// 2. 'Ared/Pred' conditions: 
//    $ \| F(x_{n-1}+ \lambda s) \| \le \| F(x_{n-1}) \| (1-\alpha(1-\eta)) 
//    $ where $ \eta $ is the linear solve tolerance in the inexact Newton method. 
params_nox = add_param(params_nox,'ls_mt_sufficient_decrease_condition',1);

// 'ls_mt_sufficient_decrease':
// The ftol in the sufficient decrease condition (defaults to 1.0e-4)
params_nox = add_param(params_nox,'ls_mt_sufficient_decrease',1.0e-4);

// 'ls_mt_curvature_condition':
// The gtol in the curvature condition (defaults to 0.9999)
params_nox = add_param(params_nox,'ls_mt_curvature_condition',0.9999);

// 'ls_mt_optimize_slope_computation':
// If set to true the value of $ s^TJ^TF $ is estimated using a directional derivative in a call to 
// NOX::LineSearch::Common::computeSlopeWithOutJac. If false the slope computation is computed with the 
// NOX::LineSearch::Common::computeSlope method. Setting this to true eliminates having to compute the Jacobian at each
// inner iteration of the More'-Thuente line search (defaults to false).
params_nox = add_param(params_nox,'ls_mt_optimize_slope_computation',0);

// 'ls_mt_interval_width':
// The maximum width of the interval containing the minimum of the modified function (defaults to 1.0e-15)
params_nox = add_param(params_nox,'ls_mt_interval_width',1.0e-15);

// 'ls_mt_maximum_step':
// maximum allowable step length (defaults to 1.0e6)
params_nox = add_param(params_nox,'ls_mt_maximum_step',1.0e6);

// 'ls_mt_minimum_step':
// minimum allowable step length (defaults to 1.0e-12)
params_nox = add_param(params_nox,'ls_mt_minimum_step',1.0e-12);

// 'ls_mt_max_iters':
// maximum number of right-hand-side and corresponding Jacobian evaluations (defaults to 20)
params_nox = add_param(params_nox,'ls_mt_max_iters',20);

// 'ls_mt_default_step':
// starting step length (defaults to 1.0)
params_nox = add_param(params_nox,'ls_mt_default_step',1.0);

// 'ls_mt_recovery_step_type':
// Determines the step size to take when the line search fails. Choices are:
// * 'Constant' [default] - Uses a constant value set in 'Recovery Step'.
// * 'Last Computed Step' - Uses the last value computed by the line search algorithm.
params_nox = add_param(params_nox,'ls_mt_recovery_step_type',1);

// 'ls_mt_recovery_step':
// The value of the step to take when the line search fails. Only used if the 'Recovery Step Type' is set to 'Constant'. 
// Defaults to value for 'Default Step'.
params_nox = add_param(params_nox,'ls_mt_recovery_step',1.0);

///////////////////////////////////////
// Selection of the direction Method //
///////////////////////////////////////

// 'dir_method':
// 1 - Broyden
// 2 - Newton
// 3 - Steepest Descent
// 4 - Nonlinear CG
params_nox = add_param(params_nox,'dir_method',1);

///////////////////////////////////////
// Parameters for the Broyden Method //
///////////////////////////////////////

// 'dir_broyden_restart_frequency':
// How often the Jacobian should be refreshed. A value of 5, for example, means that the Jacobian should be updated every 5 iterations.
// Defaults to 10.
params_nox = add_param(params_nox,'dir_broyden_restart_frequency',10);

// 'dir_broyden_max_convergence_rate':
// Maximum convergence rate allowed when reusing the Jacobian. The Jacobian will be refreshed if the convergence rate, $ \alpha $, 
// is larger than this value. The convergence rate is calculated by 
// $ \alpha = \frac{\| F_k \| }{\| F_{k-1} \|} 
// $ where F is the nonlinear residual and $ k $ is the nonlinear iteration. Defaults to 1.0.
params_nox = add_param(params_nox,'dir_broyden_max_convergence_rate',1.0);

// 'dir_broyden_memory':
// The maximum number of past updates that can be saved in memory. Defaults to the value of 'Restart Frequency'.
params_nox = add_param(params_nox,'dir_broyden_memory',10);

//////////////////////////////////////
// Parameters for the Newton Method //
//////////////////////////////////////

// 'dir_newton_forcing_term_method':
// Method to compute the forcing term, i.e., the tolerance for the linear solver.
// see http://trilinos.sandia.gov/packages/docs/r4.0/packages/nox/doc/html/classNOX_1_1Direction_1_1Newton.html
params_nox = add_param(params_nox,'dir_newton_forcing_term_method',1);

// 'dir_newton_forcing_term_initial_tolerance':
// $\eta_0$ (initial linear solver tolerance). Defaults to 0.1
params_nox = add_param(params_nox,'dir_newton_forcing_term_initial_tolerance',0.1);

// 'dir_newton_forcing_term_minimum_tolerance':
// $\eta_{\min}$. Defaults to 1.0e-6.
params_nox = add_param(params_nox,'dir_newton_forcing_term_minimum_tolerance',1.0e-6);

// 'dir_newton_forcing_term_maximum_tolerance':
// $\eta_{\max}$. Defaults to 0.01.
params_nox = add_param(params_nox,'dir_newton_forcing_term_maximum_tolerance',0.01);

// 'dir_newton_forcing_term_alpha':
// $\alpha$ (used only by 'Type 2'). Defaults to 1.5.
params_nox = add_param(params_nox,'dir_newton_forcing_term_alpha',1.5);

// 'dir_newton_forcing_term_gamma':
// $\gamma$ (used only by 'Type 2'). Defaults to 0.9
params_nox = add_param(params_nox,'dir_newton_forcing_term_gamma',1.5);

// 'dir_newton_forcing_term_method':
// see http://trilinos.sandia.gov/packages/docs/r6.0/packages/nox/doc/html/classNOX_1_1Direction_1_1SteepestDescent.html
params_nox = add_param(params_nox,'dir_newton_forcing_term_method',1);

////////////////////////////////////////////
// Parameters for the Nonlinear CG Method //
////////////////////////////////////////////

// 'dir_nonlinear_cg_restart_frequency':
// An integer specification of the number of nonlinear iterations between restarts [default = 10]. 
// Restart corresponds to setting $\beta = 0$. A good heuristic is to limit this value to the number of problem degrees of freedom. 
// Setting this value to 1 forces $ \beta = 0 $ for every nonlinear iteration which corresponds to suppressing orthogonalization 
// against the previous search direction.
params_nox = add_param(params_nox,'dir_nonlinear_cg_restart_frequency',10);

// 'dir_nonlinear_cg_precondition':
// can be either 'On' or 'Off' [default]: determines whether or not to compute and apply preconditioner $ M $. 
// If 'Off' is selected, no preconditioner is computed and the behavior is equivalent to $ M = I $ where $ I $ is the identity matrix. 
// If 'On', $ M $ is computed and applied as determined by the underlying implementation of the 'applyRightPreconditioning' method in the Group.
params_nox = add_param(params_nox,'dir_nonlinear_cg_precondition',1);

// 'dir_nonlinear_cg_orthogonalize':
// see http://trilinos.sandia.gov/packages/docs/r7.0/packages/nox/doc/html/classNOX_1_1Direction_1_1NonlinearCG.html
params_nox = add_param(params_nox,'dir_nonlinear_cg_orthogonalize',1);

//////////////////////////////
// Various other parameters //
//////////////////////////////

// 'newt_ls_aztec_solver':
// Determine the iterative technique used in the solve. The following options are valid:
// * 1 - 'GMRES'    - Restarted generalized minimal residual (default).
// * 2 - 'CG'       - Conjugate gradient.
// * 3 - 'CGS'      - Conjugate gradient squared.
// * 4 - 'TFQMR'    - Transpose-free quasi-minimal reasidual.
// * 5 - 'BiCGStab' - Bi-conjugate gradient with stabilization.
// * 6 - 'LU'       - Sparse direct solve (single processor only).
params_nox = add_param(params_nox,'newt_ls_aztec_solver',1);

// 'newt_ls_size_of_krylov_subspace':
// When using restarted GMRES this sets the maximum size of the Krylov subspace (defaults to 300).
params_nox = add_param(params_nox,'newt_ls_size_of_krylov_subspace',300);

// 'newt_ls_orthogonalization':
// The orthogonalization routine used for the Gram-Schmidt orthogonalization procedure in Aztec. The following options are valid:
// * 'Classical' - (default).
// * 'Modified'
params_nox = add_param(params_nox,'newt_ls_orthogonalization',1);

// 'newt_ls_convergence_test':
// Algorithm used to calculate the residual that is used for determining the convergence of the linear solver. 
// See the Aztec 2.1 manual for more information. The following options are valid:
// * 1 - 'r0' - (default)
// * 2 - 'rhs'
// * 3 - 'norm'
// * 4 - 'no scaling'
// * 5 - 'sol'
params_nox = add_param(params_nox,'newt_ls_convergence_test',1);

// 'newt_ls_tolerance':
// Tolerance used by AztecOO to determine if an iterative linear solve has converged.
params_nox = add_param(params_nox,'newt_ls_tolerance',1.0e-4);

// 'newt_ls_ill_conditioning_threshold':
// If the upper hessenberg matrix during GMRES generates a condition number greater than this parameter value, 
// aztec will exit the linear solve returning the it's current solution. The default is 1.0e11.
params_nox = add_param(params_nox,'newt_ls_ill_conditioning_threshold',1.0e11);

// 'newt_ls_preconditioner_iterations':
// Number of iterations an AztecOO_Operator should take when solving the preconditioner. 
// This is only used if an AztecOO preconditioner is used and the solver makes a call to NOX::Epetra::Group::applyRightPreconditioning(). 
// This is NOT a recomended approach.
params_nox = add_param(params_nox,'newt_ls_preconditioner_iterations',400);

// 'newt_ls_max_iterations':
// maximum number of iterations in the linear solve. Default is 400.
params_nox = add_param(params_nox,'newt_ls_max_iterations',400);

// 'newt_ls_zero_initial_guess':
// Zero out the initial guess for linear solves performed through applyJacobianInverse calls 
// (i.e. zero out the result vector before the linear solve). Defaults to false.
params_nox = add_param(params_nox,'newt_ls_zero_initial_guess',0);

// 'newt_ls_throw_error_on_prec_failure':
// If set to true, an exception will be thrown if the preconditioner fails to
// initialize or recompute/refactor. If set to false, a warning will br printed if the NOX::Utils::Warning is enabled 
// in the printing utilities (NOX::Utils). Defaults to true.
params_nox = add_param(params_nox,'newt_ls_throw_error_on_prec_failure',1);

// 'newt_ls_output_frequency':
// number of linear solve iterations between output of the linear solve residual. Takes an integer, or one of 
// the AztecOO flags: AZ_none, AZ_last, or AZ_all as a value. Defaults to AZ_last.
params_nox = add_param(params_nox,'newt_ls_output_frequency',2);

// 'newt_ls_jacobian_operator':
// If a constructor is used that does not supply a Jacobian operator, nox will create an internal Jacobian operator. 
// This flag is ONLY valid in such cases. This will determine which Operator is used:
// * 'Matrix-Free' - Create a NOX::Epetra::MatrixFree object.
// * 'Finite Difference' - Create a NOX::Epetra::FiniteDifference object.
params_nox = add_param(params_nox,'newt_ls_jacobian_operator',1);

// 'newt_ls_preconditioner':
// Sets the choice of the preconditioner to use during linear solves. The validity of the choice of preconditioner 
// will depend on the types of operators that are available for the Jacobian and preconditioner. 
// NOTE: This flag will override any constructor details. For example, if you supply a preconditioner operator in the constructor, 
// it will not be used if this flag is set to 'None'. If you supply an Epetra_Operator for the preconditioner but the 'Preconditioner' 
// flag is set to 'AztecOO' (this requires an Epetra_RowMatrix for the preconditioner operator), this object will exit with a failure. 
// The valid options and any requirements on the operator type are listed below:
// * 'None' - No preconditioning. (default)
// * 'AztecOO' - AztecOO internal preconditioner. This requires a preconditioner operator that derives from the Epetra_RowMatrix class.
// * 'Ifpack' - Ifpack internal preconditioner. This requires a preconditioner object that derives from the Epetra_RowMatrix class
//              or it can use a Jacobian if the Jacobian derives from an Epetra_RowMatrix. This option is deprecated. 
//              Please use 'New Ifpack'.
// * 'New Ifpack' - Ifpack internal preconditioner. This requires a preconditioner object that derives from the Epetra_RowMatrix 
//                  class or it can use a Jacobian if the Jacobian derives from an Epetra_RowMatrix.
params_nox = add_param(params_nox,'newt_ls_preconditioner',1);

// 'newt_ls_preconditioner_operator':
// If a constructor is used that does not supply a preconditioner operator, nox will create an internal
//  preconditioner operator. This flag is ONLY valid in such cases. This will determine which Operator is used:
// * 'Use Jacobian' - Use the Jacobian Operator (it must be an Epetra_RowMatrix derived object).
// * 'Finite Difference' - Create a NOX::Epetra::FiniteDifference object.
params_nox = add_param(params_nox,'newt_ls_preconditioner_operator',1);

// 'newt_ls_aztec_preconditioner':
// If the 'Preconditioner' flag is set to 'AztecOO' then the specific AztecOO preconditioner is specified with this flag. 
// Currently supported preconditioners and their corresponding parameters that can be set are shown below
// (See the Aztec 2.1 manual for more information):
// * 'ilu' - ilu preconditioning. This choice allows the following additional parameters to be specified:
//   o 'Overlap' - defaults to 0
//   o 'Graph Fill' - defaults to 0
// * 'ilut' - ilut preconditioning. This choice allows the following additional parameters to be specified:
//   o 'Overlap' - defaults to 0
//   o 'Fill Factor' - defaults to 1.0
//   o 'Drop Tolerance' - defaults to 1.0e-12
// * 'Jacobi' - k step Jacobi where k is set by the 'Steps' flag:
//   o 'Steps' - defaults to 3.
// * 'Symmetric Gauss-Siedel' - Non-overlapping domain decomposition k step symmetric Gauss-Siedel where k is set by the 'Steps' flag:
//   o 'Steps' - defaults to 3.
// * 'Polynomial' - Neumann polynomial with order set by the parameter:
//   o 'Polynomial Order' - defaults to 3.
// * 'Least-squares Polynomial' - Least-squares polynomial with order set by the parameter:
//   o 'Polynomial Order' - defaults to 3.
params_nox = add_param(params_nox,'newt_ls_aztec_preconditioner',1);

// 'newt_ls_preconditioner_reuse_policy':
// Allows the user to set how and when the preconditioner should be computed. 
// This flag supports native Aztec, Ifpack and ML preconditioners. There are three options:
// * 'Rebuild' - The 'Rebuild' option always completely destroys and then rebuilds the preconditioner each time a linear solve is requested.
// * 'Reuse' - The group/linear solver will not recompute the preconditioner even if the group's solution vector changes. 
//             It just blindly reuses what has been constructed. This turns off control of preconditioner recalculation. 
//             This is a dangerous condition but can really speed up the computations if the user knows what they are doing. 
//             We don't recommend users trying this.
// * 'Recompute' - Recomputes the preconditioner, but will try to efficiently reuse any objects that don't need to be destroyed. 
//                 How efficient the 'Recompute' option is depends on the type of preconditioner. For example if we are using ILU from the
//                 Ifpack library, we would like to not destroy and reallocate the graph each solve. With this option, we tell Ifpack to reuse
//                 the graph from last time - e.g the sparisty pattern has not changed between applications of the preconditioner.
params_nox = add_param(params_nox,'newt_ls_preconditioner_reuse_policy',1);

// 'newt_ls_max_age_of_prec':
// If the 'Preconditioner Reuse Policy' is set to 'Reuse', this integer tells the linear system how many 
// times to reuse the preconditioner before rebuilding it. Defaults to 1.
params_nox = add_param(params_nox,'newt_ls_max_age_of_prec',1);

// 'newt_ls_rcm_reordering':
// Enables RCM reordering in conjunction with domain decomp incomplete factorization preconditioning. The following options are valid:
// * 'Disabled' - (default).
// * 'Enabled'
params_nox = add_param(params_nox,'newt_ls_rcm_reordering',0);

// 'newt_ls_output_solver_details':
// Write the output sublist below to the parameter list after each linear solve. default is true.
params_nox = add_param(params_nox,'newt_ls_output_solver_details',1);

//////////////////////////////////////////////////
// Selection of the Jacobian computation method //
//////////////////////////////////////////////////

params_nox = add_param(params_nox,'use_analytic',1);
params_nox = add_param(params_nox,'use_matrix_free',1);

///////////////////////
// Finite Difference //
///////////////////////

params_nox = add_param(params_nox,'use_finite_difference',1);
params_nox = add_param(params_nox,'use_finite_difference_alpha',1.0e-4);
params_nox = add_param(params_nox,'use_finite_difference_beta',1.0e-6);
params_nox = add_param(params_nox,'use_finite_difference_type',1);

///////////////////////
// Convergence tests //
///////////////////////

// Various convergence tests based on the norm of the residual. the Absolute tolerance type. 
params_nox = add_param(params_nox,'test_absresid',1.0e-8);

// Various convergence tests based on the norm of the residual. the Relative tolerance type. 
params_nox = add_param(params_nox,'test_relresid',1.0e-2);

// Various convergence tests based on the norm of the change in the solution vector, x, between outer iteration
// the Absolute ToleranceType and TWO NormType. 
params_nox = add_param(params_nox,'test_update',1.0e-5);

// Convergence test based on the weighted root mean square norm fo the solution update between iterations.
params_nox = add_param(params_nox,'test_wrms_rtol',1.0e-5);
params_nox = add_param(params_nox,'test_wrms_atol',1.0e-5);

// Failure test based on the maximum number of nonlinear solver iterations. 
params_nox = add_param(params_nox,'test_maxiters',20);

// Failure test based on whether the norm of a vector has a finite value. 
params_nox = add_param(params_nox,'test_finite_value',0);

// Failure test based on a threshold value of the norm of F. 
params_nox = add_param(params_nox,'test_divergence',1.0e6);

// Failure test based on the convergence rate between nonlinear iterations. 
params_nox = add_param(params_nox,'test_stagnation_threshold',1.0);
params_nox = add_param(params_nox,'test_stagnation_iterations',50);

// Status of the solver:
// -2 Unevaluated Unevaluated.
//  0 Unconverged Neither Converged nor Failed.
//  1 Converged   Converged.
// -1 Failed      Failed. 

printf('begin to set params_nox\n');

params_nox = init_param();
params_nox = add_param(params_nox,'printing_warning',1);
params_nox = add_param(params_nox,'printing_outeriterations',1);
params_nox = add_param(params_nox,'printing_inneriterations',1);
params_nox = add_param(params_nox,'printing_parameters',1);
params_nox = add_param(params_nox,'printing_details',1);
params_nox = add_param(params_nox,'test_absresid',1.0e-8);
params_nox = add_param(params_nox,'test_relresid',1.0e-2);
params_nox = add_param(params_nox,'test_update',1.0e-5);
params_nox = add_param(params_nox,'test_wrms_rtol',1.0e-5);
params_nox = add_param(params_nox,'test_wrms_atol',1.0e-5);
params_nox = add_param(params_nox,'test_maxiters',200);
params_nox = add_param(params_nox,'test_finite_value',0);
params_nox = add_param(params_nox,'test_divergence',1.0e6);
params_nox = add_param(params_nox,'test_stagnation_threshold',1.0);
params_nox = add_param(params_nox,'test_stagnation_iterations',50);
//  params_nox = add_param(params_nox,'newton_linear_jacobian_operator',2);

//params_nox = add_param(params_nox,'non_linear_solver',1); // 1: line search based
//params_nox = add_param(params_nox,'dir_method',2);
//params_nox = add_param(params_nox,'ls_method',2);
//params_nox = add_param(params_nox,'newt_ls_aztec_solver',1); // 1: GMRES
//params_nox = add_param(params_nox,'newt_ls_output_frequency',3); // 3 - AZ_all

params_nox = add_param(params_nox,'use_analytic',1);
//params_nox = add_param(params_nox,'use_finite_difference',1);
//params_nox = add_param(params_nox,'use_finite_difference_type',1);

//params_nox = add_param(params_nox,'newt_ls_jacobian_operator',1);
//params_nox = add_param(params_nox,'newt_ls_compute_scaling_manually',0);
//params_nox = add_param(params_nox,'newt_ls_aztec_preconditioner',1);

////////////////
// Resolution //
////////////////

NbVar = 20; // 2 for petsc_pb_1
x0 = 4*ones(NbVar,1).*rand(NbVar,1) - 2*ones(NbVar,1);

// With analytical Jacobian and with fsolve
tic();
//[x_opt, status] = fsolver_nox(petsc_pb_1,x0,params_nox);
//[x_opt, status] = fsolver_nox(genrosenbrockdense,x0,params_nox);
[x_opt, status] = fsolver_nox(genrosenbrock,x0,params_nox);
t = toc();

printf('with Jacobian - result: \n');
printf('elapsed time = %f\n' ,t);
printf('number of Jacobian evaluation = %d\n', jac_ev);
printf('number of function evaluation = %d\n', fun_ev);
printf('info_nox     = %d\n', status('info_nox'));
printf('iterations   = %d\n', status('iterations'));

funcprot(funcprot_old);
