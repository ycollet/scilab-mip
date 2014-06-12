lines(0);

stacksize('max');

path = get_absolute_file_path('bonmin_optim.sce');

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

// 26 - definition domain worng ?
// 27 - convergence to a local infeasible point
// 28 - Optim solution found (long)
// 29 - Relaxed pb not feasible
// 30 - Too big
// 31 - Relaxed pb not feasible
// 32 - Optim solution found
// 33 - Optim solution found
// 34 - Optim solution found
// 35 - Relaxed pb not feasible
// 36 - Relaxed pb not feasible
// 37 - Must add a delta to x0 - memory leak make scilab hangs after a while
// 38 - Must add a delta to x0 - not feasible after a while
// 39 - Relaxed pb not feasible
// 40 - Must add a delta to x0 - not feasible after a while
// 41 - Must add a delta to x0 - not feasible after a while
// 42 - Must add a delta to x0 - not feasible after a while
nl_index = 24;// 18 // Bug avec bonmin_optim.sce en premier + pb 18 de macminlp; 13, 14 17 (OK with ipopt NOK with bonmin): pb 
              // CoinOR + Index 13  = hs100 
              // ModNL: 52 - 62 - 462 - 489 ??
              // macminlp: 8 - 24 - 26 - [37 - 43]??

max_iter   = 100;
time_limit = 3600;

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

/////////////////////////

tmp_var_type = ampl_get_type(asl);
printf('variables type: %s\n', tmp_var_type);

var_type = zeros(1,length(tmp_var_type));

var_type(strindex(tmp_var_type,'b')) = 1;
var_type(strindex(tmp_var_type,'i')) = 2; 

INT_MAX = 100000;

/////////////////////////

// Define type of variables and constraints
nb_constr = length(constr_rhs);
nb_var    = length(x0);

var_lin_type    = ones(nb_var,1);    // 0 Linear - 1 Non-Linear
constr_lin_type = ones(nb_constr,1); // 0 Linear - 1 Non-Linear

printf('\n');
printf('number of variables:         %d\n', nb_var);
printf('number of integer variables: %d\n', length(strindex(tmp_var_type,'i')));
printf('number of binary variables:  %d\n', length(strindex(tmp_var_type,'b')));
printf('number of constraints:       %d\n', nb_constr);

////////////////////////////////////////////////////////////////////////

params = init_param();

///////////////////
// Set algorithm //
///////////////////

// Possible choices:
// B-BB  simple branch-and-bound algorithm
// B-OA  OA Decomposition algorithm
// B-QG  Quesada and Grossmann branch-and-cut algorithm
// B-Hyb hybrid outer approximation based branch-and-cut
// B-Ecp ecp cuts based branch-and-cut a la FilMINT
//params = add_param(params, 'algorithm', 'B-BB');
//params = add_param(params, 'algorithm', 'B-OA');
params = add_param(params, 'algorithm', 'B-QG');
//params = add_param(params, 'algorithm', 'B-Hyb');
//params = add_param(params, 'algorithm', 'B-Ecp');

////////////////////////////////////////
//  Bonmin ecp based strong branching //
////////////////////////////////////////

// Set the relative termination tolerance for ECP rounds in strong branching.
params = add_param(params, 'ecp_abs_tol_strong', 1e-6);

// Set the absolute termination tolerance for ECP rounds in strong branching.
params = add_param(params, 'ecp_max_rounds_strong', 0);

params = add_param(params, 'ecp_rel_tol_strong', 0.1);

// Set the relative termination tolerance for ECP rounds in strong branching.
// Choose method to use for warm starting lp in strong branching:
// - Basis: Use optimal basis of node
// - Clone: Clone optimal problem of node (Advanced stuff)
params = add_param(params, 'lp_strong_warmstart_method', 'Basis');

///////////////////////////////
//  Branch-and-bound options //
///////////////////////////////

// Specify the value of relative gap under which the algorithm stops.
// Stop the tree search when the gap between the objective value of the best known solution and the best bound on the objective of any solution is less than this
// fraction of the absolute value of the best known solution value.
params = add_param(params, 'allowable_fraction_gap', 0);

// Specify the value of absolute gap under which the algorithm stops.
// Stop the tree search when the gap between the objective value of the best known solution and the best bound on the objective of any solution is less than this.
params = add_param(params, 'allowable_gap', 0);

// Specify cutoff value.
// cutoff should be the value of a feasible solution known by the user (if any). The algorithm will only look for solutions better than cutoof.
params = add_param(params, 'cutoff', 1e100);

// Specify cutoff decrement.
// Specify the amount by which cutoff is decremented below a new best upper-bound (usually a small positive value but in non-convex problems it may be a negative value).
params = add_param(params, 'cutoff_decr', 1e-5);

// Set integer tolerance
params = add_param(params, 'integer_tolerance', 1e-6);

// Set the cumulated maximum number of iteration in the algorithm used to process nodes continuous relaxations in the branch-and-bound
params = add_param(params, 'iteration_limit', INT_MAX);
//params = add_param(params, 'iteration_limit', 1000); // YC: verifier cette options

// Set the behavior when an NLP or a series of NLP are unsolved by Ipopt (we call unsolved an NLP for which Ipopt is not able to guarantee
// optimality within the specified tolerances).
// - stop: Stop when failure happens.
// - fathom: Continue when failure happens.
// If set to 'fathom', the algorithm will fathom the node when Ipopt fails to find a solution to the nlp at that node whithin the specified tolerances.
// The algorithm then becomes a heuristic, and the user will be warned that the solution might not be optimal.
params = add_param(params, 'nlp_failure_behavior', 'Stop');

// Choose the node selection strategy.
// - best-bound: choose node with the smallest bound
// - depth-first: Perform depth first search
// - breadth-first: Perform breadth first search
// - dynamic: Cbc dynamic strategy (starts with a depth first search and turn to best bound after 3 integer feasible solutions have been found)
// - best-guess: choose node with smallest guessed integer solution
//params = add_param(params, 'node_comparison', 'dynamic');
params = add_param(params, 'node_comparison', 'best-bound');

// Set the maximum number of nodes explored in the branch-and-bound search
params = add_param(params, 'node_limit', INT_MAX);

// Set the maximum number of cut passes at regular nodes of the branch-and-cut.
params = add_param(params, 'num_cut_passes', 1);

// Set the maximum number of cut passes at regular nodes of the branch-and-cut.
params = add_param(params, 'num_cut_passes_at_root', 20);

// Set the number of branches on a variable before its pseudo costs are to be believed in dynamic strong branching.
params = add_param(params, 'number_before_trust', 8);
params = add_param(params, 'number_strong_branch', 20);

// Abort after that much integer feasible solution have been found by algorithm
params = add_param(params, 'solution_limit', INT_MAX);

// Wether or not to activate SOS constraints.
// - enable
// - disable
// BE CAREFUL: maybe this option can make bonmin hangs on certain problems. Try not tu use it ...
//params = add_param(params, 'sos_constraints', 'enable');
params = add_param(params, 'sos_constraints', 'disable');

// Set the global maximum computation time (in secs) for the algorithm.
//params = add_param(params, 'time_limit', 1e10);
params = add_param(params, 'time_limit', time_limit);

// Pick a strategy for traversing the tree
// - probed-dive:
// - top-node: Always pick the top node as sorted by the node comparison function
// - dive: Dive in the tree if possible, otherwise pick top node as sorted by the tree comparison function.
// - probed-dive: Dive in the tree exploring two childs before continuing the dive at each level.
// - dfs-dive: Dive in the tree if possible doing a depth first search.
// Backtrack on leaves or when a prescribed depth is attained or when estimate of best possible integer feasible solution in subtree
// is worst than cutoff. Once a prescribed limit of backtracks is attained pick top node as sorted by the tree comparison function
// - dfs-dive-dynamic: Same as dfs-dive but once enough solution are found switch to best-bound and if too many nodes switch to depth-first.
// All strategies can be used in conjunction with any of the node comparison functions. Options which affect dfs-dive are max-backtracks-in-dive and max-dive-depth.
// The dfs-dive won't work in a non-convex problem where objective does not decrease down branches.
params = add_param(params, 'tree_search_strategy', 'top-node');

// Chooses variable selection strategy
// - most-fractional: Choose most fractional variable
// - strong-branching: Perform strong branching
// - reliability-branching: Use reliability branching
// - curvature-estimator: Use curvature estimation to select branching variable
// - qp-strong-branching: Perform strong branching with QP approximation
// - lp-strong-branching: Perform strong branching with LP approximation
// - nlp-strong-branching: Perform strong branching with NLP approximation
// - osi-simple: Osi method to do simple branching
// - osi-strong: Osi method to do strong branching
// - random: Method to choose branching variable randomly
params = add_param(params, 'variable_selection', 'strong-branching');
//params = add_param(params, 'variable_selection', 'most-fractional');

/////////////////////
//  Diving options //
/////////////////////

// Set the number of backtracks in a dive when using dfs-dive tree search strategy.
params = add_param(params, 'max_backtracks_in_dive', 5);

// When using dfs-dive search. Maximum depth to go to from the diving board (node where the diving started.
params = add_param(params, 'max_dive_depth', INT_MAX);

// Flag indicating whether we stop diving based on guessed feasible objective and the current cutoff
// - no
// - yes
params = add_param(params, 'stop_diving_on_cutoff', 'no');

/////////////////////////////////
//  Feasibility pump heuristic //
/////////////////////////////////

//params = add_param(params, 'feasibility_pump_objective_norm', 1); // NOT AVAILABLE
//params = add_param(params, 'heuristic_feasibility_pump', 'no'); // NOT AVAILABLE

//////////////////////////////////////
//  Fractional diving MIP heuristic //
//////////////////////////////////////

//params = add_param(params, 'heuristic_dive_MIP_fractional', 'no'); // NOT AVAILABLE

//////////////////////////////////
//  Fractional diving heuristic //
//////////////////////////////////

//params = add_param(params, 'heuristic_dive_fractional', 'no'); // NOT AVAILABLE

////////////////////////////////////
//  Local search based heuristics //
////////////////////////////////////

//params = add_param(params, 'dummy_pump_heuristic', 'no'); // NOT AVAILABLE
//params = add_param(params, 'fix_and_solve_heuristic', 'no'); // NOT AVAILABLE
//params = add_param(params, 'heuristic_RINS', 'no'); // NOT AVAILABLE
//params = add_param(params, 'heuristic_local_branching', 'no'); // NOT AVAILABLE
//params = add_param(params, 'local_search_node_limit', 1000); // NOT AVAILABLE
//params = add_param(params, 'local_search_solution_limit', 5); // NOT AVAILABLE
//params = add_param(params, 'local_search_time_limit', 60); // NOT AVAILABLE

////////////////////////////////////
//  MILP cutting planes in hybrid //
////////////////////////////////////

// Frequency k (in terms of nodes) for generating Gomory cuts in branch-and-cut.
// If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but 
// Cbc may decide to stop generating cuts, if not enough are generated at the root node
// if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts
// BE CAREFUL: you must add at least one cut (one of the priority must be less than 0 strictly)
params = add_param(params, '2mir_cuts', -5); // 0
params = add_param(params, 'Gomory_cuts', -5); // 0
params = add_param(params, 'clique_cuts', -5); // -5
params = add_param(params, 'cover_cuts', -5);
params = add_param(params, 'flow_covers_cuts', -5);
params = add_param(params, 'lift_and_project_cuts', 0);
params = add_param(params, 'mir_cuts', -5);
//params = add_param(params, 'probing_cuts', -5); // NOT AVAILABLE
params = add_param(params, 'reduce_and_split_cuts', 0);

//////////////////////////////
//  Nlp solution robustness //
//////////////////////////////

params = add_param(params, 'max_consecutive_failures', 10);
params = add_param(params, 'max_random_point_radius', 100000);
params = add_param(params, 'num_iterations_suspect', -1);
params = add_param(params, 'num_retry_unsolved_random_point', 0);
params = add_param(params, 'random_point_perturbation_interval', 1);
params = add_param(params, 'random_point_type', 'Jon');

/////////////////////////////////
//  Nlp solve options in B-Hyb //
/////////////////////////////////

params = add_param(params, 'nlp_solve_frequency', 10);
params = add_param(params, 'nlp_solve_max_depth', 10);
params = add_param(params, 'nlp_solves_per_depth', 1e30);

/////////////////////////////////////////////////////
//  Options for MILP subsolver in OA decomposition //
/////////////////////////////////////////////////////

params = add_param(params, 'milp_log_level', 0);
params = add_param(params, 'milp_subsolver', 'Cbc_D');

///////////////////////////////////
//  Options for OA decomposition //
///////////////////////////////////

params = add_param(params, 'oa_dec_time_limit', 30);
params = add_param(params, 'oa_log_frequency', 100);
params = add_param(params, 'oa_log_level', 1);

//////////////////////////////////////
//  Options for ecp cuts generation //
//////////////////////////////////////

params = add_param(params, 'ecp_abs_tol', 1e-6);
params = add_param(params, 'ecp_max_rounds', 5);
params = add_param(params, 'ecp_propability_factor', 1000);
params = add_param(params, 'ecp_rel_tol', 0);
params = add_param(params, 'filmint_ecp_cuts', 0);

//////////////////////////////////////
//  Options for non-convex problems //
//////////////////////////////////////

params = add_param(params, 'max_consecutive_infeasible', 3); // 0
params = add_param(params, 'num_resolve_at_infeasibles', 5); // 0
params = add_param(params, 'num_resolve_at_node', 3); // 0
params = add_param(params, 'num_resolve_at_root', 3); // 0

//////////////////////////////////////////
//  Outer Approximation cuts generation //
//////////////////////////////////////////

params = add_param(params, 'add_only_violated_oa', 'no');
params = add_param(params, 'cut_strengthening_type', 'none');
params = add_param(params, 'disjunctive_cut_type', 'none');
params = add_param(params, 'oa_cuts_log_level', 0);
params = add_param(params, 'oa_cuts_scope', 'global');
params = add_param(params, 'tiny_element', 1e-8);
params = add_param(params, 'very_tiny_element', 1e-17);

////////////////////////////////////
//  Output ond log-levels options //
////////////////////////////////////

// bb_log_interval description:
// Set the interval (in terms of number of nodes) at which a log on node resolutions (consisting of lower and upper bounds) is given.
params = add_param(params, 'bb_log_interval', 1);

// bb_log_level description:
// Set the level of output of the branch-and-bound :
// 0 - none, 1 - minimal, 2 - normal low, 3 - normal high
params = add_param(params, 'bb_log_level', 1);

// lp_log_level description:
// Set the level of output of the linear programming sub-solver in B-Hyb or B-QG :
// 0 - none, 1 - minimal, 2 - normal low, 3 - normal high, 4 - verbose
params = add_param(params, 'lp_log_level', 1); // 0

/////////////////////////////
//  Strong branching setup //
/////////////////////////////

params = add_param(params, 'candidate_sort_criterion', 'best-ps-cost');
params = add_param(params, 'maxmin_crit_have_sol', 0.1);
params = add_param(params, 'maxmin_crit_no_sol', 0.7);

// Choose the maximum number of variables considered for strong branching.
params = add_param(params, 'min_number_strong_branch', 0);
params = add_param(params, 'number_before_trust_list', 0);
params = add_param(params, 'number_look_ahead', 0);
params = add_param(params, 'number_strong_branch_root', INT_MAX);
params = add_param(params, 'setup_pseudo_frac', 0.5);
params = add_param(params, 'trust_strong_branching_for_pseudo_cost', 'yes');

////////////////////////////////////////
//  VectorLength diving MIP heuristic //
////////////////////////////////////////

//params = add_param(params, 'heuristic_dive_MIP_vectorLength', 'no'); // NOT AVAILABLE

////////////////////////////////////
//  VectorLength diving heuristic //
////////////////////////////////////

//params = add_param(params, 'heuristic_dive_vectorLength', 'no'); // NOT AVAILABLE

///////////////////////////
//  nlp interface option //
///////////////////////////

params = add_param(params, 'file_solution', 'no');
params = add_param(params, 'nlp_log_level', 2);
params = add_param(params, 'nlp_solver', 'Ipopt');
params = add_param(params, 'warm_start', 'none');

// Linear solver used for step computations. 
// Determines which linear algebra package is to be used for the solution of the augmented linear system (for obtaining the search directions). 
// Note, the code must have been compiled with the linear solver you want to choose. Depending on your Ipopt installation, not all options are available. 
// The default value for this string option is 'ma27'.
// Possible values:
//  * ma27: use the Harwell routine MA27
//  * ma57: use the Harwell routine MA57
//  * pardiso: use the Pardiso package
//  * wsmp: use WSMP package
//  * mumps: use MUMPS package
//  * custom: use custom linear solver
params = add_param(params,'ipopt.linear_solver','ma27'); // YC trouver avec quel solveur on a compile - normalement, ma27
//params = add_param(params,'ipopt.linear_solver','mumps'); // YC: linux - pour l'instant 

/////////////////////////
// Hessian information //
/////////////////////////
// 0 - exact
// 1 - limited memory
params = add_param(params, 'hessian_approximation', 'limited-memory');

/////////////////////////////////////
// Messages for the scilab console //
/////////////////////////////////////

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

params = add_param(params, 'journal_level', 5); // 3 for a summary at the end of the run and 5 for a summary after each iterations

//////////////////////////////////////////////////////////////
//                      Ipopt options                       //
// copy paster from sciipopt + remove the duplicate options //
//////////////////////////////////////////////////////////////

// Maximum number of iterations
params = add_param(params,'max_iter', max_iter);

/////////////////
// Termination //
/////////////////
// 'Acceptance' stopping criterion based on objective function change.
// If the relative change of the objective function (scaled by Max(1,|f(x)|)) is less than this value, this part of the acceptable tolerance termination 
// is satisfied; see also acceptable_tol. This is useful for the quasi-Newton option, which has trouble to bring down the dual infeasibility.
params = add_param(params,'acceptable_obj_change_tol', 1e20);
// Number of 'acceptable' iterates before triggering termination.
// If the algorithm encounters this many successive 'acceptable' iterates (see 'acceptable_tol'), it terminates, assuming that the problem 
// has been solved to best possible accuracy given round-off. If it is set to zero, this heuristic is disabled.
params = add_param(params,'acceptable_iter', 15);

/////////
// NLP //
/////////
// Indicates whether it is desired to check for Nan/Inf in derivative matrices
// * no: Don't check (faster).
// * yes: Check Jacobians and Hessian for Nan and Inf.
// Activating this option will cause an error if an invalid number is detected in the constraint Jacobians or the Lagrangian Hessian. 
// If this is not activated, the test is skipped, and the algorithm might proceed with invalid numbers and fail.
params = add_param(params,'check_derivatives_for_naninf', 'no');
// any bound less or equal this value will be considered -inf (i.e. not lower bounded).
params = add_param(params,'nlp_lower_bound_inf', -1e19);
// any bound greater or this value will be considered +inf (i.e. not upper bounded).
params = add_param(params,'nlp_upper_bound_inf', 1e19);
// Determines how fixed variables should be handled.
// * make_parameter: Remove fixed variable from optimization variables
// * make_constraint: Add equality constraints fixing variables
// * relax_bounds: Relax fixing bound constraints
// The main difference between those options is that the starting point in the 'make_constraint' case still has the fixed variables at 
// their given values, whereas in the case 'make_parameter' the functions are always evaluated with the fixed values for those variables.  
// Also, for 'relax_bounds', the fixing bound constraints are relaxed (according to 'bound_relax_factor'). For both 'make_constraints'
//  and 'relax_bounds', bound multipliers are computed for the fixed variables.
params = add_param(params,'fixed_variable_treatment', 'make_parameter');
//Indicates which linear solver should be used to detect linearly dependent equality constraints.
// none: don't check; no extra work at beginning
// mumps: use MUMPS
// wsmp: use WSMP
// ma28: use MA28
// The default and available choices depend on how Ipopt has been compiled. This is experimental and does not work well.
params = add_param(params,'dependency_detector', 'none');
// Indicates if the right hand sides of the constraints should be considered during dependency detection
// no: only look at gradients
// yes: also consider right hand side
params = add_param(params,'dependency_detection_with_rhs', 'no');
// Number of linear variables
// When the Hessian is approximated, it is assumed that the first num_linear_variables variables are linear. The Hessian is then not 
// approximated in this space. If the get_number_of_nonlinear_variables method in the TNLP is implemented, this option is ignored.
params = add_param(params,'num_linear_variables', 0);
// Indicates whether all equality constraints are linear
// * no: Don't assume that all equality constraints are linear
// * yes: Assume that equality constraints Jacobian are constant
// Activating this option will cause Ipopt to ask for the Jacobian of the equality constraints only once from the NLP and reuse this information later.
params = add_param(params,'jac_c_constant', 'no');
// Indicates whether all inequality constraints are linear
// * no: Don't assume that all inequality constraints are linear
// * yes: Assume that equality constraints Jacobian are constant
// Activating this option will cause Ipopt to ask for the Jacobian of the inequality constraints only once from the NLP and reuse this information later.
params = add_param(params,'jac_d_constant', 'no');
// Indicates whether the problem is a quadratic problem
// * no: Assume that Hessian changes
// * yes: Assume that Hessian is constant
// Activating this option will cause Ipopt to ask for the Hessian of the Lagrangian function only once from the NLP and reuse this information later.
params = add_param(params,'hessian_constant', 'no');

////////////////////
// Initialization //
////////////////////
// Initialization method for bound multipliers
// * constant: set all bound multipliers to the value of bound_mult_init_val
// * mu-based: initialize to mu_init/x_slack
// This option defines how the iterates for the bound multipliers are initialized.  If 'constant' is chosen, then all bound multipliers 
// are initialized to the value of 'bound_mult_init_val'.  If 'mu-based' is chosen, the each value is initialized to the the value 
// of 'mu_init' divided by the corresponding slack variable. This latter option might be useful if the starting point is close to the
// optimal solution.
params = add_param(params,'bound_mult_init_method', 'constant');
// Least square initialization of the primal variables
// no: take user-provided point
// yes: overwrite user-provided point with least-square estimates
// If set to yes, Ipopt ignores the user provided point and solves a least square problem for the primal variables (x and s), to fit the 
// linearized equality and inequality constraints.  This might be useful if the user doesn't know anything about the starting point, or for 
// solving an LP or QP.
params = add_param(params,'least_square_init_primal', 'no');
// Least square initialization of all dual variables
// no: use bound_mult_init_val and least-square equality constraint multipliers
// yes: overwrite user-provided point with least-square estimates
// If set to yes, Ipopt tries to compute least-square multipliers (considering ALL dual variables). If successful, the bound 
// multipliers are possibly corrected to be at least bound_mult_init_val. This might be useful if the user doesn't know anything 
// about the starting point, or for solving an LP or QP. This overwrites option 'bound_mult_init_method'.
params = add_param(params,'least_square_init_duals', 'no');

///////////////////////
// Barrier Parameter //
///////////////////////
// Indicates if we want to do Mehrotra's algorithm.
// * no: Do the usual Ipopt algorithm.
// * yes: Do Mehrotra's predictor-corrector algorithm.
// If set to yes, Ipopt runs as Mehrotra's predictor-corrector algorithm. This works usually very well for LPs and convex QPs. This
// automatically disables the line search, and chooses the (unglobalized) adaptive mu strategy with the 'probing' oracle, and uses 
// 'corrector_type=affine' without any safeguards; you should not set any of those options explicitly in addition. Also, unless 
// otherwise specified, the values of 'bound_push', 'bound_frac', and 'bound_mult_init_val' are set more aggressive, and sets 
// 'alpha_for_y=bound_mult'.
params = add_param(params,'mehrotra_algorithm', 'no');
params = add_param(params,'quality_function_max_section_steps', 8); // YC
// Oracle for the barrier parameter when switching to fixed mode.
// * probing: Mehrotra's probing heuristic
// * loqo: LOQO's centrality rule
// * quality-function: minimize a quality function
// * average_compl: base on current average complementarity
// Determines how the first value of the barrier parameter should be computed when switching to the 'monotone mode' in the adaptive strategy. 
// (Only considered if 'adaptive' is selected for option 'mu_strategy')
params = add_param(params,'fixed_mu_oracle', 'average_compl');
// Factor for initialization of maximum value for barrier parameter.
// This option determines the upper bound on the barrier parameter. This upper bound is computed as the average complementarity at the initial 
// point times the value of this option. (Only used if option 'mu_strategy' is chosen as 'adaptive'.)
params = add_param(params,'mu_max_fact', 1000);
// Minimum value for barrier parameter.
// This option specifies the lower bound on the barrier parameter in the adaptive mu selection mode. By default, it is set to the minimum of 1e-11 and 
// min('tol','compl_inf_tol')/('barrier_tol_factor'+1), which should be a reasonable value. (Only used if option 'mu_strategy' is chosen as 'adaptive'.)
params = add_param(params,'mu_min', 1e-11);
// Factor for mu in barrier stop test.
// The convergence tolerance for each barrier problem in the monotone mode is the value of the barrier parameter times 'barrier_tol_factor'.
// This option is also used in the adaptive mu strategy during the monotone mode. (This is kappa_epsilon in implementation paper).
params = add_param(params,'barrier_tol_factor', 10);
// Determines linear decrease rate of barrier parameter.
// For the Fiacco-McCormick update procedure the new barrier parameter mu is obtained by taking the minimum of mu*'mu_linear_decrease_factor'
// and mu^'superlinear_decrease_power'. (This is kappa_mu in implementation paper.) This option is also used in the adaptive mu 
// strategy during the monotone mode.
params = add_param(params,'mu_linear_decrease_factor', 0.2);
// Determines superlinear decrease rate of barrier parameter.
// For the Fiacco-McCormick update procedure the new barrier parameter mu is obtained by taking the minimum of mu*'mu_linear_decrease_factor'
// and mu^'superlinear_decrease_power'. (This is theta_mu in implementation paper.) This option is also used in the adaptive mu 
// strategy during the monotone mode.
params = add_param(params,'mu_superlinear_decrease_power', 1.5);
// Allow skipping of barrier problem if barrier test is already met.
// no: Take at least one iteration per barrier problem
// yes: Allow fast decrease of mu if barrier test it met
// If set to 'no', the algorithm enforces at least one iteration per barrier problem, even if the barrier test is already met for the updated barrier parameter.
params = add_param(params,'mu_allow_fast_monotone_decrease', 'yes');
// Lower bound on fraction-to-the-boundary parameter tau
// (This is tau_min in the implementation paper.)  This option is also used in the adaptive mu strategy during the monotone mode.
params = add_param(params,'tau_min', 0.99);
// Factor limiting the deviation of dual variables from primal estimates.
// If the dual variables deviate from their primal estimates, a correction is performed. (See Eqn. (16) in the implementation paper.) 
// Setting the value to less than 1 disables the correction.
params = add_param(params,'kappa_sigma', 1e10);
// Globalization method used in backtracking line search
// filter: Filter method
// cg-penalty: Chen-Goldfarb penalty function
// penalty: Standard penalty function
//params = add_param(params,'line_search_method', 'cg-penalty');
params = add_param(params,'line_search_method', 'penalty'); // YC: warning, this options must be chosen carefully when we use ipopt under bonmin
// Globalization strategy for the adaptive mu selection mode.
// kkt-error: nonmonotone decrease of kkt-error
// obj-constr-filter: 2-dim filter for objective and constraint violation
// never-monotone-mode: disables globalization
// To achieve global convergence of the adaptive version, the algorithm has to switch to the monotone mode (Fiacco-McCormick approach) when 
// convergence does not seem to appear.  This option sets the criterion used to decide when to do this switch. (Only used if option 
// 'mu_strategy' is chosen as 'adaptive'.)
params = add_param(params,'adaptive_mu_globalization', 'obj-constr-filter');
// Maximum number of iterations requiring sufficient progress.
// For the 'kkt-error' based globalization strategy, sufficient progress must be made for 'adaptive_mu_kkterror_red_iters
// iterations. If this number of iterations is exceeded, the globalization strategy switches to the monotone mode.
params = add_param(params,'adaptive_mu_kkterror_red_iters', 4);
// Sufficient decrease factor for 'kkt-error' globalization strategy.
// For the 'kkt-error' based globalization strategy, the error must decrease by this factor to be deemed sufficient decrease.
params = add_param(params,'adaptive_mu_kkterror_red_fact', 0.9999);
// Factor determining width of margin for obj-constr-filter adaptive globalization strategy.
// When using the adaptive globalization strategy, 'obj-constr-filter' sufficient progress for a filter entry is defined as 
// follows: (new obj) < (filter obj) - filter_margin_fact*(new constr-viol) OR (new constr-viol) < (filter constr-viol) - 
// filter_margin_fact*(new constr-viol).  For the description of the 'kkt-error-filter' option see 'filter_max_margin.
params = add_param(params,'filter_margin_fact', 1e-5);
// Maximum width of margin in obj-constr-filter adaptive globalization strategy.
params = add_param(params,'filter_max_margin', 1.0);
// Indicates if the previous iterate should be restored if the monotone mode is entered.
// no: don't restore accepted iterate
// yes: restore accepted iterate
// When the globalization strategy for the adaptive barrier algorithm switches to the monotone mode, it can either start 
// from the most recent iterate (no), or from the last iterate that was accepted (yes).
params = add_param(params,'adaptive_mu_restore_previous_iterate', 'no');
// Determines the initial value of the barrier parameter when switching to the monotone mode.
// When the globalization strategy for the adaptive barrier algorithm switches to the monotone mode and fixed_mu_oracle is chosen as 
// 'average_compl', the barrier parameter is set to the current average complementarity times the value of 'adaptive_mu_monotone_init_factor'
params = add_param(params,'adaptive_mu_monotone_init_factor', 0.8);
// Norm used for the KKT error in the adaptive mu globalization strategies.
// 1-norm: use the 1-norm (abs sum)
// 2-norm-squared: use the 2-norm squared (sum of squares)
// max-norm: use the infinity norm (max)
// 2-norm: use 2-norm
// When computing the KKT error for the globalization strategies, the norm to be used is specified with this option. Note, this options is also used 
// in the QualityFunctionMuOracle.
params = add_param(params,'adaptive_mu_kkt_norm_type', '2-norm-squared');

////////////////////////
// Multiplier Updates //
////////////////////////
// Tolerance for switching to full equality multiplier steps
// This is only relevant if 'alpha_for_y' is chosen 'primal-and-full' or 'dual-and-full'. The step size for the equality constraint 
// multipliers is taken to be one if the max-norm of the primal step is less than this tolerance
params = add_param(params,'alpha_for_y_tol', 10.0);
// Tolerance for detecting numerically insignificant steps
// If the search direction in the primal variables (x and s) is, in relative terms for each component, less than this value, the 
// algorithm accepts the full step without line search. If this happens repeatedly, the algorithm will terminate with a corresponding exit 
// message. The default value is 10 times machine precision.  
params = add_param(params,'tiny_step_tol', 10.0*%eps);
// Tolerance for quitting because of numerically insignificant steps.
// If the search direction in the primal variables (x and s) is, in relative terms for each component, repeatedly less than tiny_step_tol,
// and the step in the y variables is smaller than this threshold, the algorithm will terminate.
params = add_param(params,'tiny_step_y_tol', 1e-2);
// Tells the algorithm to recalculate the equality and inequality multipliers as least square estimates.
// no: use the Newton step to update the multipliers
// yes: use least-square multiplier estimates
// This asks the algorithm to recompute the multipliers, whenever the current infeasibility is less than recalc_y_feas_tol. 
// Choosing yes might be helpful in the quasi-Newton option.  However, each recalculation requires an extra factorization of the linear 
// system.  If a limited memory quasi-Newton option is chosen, this is used by default.
params = add_param(params,'recalc_y', 'no');
// Feasibility threshold for recomputation of multipliers.
// If recalc_y is chosen and the current infeasibility is less than this value, then the multipliers are recomputed.
params = add_param(params,'recalc_y_feas_tol', 1e-6);

/////////////////
// Line Search //
/////////////////
// Maximum number of watchdog iterations.
// This option determines the number of trial iterations allowed before the watchdog procedure is aborted and the algorithm returns to the stored point.
params = add_param(params,'watchdog_trial_iter_max', 3);
// * none: no corrector
// * affine: corrector step towards mu=0
// * primal-dual: corrector step towards current mu
//params = add_param(params,'corrector_type', 'none');
params = add_param(params,'corrector_type', 'affine');

///////////////////
// Line Search 2 //
///////////////////
// Maximal number of second order correction trial steps
params = add_param(params,'max_soc', 4);
// Trigger counter for watchdog procedure
//params = add_param(params,'watchdog_shortened_iter_trigger', 10);
params = add_param(params,'watchdog_shortened_iter_trigger', 0); // YC: bonmin imposes this value to be equal to 0 !
// Update strategy for barrier parameter
// * monotone: use the monotone (Fiacco-McCormick) strategy
// * adaptive: use the adaptive update strategy
//params = add_param(params,'mu_strategy', 'monotone');
params = add_param(params,'mu_strategy', 'adaptive'); // BONMIN set this option
// Oracle for a new barrier parameter in the adaptive strategy
// * probing: Mehrotra's probing heuristic
// * loqo: LOQO's centrality rule
// * quality-function: minimize a quality function
//params = add_param(params,'mu_oracle', 'quality-function');
params = add_param(params,'mu_oracle', 'probing');
// Initial value for the barrier parameter
params = add_param(params,'mu_init', 0.1);
// Maximal value for barrier parameter for adaptive strategy
params = add_param(params,'mu_max', 100000);
// Select the technique used for scaling the NLP
//  none: no problem scaling will be performed
//  user-scaling: scaling parameters will come from the user
//  gradient-based: scale the problem so the maximum gradient at the starting point is scaling_max_gradient
//  equilibration-based: scale the problem so that first derivatives are of order 1 at random points (only available with MC19)
params = add_param(params,'nlp_scaling_method', 'gradient-based');
// Maximum gradient after scaling
params = add_param(params,'nlp_scaling_max_gradient', 100);
// Scaling factor for the objective function
params = add_param(params,'obj_scaling_factor', 1);

////////////////
// Warm Start //
////////////////
params = add_param(params,'warm_start_bound_frac', 1e-3);
params = add_param(params,'warm_start_slack_bound_frac', 1e-3);
params = add_param(params,'warm_start_slack_bound_push', 1e-3);
params = add_param(params,'warm_start_mult_init_max', 1e6);

///////////////////////
// Restoration Phase //
///////////////////////
// Threshold for disabling 'expect_infeasible_problem' option.
// If the constraint violation becomes smaller than this threshold, the 'expect_infeasible_problem' heuristics in the filter line 
// search are disabled. If the problem is square, this options is set to 0.
params = add_param(params,'expect_infeasible_problem_ctol', 1e-3);
// Tells algorithm to switch to restoration phase in first iteration.
// * no: don't force start in restoration phase
// * yes: force start in restoration phase
// Setting this option to 'yes' forces the algorithm to switch to the feasibility restoration phase in the first iteration. If the initial 
// point is feasible, the algorithm will abort with a failure.
params = add_param(params,'start_with_resto', 'no');
//params = add_param(params,'start_with_resto', 'yes');
// Required reduction in primal-dual error in the soft restoration phase.
// The soft restoration phase attempts to reduce the primal-dual error with regular steps. If the damped 
// primal-dual step (damped only to satisfy the fraction-to-the-boundary rule) is not decreasing the primal-dual error 
// by at least this factor, then the regular restoration phase is called. Choosing '0' here disables the soft restoration phase.
params = add_param(params,'soft_resto_pderror_reduction_factor', 0.9999);
// Maximum number of iterations performed successively in soft restoration phase.
// If the soft restoration phase is performed for more than so many iterations in a row, the regular restoration phase is called.
params = add_param(params,'max_soft_resto_iters', 10);
// Maximum number of successive iterations in restoration phase.
// The algorithm terminates with an error message if the number of iterations successively taken in the restoration phase exceeds this number.
params = add_param(params,'max_resto_iter', 3000000);
// Threshold for resetting bound multipliers after the restoration phase.
// After returning from the restoration phase, the bound multipliers are updated with a Newton step for complementarity. Here, the
// change in the primal variables during the entire restoration phase is taken to be the corresponding primal Newton step. 
// However, if after the update the largest bound multiplier exceeds the threshold specified by this option, the multipliers are all reset to 1.
params = add_param(params,'bound_mult_reset_threshold', 1000);
// Threshold for resetting equality and inequality multipliers after restoration phase.
// After returning from the restoration phase, the constraint multipliers are recomputed by a least square estimate.  This option triggers when
// those least-square estimates should be ignored.
params = add_param(params,'constr_mult_reset_threshold', 0);
// Determines if the original objective function should be evaluated at restoration phase trial points.
// Setting this option to 'yes' makes the restoration phase algorithm evaluate the objective function of the original problem at every trial 
// point encountered during the restoration phase, even if this value is not required. In this way, it is guaranteed that the original 
// objective function can be evaluated without error at all accepted iterates; otherwise the algorithm might fail at a point where the 
// restoration phase accepts an iterate that is good for the restoration phase problem, but not the original problem. On the other hand, if 
// the evaluation of the original objective is expensive, this might be costly.
// * no: skip evaluation
// * yes: evaluate at every trial point
//params = add_param(params,'evaluate_orig_obj_at_resto_trial', 'yes');
params = add_param(params,'evaluate_orig_obj_at_resto_trial', 'no');
// Penalty parameter in the restoration phase objective function.
// This is the parameter rho in equation (31a) in the Ipopt implementation paper.
params = add_param(params,'resto_penalty_parameter', 1000.0);

// Enable heuristics to quickly detect an infeasible problem
// * no: the problem probably be feasible
// * yes: the problem has a good chance to be infeasible
// BONMIN sets this option to true in BonIpoptSolver.cpp
params = add_param(params,'expect_infeasible_problem', 'yes');

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
params = add_param(params,'max_hessian_perturbation', 1e20);
// Smallest perturbation of the Hessian block.
// The size of the perturbation of the Hessian block is never selected smaller than this value, unless no perturbation is necessary. (This 
// is delta_w^min in implementation paper.)
params = add_param(params,'min_hessian_perturbation', 1e-20);
// Size of first x-s perturbation tried.
// The first value tried for the x-s perturbation in the inertia correction scheme.
// (This is delta_0 in the implementation paper.)
params = add_param(params,'first_hessian_perturbation', 1e-4);
// Increase factor for x-s perturbation for very first perturbation.
// The factor by which the perturbation is increased when a trial value was not sufficient - this value is used for the computation of the 
// very first perturbation and allows a different value for for the first perturbation than that used for the remaining perturbations. 
// (This is bar_kappa_w^+ in the implementation paper.)
params = add_param(params,'perturb_inc_fact_first', 100);
// Increase factor for x-s perturbation.
// The factor by which the perturbation is increased when a trial value was not sufficient - this value is used for the computation of 
// all perturbations except for the first. (This is kappa_w^+ in the implementation paper.)
params = add_param(params,'perturb_inc_fact', 8);
// Decrease factor for x-s perturbation.
// The factor by which the perturbation is decreased when a trial value is deduced from the size of the most recent successful perturbation. 
// (This is kappa_w^- in the implementation paper.)
params = add_param(params,'perturb_dec_fact', 0.333333);
// Size of the regularization for rank-deficient constraint Jacobians.
// (This is bar delta_c in the implementation paper.)
params = add_param(params,'jacobian_regularization_value', 1e-8);
// Exponent for mu in the regularization for rank-deficient constraint Jacobians.
// (This is kappa_c in the implementation paper.)
params = add_param(params,'jacobian_regularization_exponent', 0.25);
// Active permanent perturbation of constraint linearization.
// no: perturbation only used when required
// yes: always use perturbation
// This options makes the delta_c and delta_d perturbation be used for the computation of every search direction.  Usually, it is only used 
// when the iteration matrix is singular.
params = add_param(params,'perturb_always_cd', 'no'); 

//////////////////
// Quasi-Newton //
//////////////////
// Indicates in which subspace the Hessian information is to be approximated.
// nonlinear-variables: only in space of nonlinear variables.
// all-variables: in space of all variables (without slacks)
params = add_param(params,'hessian_approximation_space', 'nonlinear-variables');
// Maximum size of the history for the limited quasi-Newton Hessian approximation.
// This option determines the number of most recent iterations that are taken into account for the limited-memory quasi-Newton approximation.
params = add_param(params,'limited_memory_max_history', 6);
// Threshold for successive iterations where update is skipped.
// If the update is skipped more than this number of successive iterations, we quasi-Newton approximation is reset.
params = add_param(params,'limited_memory_max_skipping', 2);
// Quasi-Newton update formula for the limited memory approximation.
// bfgs: BFGS update (with skipping)
// sr1:  SR1 (not working well)
// Determines which update formula is to be used for the limited-memory quasi-Newton approximation.
params = add_param(params,'limited_memory_update_type', 'bfgs');
// Initialization strategy for the limited memory quasi-Newton approximation.
// scalar1: sigma = s^Ty/s^Ts
// scalar2: sigma = y^Ty/s^Ty
// constant: sigma = limited_memory_init_val
// Determines how the diagonal Matrix B_0 as the first term in the limited memory approximation should be computed.
params = add_param(params,'limited_memory_initialization', 'scalar1');
// Upper bound on value for B0 in low-rank update.
// The starting matrix in the low rank update, B0, is chosen to be this multiple of the identity in the first iteration (when no updates have 
// been performed yet), and is constantly chosen as this value, if 'limited_memory_initialization' is 'constant'.
params = add_param(params,'limited_memory_init_val_max', 1e8);
// Lower bound on value for B0 in low-rank update.
// The starting matrix in the low rank update, B0, is chosen to be this multiple of the identity in the first iteration (when no updates have 
// been performed yet), and is constantly chosen as this value, if 'limited_memory_initialization' is 'constant'.
params = add_param(params,'limited_memory_init_val_min', 1e-8);
// Value for B0 in low-rank update
// The starting matrix in the low rank update, B0, is chosen to be this multiple of the identity in the first iteration (when no updates have
// 'been performed yet), and is constantly chosen as this value, if 'limited_memory_initialization' is 'constant'.
params = add_param(params,'limited_memory_init_val', 1);

/////////////////////
// Derivative Test //
/////////////////////
// Enable derivative checker
// * none: do not perform derivative test
// * first-order: perform test of first derivatives at starting point
// * second-order: perform test of first and second derivatives at starting point
// If this option is enabled, a (slow) derivative test will be performed before the optimization. The test is performed at the user provided 
// starting point and marks derivative values that seem suspicious
params = add_param(params,'derivative_test', 'none');
// Size of the finite difference perturbation in derivative test.
// This determines the relative perturbation of the variable entries.
params = add_param(params,'derivative_test_perturbation', 1e-8);
// Threshold for indicating wrong derivative.
// If the relative deviation of the estimated derivative from the given one is larger than this value, the corresponding derivative is marked as wrong.
params = add_param(params,'derivative_test_tol', 1e-4);
// Indicates whether information for all estimated derivatives should be printed.
// * no: Print only suspect derivatives
// * yes: Print all derivatives
params = add_param(params,'derivative_test_print_all', 'no');
// Maximal perturbation of an evaluation point.
// If a random perturbation of a points is required, this number indicates the maximal perturbation. This is for example used when 
// determining the center point at which the finite difference derivative test is executed.
params = add_param(params,'point_perturbation_radius', 10);
// Specifies technique to compute constraint Jacobian
// exact: user-provided derivatives
// finite-difference-values: user-provided structure, values by finite differences
//params = add_param(params,'jacobian_approximation', 'exact');
params = add_param(params,'jacobian_approximation', 'finite-difference-values');
// Size of the finite difference perturbation for derivative approximation.
// This determines the relative perturbation of the variable entries.
params = add_param(params,'findiff_perturbation', 1e-7);

////////////////////////
// MA27 Linear Solver //
////////////////////////
// Integer workspace memory for MA27.
// The initial integer workspace memory = liw_init_factor * memory required by unfactored system. Ipopt will increase the workspace 
// size by meminc_factor if required.  This option is only available if Ipopt has been compiled with MA27.
params = add_param(params,'ma27_liw_init_factor', 5);
// Real workspace memory for MA27.
// The initial real workspace memory = la_init_factor * memory required by unfactored system. Ipopt will increase the workspace
// size by meminc_factor if required.  This option is only available if Ipopt has been compiled with MA27.
params = add_param(params,'ma27_la_init_factor', 5);
// Increment factor for workspace size for MA27.
// If the integer or real workspace is not large enough, Ipopt will increase its size by this factor.  This option is only 
// available if Ipopt has been compiled with MA27.
params = add_param(params,'ma27_meminc_factor', 5);
// Always pretend inertia is correct.
// no: check inertia
// yes: skip inertia check
// Setting this option to 'yes' essentially disables inertia check. This option makes the algorithm non-robust and easily fail, but it 
// might give some insight into the necessity of inertia control.
params = add_param(params,'ma27_skip_inertia_check', 'no');
// Enables MA27's ability to solve a linear system even if the matrix is singular.
// no: Don't have MA27 solve singular systems
// yes: Have MA27 solve singular systems
// Setting this option to \'yes\' means that Ipopt will call MA27 to compute solutions for right hand sides, even if MA27 has detected that 
// the matrix is singular (but is still able to solve the linear system). In some cases this might be better than using Ipopt's heuristic of 
// small perturbation of the lower diagonal of the KKT matrix.
params = add_param(params,'ma27_ignore_singularity', 'no');

////////////////////////
// MA57 Linear Solver //
////////////////////////
// Safety factor for work space memory allocation for the linear solver MA57.
// If 1 is chosen, the suggested amount of work space is used.  However, choosing a larger number might avoid reallocation if the suggest values 
// do not suffice.  This option is only available if Ipopt has been compiled with MA57.
params = add_param(params,'ma57_pre_alloc', 3);

/////////////////////////
// MUMPS Linear Solver //
/////////////////////////
// Pivot tolerance for the linear solver MUMPS.
// A smaller number pivots for sparsity, a larger number pivots for stability.  This option is only available if Ipopt has been compiled with MUMPS.
params = add_param(params,'mumps_pivtol', 1e-6);
// Maximum pivot tolerance for the linear solver MUMPS.
// Ipopt may increase pivtol as high as pivtolmax to get a more accurate solution to the linear system.  This option is only available if 
// Ipopt has been compiled with MUMPS.
params = add_param(params,'mumps_pivtolmax', 0.1);
// Percentage increase in the estimated working space for MUMPS.
// In MUMPS when significant extra fill-in is caused by numerical pivoting, larger values of mumps_mem_percent may help use the workspace more efficiently.
params = add_param(params,'mumps_mem_percent', 1000);
// Controls permuting and scaling in MUMPS
// This is ICTL(6) in MUMPS.
params = add_param(params,'mumps_permuting_scaling', 7);
// Controls pivot order in MUMPS
// This is ICTL(8) in MUMPS.
params = add_param(params,'mumps_pivot_order', 7);
// Controls scaling in MUMPS
// This is ICTL(8) in MUMPS.
params = add_param(params,'mumps_scaling', 7);
// Pivot threshold for detection of linearly dependent constraints in MUMPS.
// When MUMPS is used to determine linearly dependent constraints, this is determines the threshold for a pivot to be considered zero. This is CNTL(3) in MUMPS.
params = add_param(params,'mumps_dep_tol', -1.0);

///////////////////////////
// Pardiso Linear Solver //
///////////////////////////
// Toggle for handling case when elements were perturbed by Pardiso.
// no: Always redo symbolic factorization when elements were perturbed
// yes: Only redo symbolic factorization when elements were perturbed if also the inertia was wrong
// This option is only available if Ipopt has been compiled with Pardiso.
params = add_param(params,'pardiso_redo_symbolic_fact_only_if_inertia_wrong', 'no');
// Interpretation of perturbed elements
// no: Don't assume that matrix is singular if elements were perturbed after recent symbolic factorization
// yes: Assume that matrix is singular if elements were perturbed after recent symbolic factorization
// This option is only available if Ipopt has been compiled with Pardiso.
params = add_param(params,'pardiso_repeated_perturbation_means_singular', 'no');
// Always pretent inertia is correct.
// no: check inertia
// yes: skip inertia check
// Setting this option to 'yes' essentially disables inertia check. This option makes the algorithm non-robust and easily fail, but it 
// might give some insight into the necessity of inertia control.
params = add_param(params,'pardiso_skip_inertia_check', 'no');

////////////////////////
// WSMP Linear Solver //
////////////////////////
// Number of threads to be used in WSMP
// This determines on how many processors WSMP is running on. This option is only available if Ipopt has been compiled with WSMP.
params = add_param(params,'wsmp_num_threads', 1);
// Determines how ordering is done in WSMP
// This corresponds to the value of WSSMP's IPARM(16). This option is only available if Ipopt has been compiled with WSMP.
params = add_param(params,'wsmp_ordering_option', 1);
// Pivot tolerance for the linear solver WSMP.
// A smaller number pivots for sparsity, a larger number pivots for stability.  This option is only available if Ipopt has been compiled with WSMP.
params = add_param(params,'wsmp_pivtol', 1e-4);
// Maximum pivot tolerance for the linear solver WSMP.
// Ipopt may increase pivtol as high as pivtolmax to get a more accurate solution to the linear system. This option is only available if Ipopt 
// has been compiled with WSMP.
params = add_param(params,'wsmp_pivtolmax', 0.1);
// Determines how the matrix is scaled by WSMP.
// This corresponds to the value of WSSMP's IPARM(10). This option is only available if Ipopt has been compiled with WSMP.
params = add_param(params,'wsmp_scaling', 0);

///////////////////
// Launch bonmin //
///////////////////

tic();
[x_sol, f_sol, extra] = bonmin(x0, f, df, nb_constr, g, dg, sparse_dg, dh, sparse_dh, var_type, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);
t = toc();

// Status values:
// SUCCESS 	      0
// INFEASIBLE     1
// LIMIT_EXCEEDED 2
// MINLP_ERROR    3
printf('status           = %d\n',extra('status'));

// MIP status
// FeasibleOptimal  Optimum solution has been found and its optimality proved.
// ProvenInfeasible Problem has been proven to be infeasible.
// Feasible         An integer solution to the problem has been found.
// NoSolutionKnown  No feasible solution to the problem is known. 
printf('mip status       = %d\n',extra('mip_status'));
printf('best bound       = %f\n',extra('best_bound'));
printf('number of nodes  = %d\n',extra('num_nodes'));
printf('iterations count = %d\n',extra('iter_count'));

// bit 1 : Are there a numerical difficulties?
// bit 2 : Is optimality proven?
// bit 3 : Is primal infeasiblity proven?
// bit 4 : Is dual infeasiblity proven?
// bit 5 : Is the given primal objective limit reached?
// bit 6 : Is the given dual objective limit reached?
// bit 7 : Iteration limit reached?
printf('bonmin status    = %s\n',dec2bin(extra('bonmin_status')));
printf('time needed for resolution: %f\n', t);

if (extra('status')==0) then
  printf('\nconstraints result:\n');
  printf('lhs <= g(x_sol) <= rhs\n');

  disp(string(constr_lhs) + ' <= ' + string(g(x_sol)) + ' <= ' + string(constr_rhs));

  printf('\nSolution result:\n');
  printf('x_sol - var_type\n');
  disp(string(x_sol) + ' - ' + string(var_type'));
end

// Release the memory allocated for the selected AMPL problem
ampl_free(asl);
asl = [];

