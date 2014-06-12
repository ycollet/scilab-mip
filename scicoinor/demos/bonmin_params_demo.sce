lines(0);

f  = [];
df = [];
g  = [];
dg = [];
h  = [];
dh = [];

sparse_dg = [];
sparse_dh = [];

///////////////////////////////////////////////////////////////////////////////
// a = 5, b = 3
function y = f(x,x_new,a,b)
  y = a*x(1)^2-b*x(2)^2;
endfunction

function y = df(x,x_new,a,b)
  y(1) = 2*a*x(1);
  y(2) = -2*b*x(2);
endfunction

// a = -1, b = -1
function y = g(x,x_new,a,b)
  y(1) = a*x(1);
  y(2) = b*x(2);
endfunction

function y = dg(x,x_new,a,b)
  y(1) = a;
  y(2) = b;
endfunction

sparse_dg = [1 1; ...
             2 2];

upper = [4;4];
lower = [-4;-4];
x0 = [1;1];

var_type        = [];
var_lin_type    = [];
constr_lin_type = [];
constr_rhs      = [];
constr_lsh      = [];

var_type(1) = 2; // Integer
var_type(2) = 0; // Continuous
var_lin_type(1) = 0; // Linear
var_lin_type(2) = 0; // Linear
constr_lin_type (1) = 0; // Linear
constr_lin_type (2) = 0; // Linear
constr_rhs(1) = 0;
constr_rhs(2) = 0;
constr_lhs(1) = -10000;
constr_lhs(2) = -10000;

///////////////////////////////////////////////////////////////////////////////

INT_MAX = 100000;

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
params = add_param(params, 'algorithm', 'B-QG');

// Set the cumulated maximum number of iteration in the algorithm used to process nodes continuous relaxations in the branch-and-bound
params = add_param(params, 'iteration_limit', 10000);

params = add_param(params, 'nlp_failure_behavior', 'Stop');

// Set the maximum number of nodes explored in the branch-and-bound search
params = add_param(params, 'node_limit', 10000);

////////////////////////////////////
//  MILP cutting planes in hybrid //
////////////////////////////////////

// Frequency k (in terms of nodes) for generating Gomory cuts in branch-and-cut.
// If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but 
// Cbc may decide to stop generating cuts, if not enough are generated at the root node
// if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts
params = add_param(params, '2mir_cuts', -5);
params = add_param(params, 'Gomory_cuts', -5);
params = add_param(params, 'clique_cuts', -5);

//////////////////////////////
//  Nlp solution robustness //
//////////////////////////////

params = add_param(params, 'max_consecutive_failures', 10);
params = add_param(params, 'max_random_point_radius', 100000);
params = add_param(params, 'num_iterations_suspect', -1);
params = add_param(params, 'num_retry_unsolved_random_point', 0);
params = add_param(params, 'random_point_perturbation_interval', 1);
params = add_param(params, 'random_point_type', 'Jon');

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

params = add_param(params, 'journal_level', 5); // 3 for a summary at the end of the run and 5 for a summary after each iterations

// Maximum number of iterations
params = add_param(params,'max_iter', 100);

///////////////////
// Launch bonmin //
///////////////////

[x_sol, f_sol, extra] = bonmin(x0, list(f,5,3), list(df,5,3), list(g,-1,-1), list(dg,-1,-1), sparse_dg, dh, sparse_dh, ...
                               var_type, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

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

