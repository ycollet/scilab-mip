lines(0);

f  = [];
df = [];
g  = [];
dg = [];
h  = [];
dh = [];

sparse_dg = [];
sparse_dh = [];

current_dir = pwd();

///////////////////////////////////////////////////////////////////////////////

// int objective(double * fobj, double * x, int n_size_x, double x_new);

f_C = ['#include <math.h>'
       'int f_C(double * x, double * f, int n_size_x, double x_new)'
       '{'
       '  f[0] = 5*pow(x[0],2) - 3*pow(x[1],2);'
       '  return 1;'
       '}'];

// int objective_grad(double * fobj, double * x, int n_size_x, double x_new);

df_C = ['#include <math.h>'
        'int df_C(double * x, double * f, int n_size_x, double x_new)'
        '{'
        '  f[0] = 10*x[0];'
        '  f[1] = -6*x[1];'
        '  return 1;'
        '}'];

// int constraints(double * gobj, int n_size_constr, double * x, int n_size_x, double x_new);

g_C = ['#include <math.h>'
       'int g_C(double * x, int n_size_x, double * g, int n_size_g, double x_new)'
       '{'
       '  g[0] = -x[0];'
       '  g[1] = -x[1];'
       '  return 1;'
       '}'];

// int constraint_jac(double * x, int n_size_x, double new_x, int n_size_g, int nnz_jac, 
//                    int * iRow, int * jCol, double * values, void * param);

dg_C = ['#include <math.h>'
        '#include <stdio.h>'
        'int dg_C(double * x, int n_size_x, double new_x, int n_size_g, int nnz_jac, int * iRow, int * jCol, double * values)'
        '{'
        '  if ((iRow!=NULL)&&(jCol!=NULL))'
        '  {'
        '    iRow[0] = 1;'
        '    jCol[0] = 1;'
        '    iRow[1] = 2;'
        '    jCol[1] = 2;'
        '  }'
        '  if (values!=NULL)'
        '  {'
        '     values[0] = -1;'
        '     values[1] = -1;'
        '  }'
        '  return 1;'
        '}'];

cd TMPDIR;
mputl(f_C, TMPDIR+'/ipopt_demo_f_C.c');
mputl(df_C,TMPDIR+'/ipopt_demo_df_C.c');
mputl(g_C, TMPDIR+'/ipopt_demo_g_C.c');
mputl(dg_C,TMPDIR+'/ipopt_demo_dg_C.c');

// Get the current list of linked options

current_linked_list = link('show');

// compile the C code
printf('Compilation of the f_C function\n');
f_C_ilib_handle  = ilib_for_link('f_C', 'ipopt_demo_f_C.c', [],'c');
printf('Compilation of the df_C function\n');
df_C_ilib_handle = ilib_for_link('df_C','ipopt_demo_df_C.c',[],'c');
printf('Compilation of the g_C function\n');
g_C_ilib_handle  = ilib_for_link('g_C', 'ipopt_demo_g_C.c', [],'c');
printf('Compilation of the dg_C function\n');
dg_C_ilib_handle = ilib_for_link('dg_C','ipopt_demo_dg_C.c',[],'c');

// incremental linking
f_C_link_handle  = link(f_C_ilib_handle, 'f_C', 'c');
df_C_link_handle = link(df_C_ilib_handle,'df_C','c');
g_C_link_handle  = link(g_C_ilib_handle, 'g_C', 'c');
dg_C_link_handle = link(dg_C_ilib_handle,'dg_C','c');

new_linked_list = link('show');

[tmp,link_ids] = intersect(current_linked_list, new_linked_list);
new_linked_list(link_ids) = [];

cd(current_dir);

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

// Because we use a C function for the Jacobian, we need to set this options
// to specify the number of non zeros elements in the Jacobian
params = add_param(params,"nnz_jac",2); // MANDATORY PARAMETER HERE

///////////////////
// Launch bonmin //
///////////////////

[x_sol, f_sol, extra] = bonmin(x0, 'f_C', 'df_C', 'g_C', 'dg_C', sparse_dg, dh, sparse_dh, ...
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

// unlink the optimization problem

ulink([f_C_link_handle df_C_link_handle g_C_link_handle dg_C_link_handle]);
