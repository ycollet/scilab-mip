lines(0);

////////////////////////////////////////////////////////////////////////

f  = [];
df = [];
g  = [];
dg = [];
dh = [];
sparse_dg = [];
sparse_dh = [];

current_dir = pwd();

// int sci_ipopt_objective(double * x, double * f, int n_size_x, double x_new);

f_C = ['#include <math.h>'
       'int f_C(double * x, double * f, int n_size_x, double x_new)'
       '{'
       '  f[0] = 5*pow(x[0],2) - 3*pow(x[1],2);'
       '  return 1;'
       '}'];
       
// int sci_ipopt_objective_grad(double * x, double * f, int n_size_x, double x_new);

df_C = ['#include <math.h>'
        'int df_C(double * x, double * f, int n_size_x, double x_new)'
        '{'
        '  f[0] = 10*x[0];'
        '  f[1] = -6*x[1];'
        '  return 1;'
        '}'];
        
// int sci_ipopt_constraints(double * x, int n_size_x, double * g, int n_size_g, double x_new);

g_C = ['#include <math.h>'
       'int g_C(double * x, int n_size_x, double * g, int n_size_g, double x_new)'
       '{'
       '  g[0] = -x[0];'
       '  g[1] = -x[1];'
       '  return 1;'
       '}'];
        
// int sci_ipopt_constraints_jac(double * x, int n_size_x, double new_x, int n_size_g, int nnz_jac, int * iRow, int * jCol, double * values);

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

cd(current_dir);

upper = [4;4];
lower = [-4;-4];
x0    = [1;1];
  
var_lin_type(1) = 1; // Non-Linear
var_lin_type(2) = 1; // Non-Linear
constr_lin_type (1) = 0; // Linear
constr_lin_type (2) = 0; // Linear
constr_rhs(1) = 0;
constr_rhs(2) = 0;
constr_lhs(1) = -%inf;
constr_lhs(2) = -%inf;

////////////////////////////////////////////////////////////////////////

params = init_param();

/////////////
// Journal //
/////////////

// Journal verbosity level. 
// J_INSUPPRESSIBLE -1
// J_NONE 	     0
// J_ERROR 	     1
// J_STRONGWARNING   2
// J_SUMMARY 	     3
// J_WARNING 	     4
// J_ITERSUMMARY     5
// J_DETAILED        6
// J_MOREDETAILED    7
// J_VECTOR 	     8
// J_MOREVECTOR      9
// J_MATRIX 	     10
// J_MOREMATRIX      11
// J_ALL 	     12
// J_LAST_LEVEL      13
params = add_param(params,"journal_level",5);

params = add_param(params,"tol",1e-8);
params = add_param(params,"max_iter",3000);
// Because we use a C function for the Jacobian, we need to set this options
// to specify the number of non zeros elements in the Jacobian
params = add_param(params,"nnz_jac",2); // MANDATORY PARAMETER HERE

///////////////////
// Linear Solver //
///////////////////

params = add_param(params,"linear_solver","mumps"); 
params = add_param(params,"linear_system_scaling","none"); // YC: linux - pour l'instant

//////////////////
// Quasi-Newton //
//////////////////

params = add_param(params,"hessian_approximation","limited-memory");

/////////////////////
// Derivative Test //
/////////////////////

params = add_param(params,"derivative_test","first-order");

params = add_param(params,"mu-strategy","adaptive");

[x_sol, f_sol, extra] = ipopt(x0, 'f_C', 'df_C', 'g_C', 'dg_C', sparse_dg, dh, sparse_dh, ...
                              var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

printf('lambda = '); disp(extra('lambda'))

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

// unlink the optimization problem

ulink([f_C_link_handle, df_C_link_handle, g_C_link_handle, dg_C_link_handle]);
