lines(0);

f  = [];
df = [];
g  = [];
dg = [];
dh = [];
sparse_dg = [];
sparse_dh = [];

////////////////////////////////////////////////////////////////////////

// a = 5, b = 3
function y = f(x,x_new,a,b)
  y = a*x(1)^2-b*x(2)^2;
endfunction

function y = df(x,x_new,a,b)
  y(1) = 2*a*x(1);
  y(2) = -2*b*x(2);
endfunction

// a=-1, b=-1
function y = g(x,x_new,a,b)
  y(1) = a*x(1);
  y(2) = b*x(2);
endfunction

function y = dg(x,x_new,a,b)
  y(1) = a;
  y(2) = b;
endfunction
                  
sparse_dg = [1,1; ...
             2,2];
  
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

params = add_param(params,"journal_level",5);
params = add_param(params,"tol",1e-8);
params = add_param(params,"max_iter",3000);
params = add_param(params,"hessian_approximation","limited-memory");
params = add_param(params,"derivative_test","first-order");
params = add_param(params,"mu-strategy","adaptive");

[x_sol, f_sol, extra] = ipopt(x0, list(f,5,3), list(df,5,3), list(g,-1,-1), list(dg,-1,-1), sparse_dg, dh, sparse_dh, ...
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

