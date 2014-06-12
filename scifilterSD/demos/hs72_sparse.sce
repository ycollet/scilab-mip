function f = function_hs72d(x)
  f = 1.0 + 1.0/x(1) + 1.0/x(2) + 1.0/x(3) + 1.0/x(4);
endfunction

function g = gradient_hs72d(x)
  g(1) = -1.0/x(1)^2;
  g(2) = -1.0/x(2)^2;
  g(3) = -1.0/x(3)^2;
  g(4) = -1.0/x(4)^2;
endfunction

a = list();
a(1) = [4.0 2.25 1.0 0.25];
a(2) = [0.16 0.36 0.64 0.64];

la = list();
la(1) = [1 2 3 4];
la(2) = [1 2 3 4];

ainfty = 1.0e20;
tol    = 1.0e-12;

x0           = [];
lower_bounds = [];
upper_bounds = [];

// Upper and lower bounds for the variables
for i=1:4
  x0(i)           = 1.0;
  lower_bounds(i) = 1.0/((5 - i)*1.0e5);
  upper_bounds(i) = 1.0e3;
end

// Upper and lower bounds for the constraints
lower_bounds(5) = -ainfty;
lower_bounds(6) = -ainfty;
upper_bounds(5) = 4.01e-2;
upper_bounds(6) = 1.0085e-2;

params_in = init_param();
params_in = add_param(params_in, 'rgtol',  1.0e-5);
params_in = add_param(params_in, 'ainfty', 1.0e20);
params_in = add_param(params_in, 'ubd',    1.0e5);
params_in = add_param(params_in, 'fmin',  -ainfty);
params_in = add_param(params_in, 'iprint', 1);
params_in = add_param(params_in, 'kmax',   4);
params_in = add_param(params_in, 'maxg',   5);
params_in = add_param(params_in, 'mlp',    50);
params_in = add_param(params_in, 'mode',   0);
params_in = add_param(params_in, 'mxgr',   100);
params_in = add_param(params_in, 'mxws',   3000);
params_in = add_param(params_in, 'mxlws',  3000);
params_in = add_param(params_in, 'nout',   0);

// Solution
// x(1) = 0.5170432e-2;
// x(2) = 0.5569570e-2;
// x(3) = 0.5404878e-2;
// x(4) = 0.5927444e-2;

[x_out, f_out, ifail, params_out] = glcpd_sparse(x0,function_hs72d, gradient_hs72d, a, la, lower_bounds, upper_bounds, params_in);

// ifail_out   outcome of the process
// 0   = solution obtained
// 1   = unbounded problem (f(x)<fmin has occurred: note grad is not evaluated in this case)
// 2   = bl(i) > bu(i) for some i
// 3   = infeasible problem detected in Phase 1
// 4   = line search cannot improve f (possibly increase rgtol)
// 5   = mxgr gradient calls exceeded (this test is only carried out at the start of each iteration)
// 6   = incorrect setting of m, n, kmax, maxg, mlp, mode or tol
// 7   = not enough space in ws or lws
// 8   = not enough space in lp (increase mlp)
// 9   = dimension of reduced space too large (increase kmax)
// 10  = maximum number of unsuccessful restarts taken
// >10 = possible use by later sparse matrix codes

printf('ifail = %d\n',ifail);
printf('f_out = %f\n', f_out);
printf('x_out = '); disp(x_out);
