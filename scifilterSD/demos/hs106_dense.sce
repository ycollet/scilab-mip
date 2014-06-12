function [f, c] = function_hs106d(x)
  t = 25.0e-4;

  // objective function
  f = x(1) + x(2) + x(3);
  
  // constraint functions
  c(1) = t*(x(4) + x(6)) - 1.0;
  c(2) = t*(x(5) + x(7) - x(4)) - 1.0;
  c(3) = 1.0e-2*(x(8) - x(5)) - 1.0;
  c(4) = 1.0e2*x(1) + 833.33252*x(4) - x(1)*x(6) - 83333.333;
  c(5) = 125.0e1*(x(5) - x(4)) + x(2)*(x(4) - x(7));
  c(6) = 125.0e4 - 25.0e2*x(5) + x(3)*(x(5) - x(8));
endfunction

function a = gradient_hs106d(x)
  t = 25.0e-4;

  //  elements of grad.f and the Jacobian matrix
  a(1,1) = 1;
  a(2,1) = 1;
  a(3,1) = 1;
  a(4,1) = 0;
  a(5,1) = 0;
  a(6,1) = 0;
  a(7,1) = 0;
  a(8,1) = 0;
  
  a(1,2) = 0;
  a(2,2) = 0;
  a(3,2) = 0;
  a(4,2) = t;
  a(5,2) = 0;
  a(6,2) = t;
  a(7,2) = 0;
  a(8,2) = 0;

  a(1,3) = 0;
  a(2,3) = 0;
  a(3,3) = 0;
  a(4,3) = -t;
  a(5,3) = t;
  a(6,3) = 0;
  a(7,3) = t;
  a(8,3) = 0;

  a(1,4) = 0;
  a(2,4) = 0;
  a(3,4) = 0;
  a(4,4) = 0;
  a(5,4) = -1e-2;
  a(6,4) = 0;
  a(7,4) = 0;
  a(8,4) = 1e-2;

  a(1,5) = 1.0e2 - x(6);
  a(2,5) = 0;
  a(3,5) = 0;
  a(4,5) = 833.33252;
  a(5,5) = 0;
  a(6,5) = -x(1);
  a(7,5) = 0;
  a(8,5) = 0;

  a(1,6) = 0;
  a(2,6) = x(4) - x(7);
  a(3,6) = 0;
  a(4,6) = -125.0e1 + x(2);
  a(5,6) = 125.0e1;
  a(6,6) = 0;
  a(7,6) = -x(2);
  a(8,6) = 0;

  a(1,7) = 0;
  a(2,7) = 0;
  a(3,7) = x(5) - x(8);
  a(4,7) = 0;
  a(5,7) = -25.0e2 + x(3);
  a(6,7) = 0;
  a(7,7) = 0;
  a(8,7) = -x(3);
endfunction

lower_bounds = [100,  1000, 1000, 10,  10,  10,  10,  10];
upper_bounds = [10000,10000,10000,1000,1000,1000,1000,1000];
x0           = [5000, 5000, 5000, 200, 350, 150, 225, 425];
lambda       = 1e-2*[1,1,1,1,1,1,1,1];

// Solution and multipliers are:
// x_sol = [579.3167, 1359.943, 5110.071, 182.0174,  295.5985, 217.9799, 286.4162, 395.5979];
// lambda_sol = [1964.046, 5210.645, 5110.092, 8.475914e-03, 9.578792e-03, 0.01];

ainfty = 1e20;

// Initialisation of the bounds for the constraints
lower_bounds = [lower_bounds,   -ainfty, -ainfty, -ainfty, -ainfty, -ainfty, -ainfty];
upper_bounds = [upper_bounds,    0,       0,       0,       0,       0,       0];
//lambda       = [lambda, 1.0e-2*[1, 1, 1, 1, 1, 1]];
lambda       = [lambda, 0.0*[1, 1, 1, 1, 1, 1]];

params_in = init_param();

params_in = add_param(params_in, 'ainfty', 1.0e20);

params_in = add_param(params_in, 'fmin',  -ainfty);
params_in = add_param(params_in, 'ubd',    1.0e5);
params_in = add_param(params_in, 'mlp',    50);
params_in = add_param(params_in, 'mxf',    50);

params_in = add_param(params_in, 'rho',    1.0e2);
params_in = add_param(params_in, 'htol',   1.0e-6);
params_in = add_param(params_in, 'rgtol',  1.0e-4);
params_in = add_param(params_in, 'maxit',  60);
params_in = add_param(params_in, 'iprint', 0);
params_in = add_param(params_in, 'kmax',   length(x0));
params_in = add_param(params_in, 'maxg',   5);
//params_in = add_param(params_in, 'mxgr',   1000000);

[x_out, f_out, ifail, lambda_out, params_out] = filtersd_dense(x0,function_hs106d, gradient_hs106d, lower_bounds, upper_bounds, params_in);

//  ifail   returns failure indication as follows
//               0 = successful run
//               1 = unbounded NLP (f <= fmin at an htol-feasible point)
//               2 = bounds on x are inconsistent
//               3 = local minimum of feasibility problem and h > htol
//                   (nonlinear constraints are locally inconsistent)
//               4 = initial point x has h > ubd (reset ubd or x and re-enter)
//               5 = maxit major iterations have been carried out
//               6 = termination with rho <= htol
//               7 = not enough workspace in ws or lws (see message)
//               8 = insufficient space for filter (increase mxf and re-enter)

//              >9 = unexpected fail in LCP solver (10 has been added to ifail)

//              10 = solution obtained
//              11 = unbounded problem (f(x)<fmin has occurred: note grad is not
//                    evaluated in this case)
//              12 = bl(i) > bu(i) for some i
//              13 = infeasible problem detected in Phase 1
//              14 = line search cannot improve f (possibly increase rgtol)
//              15 = mxgr gradient calls exceeded (this test is only carried
//                    out at the start of each iteration)
//              16 = incorrect setting of m, n, kmax, maxg, mlp, m0de or tol
//              17 = not enough space in ws or lws
//              18 = not enough space in lp (increase mlp)
//              19 = dimension of reduced space too large (increase kmax)
//              20 = maximum number of unsuccessful restarts taken
//             >20= possible use by later sparse matrix codes

printf('ifail = %d\n',ifail)
printf('number of free variables = %d\n',params_out('k'));
printf('number of function and gradient calls = %d, %d\n',params_out('nft'),params_out('ngt'));
printf('cstype = %s\n', params_out('cstype'));
printf('f_out = %f\n', f_out);
printf('x_out = '); disp(x_out);
printf('lambda_out = '); disp(lambda_out);

[f_sol, c_sol] = function_hs106d(x_out);

printf('value of the constraints ='); disp(c_sol);

