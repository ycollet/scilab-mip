// Copyright (C) 2010 - DIGITEO - Michael Baudin
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

// <-- JVM NOT MANDATORY -->
// <-- ENGLISH IMPOSED -->

// TODO : fix and complete this !


//
// assert_close --
//   Returns 1 if the two real matrices computed and expected are close,
//   i.e. if the relative distance between computed and expected is lesser than epsilon.
// Arguments
//   computed, expected : the two matrices to compare
//   epsilon : a small number
//
function flag = assert_close ( computed, expected, epsilon )
  if expected==0.0 then
    shift = norm(computed-expected);
  else
    shift = norm(computed-expected)/norm(expected);
  end
  if shift < epsilon then
    flag = 1;
  else
    flag = 0;
  end
  if flag <> 1 then pause,end
endfunction
//
// assert_equal --
//   Returns 1 if the two real matrices computed and expected are equal.
// Arguments
//   computed, expected : the two matrices to compare
//   epsilon : a small number
//
function flag = assert_equal ( computed , expected )
  if computed==expected then
    flag = 1;
  else
    flag = 0;
  end
  if flag <> 1 then pause,end
endfunction

function y = rosenbrock_f ( x , x_new )
    y = 100.0 *(x(2)-x(1)^2)^2 + (1-x(1))^2;
endfunction

function y = rosenbrock_df ( x , x_new )
    y(1) = -400*(x(2)-x(1)^2)*x(1) - 2*(1-x(1));
    y(2) = 200*(x(2)-x(1)^2);
endfunction

// Nonlinear inequality constraints
function y = rosenbrock_g ( x , x_new )
  y = x(1)^2 + x(2)^2 - 1.5;
endfunction

function y = rosenbrock_dg ( x , x_new )
  y(1) = 2 * x(1);
  y(2) = 2 * x(2);
endfunction

function y = rosenbrock_Hf ( x , x_new )
  y = zeros(2,2);
  y(1,1) = diag(-400*x(2) + 1200*x(1).^2 + 2);
  y(2,2) = 200;
  y = y - diag(400*x(1),1) - diag(400*x(1),-1);
endfunction

function y = rosenbrock_Hg ( x , x_new )
  y = [2 0; ...
       0 2];
endfunction

// Check derivatives of F
x = [-1.9 2.0]';
[df1 , Hf1] = derivative(list(rosenbrock_f,%t),x,H_form="blockmat");
df2 = rosenbrock_df ( x )';
Hf2 = rosenbrock_Hf ( x );
printf('Checking the derivatives of F:\n');
printf('norm(df1-df2) = %f norm(Hf1-Hf2) = %f\n', norm(df1-df2), norm(Hf1-Hf2));



// Check derivatives of G
x = [-1.9 2.0]';
[dg1 , Hg1] = derivative(rosenbrock_g,x,H_form="blockmat");
dg2 = rosenbrock_dg ( x )';
Hg2 = rosenbrock_Hg ( x );

printf('Checking the derivatives of G:\n');
printf('norm(dg1-dg2) = %f norm(Hg1-Hg2) = %f\n', norm(dg1-dg2), norm(Hg1-Hg2));

// The sparsity structure of the constraints
sparse_dg = [1 1; ..
 	     1 2];

// The Hessian of the Lagrangian
function y = rosenbrock_hessian ( x , lambda , obj_weight , x_new , lambda_new )
  Hf = rosenbrock_Hf(x);
  Hg = rosenbrock_Hg(x);
  y  = obj_weight * Hf + lambda(1) * Hg;
endfunction

// The sparsity structure of the Lagrangian
sparse_dh = [1 1; ...
             1 2; ...
             2 1; ...
             2 2];

upper = [];
lower = [];

// Not Feasible starting point
x0                  = [-1.9 2.0]';
nb_constr           = 1;
var_lin_type(1)     = 1; // Non-Linear
constr_lin_type (1) = 1; // Non-Linear
constr_rhs(1)       = 0;
constr_lhs(1)       = -%inf;

////////////////////////////////////////////////////////////////////////

params = init_param();
// We use the given Hessian
params = add_param(params,"hessian_approximation","exact");
// We use a limited-bfgs approximation for the Hessian.
//params = add_param(params,"hessian_approximation","limited-memory");

params = add_param(params,"journal_level",5);

[x_sol, f_sol, extra] = ipopt(x0, rosenbrock_f, rosenbrock_df,  ...
                                  rosenbrock_g, rosenbrock_dg, sparse_dg, ...
				  rosenbrock_hessian, sparse_dh, ...
				  var_lin_type, constr_lin_type, ...
				  constr_rhs, constr_lhs, ...
				  lower, upper, params);
