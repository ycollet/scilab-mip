// Copyright (C) 2010 - DIGITEO - Michael Baudin
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

// <-- JVM NOT MANDATORY -->
// <-- ENGLISH IMPOSED -->


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


// Definition of the optimization problem
// The objective function
function y=f(x,x_new)
  y=4*x(1) - x(2)^2 - 12
endfunction

function y=df(x,x_new)
  y(1) = 4;
  y(2) = -2*x(2);
endfunction

// The constraints
function y=g(x,x_new)
  y(1) = - 10*x(1) + x(1)^2 - 10*x(2) + x(2)^2 + 34; 
  y(2) = 20 - x(1)^2 - x(2)^2
endfunction
function y=dg(x,x_new)
  y(1) = -10 + 2*x(1);
  y(2) = -10 + 2*x(2); 
  y(3) = -2*x(1); 
  y(4) = -2*x(2);
endfunction

// The sparsity structure of the constraints
sparse_dg = [1 1; ...
             1 2; ...
             2 1; ...
             2 2];

// The Lagrangian
function y = dh(x,lambda,obj_weight,x_new,lambda_new)
  y(1) = lambda(1)*2 - lambda(2)*2;
  y(2) = -obj_weight*2 + lambda(1)*2 - lambda(2)*2;
endfunction

// The sparsity structure of the Lagrangian
sparse_dh = [1 1; ...
             2 2];

upper = [15;15];
lower = [-15;-15];
x0    = [-12;-12]; // Feasible starting point

var_lin_type(1) = 1; // Non-Linear
var_lin_type(2) = 1; // Non-Linear (this variable appears nonlinearly in the objective function or at least in one constraint)

constr_lin_type (1) = 1; // Non-Linear
constr_lin_type (2) = 1; // Non-Linear

constr_rhs(1) = 0;
constr_rhs(2) = 0;
constr_lhs(1) = -10000;
constr_lhs(2) = 0;

////////////////////////////////////////////////////////////////////////

params = init_param();
// We use the given Hessian
//params = add_param(params,"hessian_approximation","exact");
// We use a limited-bfgs approximation for the Hessian.
params = add_param(params,"hessian_approximation","limited-memory");

[x_sol, f_sol, extra] = ipopt(x0, f, df, g, dg, sparse_dg, dh, sparse_dh, ..
  var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

assert_close ( x_sol , [1.0537922  4.3462078]' , 1.e-7 );
assert_close ( f_sol , -26.67435344 , 1.e-7 );
lambda = extra('lambda')
assert_close ( lambda , [0 0] , 1.e-100 );
status = extra("status");
assert_equal ( status , 29 );
iter_count = extra("it_count");
assert_equal ( iter_count , 25 );
cpu_time = extra("cpu_time");
assert_equal ( cpu_time > 0 , %t );
fobj_eval = extra("fobj_eval");
// TODO : fix that bug : I got 673720360 !!!
assert_equal ( fobj_eval < 10000 , %t ); 
fobj_grad_eval = extra("fobj_grad_eval");
assert_equal ( fobj_grad_eval , 35 );
constr_eval = extra("constr_eval");
assert_equal ( constr_eval , 2 );
constr_jac_eval = extra("constr_jac_eval");
hess_eval = extra("hess_eval");
assert_equal ( hess_eval , 1 );
dual_inf = extra("dual_inf");
assert_close ( dual_inf , 0 , 1.e-100 );
// TODO : constr_viol is zero, but constr_viol == 0 is false !!!
constr_viol = extra("constr_viol"); 
assert_equal ( constr_viol , 0 );
complementarity = extra("complementarity"); 
// TODO : complementarity is zero, but complementarity == 0 is false !!!
assert_equal ( complementarity , 0 ); 
kkt_error = extra("kkt_error");
// TODO : kkt_error is zero, but kkt_error == 0 is false !!!
assert_equal ( kkt_error , 0 ); 


