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


// Check that we can use numerical derivatives to compute the 
// derivative of the functions

////////////////////////////////////////////////////////////////////////
// Inputs
function f = objfun ( x )
pause
  f = exp(x(1))*(4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 2*x(2) + 1);
endfunction
function [c, ceq] = confun(x)
pause
  // Nonlinear inequality constraints
  c = [
    1.5 + x(1)*x(2) - x(1) - x(2)
    -x(1)*x(2) - 10
  ];
  // Nonlinear equality constraints
  ceq = [];
endfunction
////////////////////////////////////////////////////////////////////////
// Transform the problem into a ipopt input


function y = ipopt_df ( x )
pause
  y = derivative ( objfun , x )
endfunction

// The constraints
function y = ipopt_g ( x )
pause
  [c, ceq] = confun(x)
  y = c
endfunction

function y = ipopt_dg ( x )
pause
  y = derivative ( ipopt_g , x )
  y = matrix ( y' , 4 , 1 )
endfunction

function y = ipopt_gi ( x , i )
pause
  [c, ceq] = confun(x)
  y = c(i)
endfunction

// The Lagrangian
function y = ipopt_dh ( x , lambda , obj_weight )
pause
  [Jf , Hf] = derivative(objfun,x,%eps^(1/3),1,"blockmat" )
  [Jg1 , Hg1] = derivative(list(ipopt_gi,1),x,%eps^(1/3),1,"blockmat" )
  [Jg2 , Hg2] = derivative(list(ipopt_gi,2),x,%eps^(1/3),1,"blockmat" )
  y = obj_weight * Hf + lambda(1) * Hg1 + lambda(2) * Hg2
endfunction
// Make a starting guess at the solution
x0 = [-1 1]';
xopt = [-9.5474  1.0474]';
if ( %f ) then
  // Check that the functions perform well
  y = ipopt_df ( x0 )
  y = ipopt_g ( x0 )
  y = ipopt_dg ( x0 )
  y = ipopt_dh ( x0,[1 2],1)
end
// The sparsity structure of the constraints
sparse_dg = [
  1 1
  1 2
  2 1
  2 2
];

// The sparsity structure of the Lagrangian
sparse_dh = [
  1 1
  1 2
  2 1
  2 2
];

upper = [%inf %inf]';
lower = [-%inf -%inf]';
nb_constr = 2;
var_lin_type = [1 1]'; // Non-Linear
constr_lin_type = [1 1]'; // Non-Linear
constr_rhs = [0 0]';
constr_lhs = [-%inf -%inf]';

////////////////////////////////////////////////////////////////////////
[x_sol, f_sol, extra] = ipopt(x0, objfun, ipopt_df, ipopt_g, ..
  ipopt_dg, sparse_dg, ipopt_dh, sparse_dh, var_lin_type, ..
  constr_lin_type, constr_rhs, constr_lhs, lower, upper);

// TODO : make a nicer message than "Scilab has found a critical error (EXCEPTION_ACCESS_VIOLATION)"
// TODO : make the test run again



