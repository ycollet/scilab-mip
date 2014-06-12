function [xmin,fmin,status,extra] = glpk(c,a,lhs,rhs,lb,ub,ctype,vartype,param)
// Scilab interface for the GLPK library
//
// [xopt, fmin, status, extra] = glpk (c, a, lhs, rhs, lb, ub, ctype, vartype, param)
//
// Solve an LP/MILP problem using the GNU GLPK library. Given three
// arguments, glpk solves the following standard LP:
// 
// min C'*x subject to A*x  <= b
//
// but may also solve problems of the form
// 
// [ min | max ] C'*x
// subject to
//   A*x [ "=" | "<=" | ">=" ] b
//   x >= LB
//   x <= UB
//
// Input arguments:
// c = A column array containing the objective function coefficients.
// 
// A = A matrix containing the constraints coefficients.
// 
// b = A column array containing the right-hand side value for each constraint
//     in the constraint matrix.
// 
// lb = An array containing the lower bound on each of the variables.  If
//      lb is not supplied (or an empty array) the default lower bound for the variables is
//      minus infinite.
// 
// ub = An array containing the upper bound on each of the variables.  If
//      ub is not supplied (or an empty array) the default upper bound is assumed to be
//      infinite.
// 
// ctype = An array of characters containing the sense of each constraint in the
//         constraint matrix.  Each element of the array may be one of the
//         following values
//           'F' Free (unbounded) variable (the constraint is ignored).
//           'U' Variable with upper bound ( A(i,:)*x <= b(i)).
//           'S' Fixed Variable (A(i,:)*x = b(i)).
//           'L' Variable with lower bound (A(i,:)*x >= b(i)).
//           'D' Double-bounded variable (A(i,:)*x >= -b(i) and A(i,:)*x <= b(i)).
//  
// vartype = A column array containing the types of the variables.
//               'C' Continuous variable.
//               'I' Integer variable
//               'B' Binary variable
//
// sense = If sense is 1, the problem is a minimization.  If sense is
//         -1, the problem is a maximization.  The default value is 1.
// 
// param = A structure containing the following parameters used to define the
//         behavior of solver.  Missing elements in the structure take on default
//         values, so you only need to set the elements that you wish to change
//         from the default.
// 
//         Integer parameters:
//           msglev (default: 1) 
//                  Level of messages output by solver routines:
//                   0 - No output.
//                   1 - Error messages only.
//                   2 - Normal output.
//                   3 - Full output (includes informational messages).
// 
//           scale (default: 1). Scaling option: 
//                   0 - No scaling.
//                   1 - Equilibration scaling.
//                   2 - Geometric mean scaling, then equilibration scaling.
//
//           dual (default: 0). Dual simplex option:
//                   0 - Do not use the dual simplex.
//                   1 - If initial basic solution is dual feasible, use
//                       the dual simplex.
//
//           price (default: 1). Pricing option (for both primal and dual simplex):
//                   0 - Textbook pricing.
//                   1 - Steepest edge pricing.
//   
//           round (default: 0). Solution rounding option:
// 
//                   0 - Report all primal and dual values "as is".
//                   1 - Replace tiny primal and dual values by exact zero.
//
//           itlim (default: -1). Simplex iterations limit.  
//                 If this value is positive, it is decreased by one each
//                 time when one simplex iteration has been performed, and
//                 reaching zero value signals the solver to stop the search. 
//                 Negative value means no iterations limit.
// 
//           itcnt (default: 200). Output frequency, in iterations.  
//                 This parameter specifies how frequently the solver sends 
//                 information about the solution to the standard output.
// 
//           branch (default: 2). Branching heuristic option (for MIP only):
//                   0 - Branch on the first variable.
//                   1 - Branch on the last variable.
//                   2 - Branch using a heuristic by Driebeck and Tomlin.
//
//           btrack (default: 2). Backtracking heuristic option (for MIP only):
//                   0 - Depth first search.
//                   1 - Breadth first search.
//                   2 - Backtrack using the best projection heuristic.
// 
//           presol (default: 1). If this flag is set, the routine 
//                  lpx_simplex solves the problem using the built-in LP presolver. 
//                  Otherwise the LP presolver is not used.
//        
//           usecuts (default: 1). If this flag is set, the
//                  routine lpx_intopt generates and adds cutting planes to
//                  the MIP problem in order to improve its LP relaxation
//                  before applying the branch&bound method (Only Gomory's
//                  mixed integer cuts are implemented).
// 
//           lpsolver (default: 1) Select which solver to use.
//                    If the problem is a MIP problem this flag will be ignored.
//                       1 - Revised simplex method.
//                       2 - Interior point method.
//
//           save (default: 0). If this parameter is nonzero save a copy of 
//                the original problem to file. You can specify the 
//                file name and format by using the 'savefilename' and 'savefiletype' 
//                parameters (see in String Parameters Section here below).
//                If previous parameters are not defined the original problem 
//                is saved with CPLEX LP format in the default file "outpb.lp".
//
//           mpsinfo (default: 1) If this is set, 
//                   the interface writes to file several comment cards, 
//                   which contains some information about the problem. 
//                   Otherwise the routine writes no comment cards. 
//
//           mpsobj ( default: 2) This parameter tells the 
//                  routine how to output the objective function row: 
//                    0 - never output objective function row 
//                    1 - always output objective function row 
//                    2 - output objective function row if the problem has 
//                        no free rows 
//
//           mpsorig (default: 0) If this is set, the 
//                   routine uses the original symbolic names of rows and 
//                   columns. Otherwise the routine generates plain names 
//                   using ordinal numbers of rows and columns.
//
//           mpswide (default: 1) If this is set, the 
//                   routine uses all data fields. Otherwise the routine 
//                   keeps fields 5 and 6 empty. 
//
//           mpsfree (default: 0) If this is set, the routine 
//                   omits column and vector names every time when possible 
//                   (free style). Otherwise the routine never omits these 
//                   names (pedantic style). 
//
// 
//         Real parameters:
//           relax (default: 0.07). Relaxation parameter used 
//                 in the ratio test. If it is zero, the textbook ratio test 
//                 is used. If it is non-zero (should be positive), Harris'
//                 two-pass ratio test is used. In the latter case on the 
//                 first pass of the ratio test basic variables (in the case 
//                 of primal simplex) or reduced costs of non-basic variables 
//                 (in the case of dual simplex) are allowed to slightly violate 
//                 their bounds, but not more than relax*tolbnd or relax*toldj 
//                 (thus, relax is a percentage of tolbnd or toldj).
// 
//           tolbnd (default: 10e-7). Relative tolerance used 
//                  to check ifthe current basic solution is primal feasible.
//                  It is not recommended that you change this parameter 
//                  unless you have a detailed understanding of its purpose.
// 
//           toldj (default: 10e-7). Absolute tolerance used to 
//                 check if the current basic solution is dual feasible.  It 
//                 is not recommended that you change this parameter unless 
//                 you have a detailed understanding of its purpose.
// 
//           tolpiv (default: 10e-9). Relative tolerance used 
//                  to choose eligible pivotal elements of the simplex table.
//                  It is not recommended that you change this parameter 
//                  unless you have a detailed understanding of its purpose.
// 
//           objll ( default: -DBL_MAX). Lower limit of the 
//                 objective function. If on the phase II the objective
//                 function reaches this limit and continues decreasing, the
//                 solver stops the search. This parameter is used in the 
//                 dual simplex method only.
// 
//           objul (default: +DBL_MAX). Upper limit of the 
//                 objective function. If on the phase II the objective
//                 function reaches this limit and continues increasing, 
//                 the solver stops the search. This parameter is used in 
//                 the dual simplex only.
// 
//           tmlim (default: -1.0). Searching time limit, in 
//                 seconds. If this value is positive, it is decreased each 
//                 time when one simplex iteration has been performed by the
//                 amount of time spent for the iteration, and reaching zero 
//                 value signals the solver to stop the search. Negative 
//                 value means no time limit.
// 
//           outdly (default: 0.0). Output delay, in seconds. 
//                  This parameter specifies how long the solver should 
//                  delay sending information about the solution to the standard
//                  output. Non-positive value means no delay.
// 
//           tolint (default: 10e-5). Relative tolerance used 
//                  to check if the current basic solution is integer
//                  feasible. It is not recommended that you change this 
//                  parameter unless you have a detailed understanding of 
//                  its purpose.
// 
//           tolobj (default: 10e-7). Relative tolerance used 
//                  to check if the value of the objective function is not 
//                  better than in the best known integer feasible solution.  
//                  It is not recommended that you change this parameter 
//                  unless you have a detailed understanding of its purpose.
//
//         String Parameters:
//           savefilename (default: "outpb"). Specify the name to use to 
//                        save the original problem. MEX interface looks for 
//                        this parameter if 'save' parameter is set to 1. If 
//                        no name is provided "outpb" will be used.
//           savefiletype (default: CPLEX format). Specify the format type
//                        used to save the file. Only the following options
//                        are allowed:
//                          'fixedmps' - fixed MPS format (.mps).
//                          'freemps'  - free MPS format (.mps). 
//                          'cplex'    - CPLEX LP format (.lp).
//                          'plain'    - plain text (.txt).
// 
// Output values:
// xopt = The optimizer (the value of the decision variables at the optimum).
//
// fopt = The optimum value of the objective function.
//
// status = Status of the optimization.
//               1 solution is undefined
//               2 solution is feasible
//               3 solution is infeasible
//               4 no feasible solution exists
//               5 solution is optimal
//               6 solution is unbounded
//
//          If an error occurs, status will contain one of the following
//          codes.
//          Simplex method:
//               101  invalid basis
//               102  singular matrix
//               103  ill-conditioned matrix
//               104  invalid bounds
//               105  solver failed
//               106  objective lower limit reached
//               107  objective upper limit reached
//               108  iteration limit exceeded
//               109  time limit exceeded
//               110  no primal feasible solution
//
//          Interior point method, mixed integer problem:
//               204  Unable to start the search.
//               205  Objective function lower limit reached.
//               206  Objective function upper limit reached.
//               207  Iterations limit exhausted.
//               208  Time limit exhausted.
//               209  No feasible solution.
//               210  Numerical instability.
//               211  Problems with basis matrix.
//               212  No convergence (interior).
//               213  No primal feasible solution (LP presolver).
//               214  No dual feasible solution (LP presolver).
// 
// extra = A data structure containing the following fields:
//           lambda - Dual variables.
//           redcosts - Reduced Costs.
//           time - Time (in seconds) used for solving LP/MIP problem.
//           mem - Memory (in Kbytes) used for solving LP/MIP problem.
// 
// Example:
// 
// c = [10, 6, 4]';
// a = [ 1, 1, 1;
//      10, 4, 5;
//       2, 2, 6];
// b = [100, 600, 300]';
// lb = [0, 0, 0]';
// ub = [];
// ctype = "UUU";
// vartype = "CCC";
// s = -1;
// 
// param.msglev = 1;
// param.itlim = 100;
// 
// [xmin, fmin, status, extra] = glpk(c, a, b, lb, ub, ctype, vartype, s, param);
//
// Copyright 2008 Yann Collette
// Copyright 2005-2007 Nicolo' Giorgetti
// Email: Nicolo' Giorgetti <giorgetti __at __ ieee.org>
// Email: Yann Collette <ycollet __at __ freesurf.fr>

// This file is part of scicoinor.
//
// The GLPK Scilab interface is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// This part of code is distributed with the FURTHER condition that it 
// can be linked to the Matlab libraries and/or use it inside the Matlab 
// environment.
//
// The GLPK Scilab interface is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with SCIGLPK; see the file COPYING.  If not, write to the Free
// Software Foundation, 59 Temple Place - Suite 330, Boston, MA
// 02111-1307, USA.

[nargout,nargin] = argn();

// If there is no input output the version and syntax
if (nargin < 3 | nargin > 9) then
  printf('warning: some parameters must be defined\n');
  printf('Syntax: [xmin,fmin,status,extra]=glpk(c,a,lhs,rhs,lb,ub,ctype,vartype,param)\n');
  return;
end

if (and(size(c) > 1) | ~isreal(c) | type(c)==10) then
  error('C must be a real vector');
end
nx = length (c);
// Force column vector.
c = c(:);

// 2) Matrix constraint
if (isempty(a)) then
  error('A cannot be an empty matrix');
end
[nc, nxa] = size(a);
if (~isreal(a) | nxa ~= nx) then
  error(sprintf('A must be a real valued %d by %d matrix', nc, nx));
end

// 3) LHS
if (isempty(lhs)) then
  error('LHS cannot be an empty vector');
end
if (~isreal(lhs) | length(lhs) ~= nc) then
  error(sprintf('LHS must be a real valued %d by 1 vector', nc));
end

// 4) RHS
if (isempty(rhs)) then
  error('RHS cannot be an empty vector');
end
if (~isreal(rhs) | length(rhs) ~= nc) then
  error(sprintf('RHS must be a real valued %d by 1 vector', nc));
end

// 5) Vector with the lower bound of each variable
if (nargin>=5) then
  if (isempty(lb))
    lb = -%inf*ones(nx, 1); 
  elseif (~isreal(lb) | and(size(lb) > 1) | length(lb) ~= nx) then
    error(sprintf('LB must be a real valued %d by 1 column vector', nx));
  end
else
  lb = -%inf*ones(nx, 1);
end

// 6) Vector with the upper bound of each variable
if (nargin>=6) then
  if (isempty(ub))
    ub = %inf*ones(nx, 1);
  elseif (~isreal(ub) | and(size(ub) > 1) | length(ub) ~= nx) then
    error(sprintf('UB must be a real valued %d by 1 column vector', nx));
  end
else
  ub = %inf*ones(nx, 1);
end

// 7) Sense of each constraint
if (nargin>=7) then
  if max(size(ctype))~=1 then
    ctype = strcat(ctype);
  end
  if (isempty (ctype)) then
    ctype = asciimat(ascii('U')*ones(nc,1));
  elseif (type(ctype)~=10 | length(ctype)~=nc) then
    error(sprintf('CTYPE must be a char valued vector of length %d', nc));
  else
    for i=1:nc
      if ~or(part(ctype,i)==['l','L','e','E','g','G','r','R','n','N']) then
        error(sprintf('CTYPE must contain only L, E, G, R, or N'));
      end
    end
  end
else
  ctype = ascii(ascii('U')*ones(nc,1));
end

// 8) Vector with the type of variables
if (nargin>=8) then
  if max(size(vartype))~=1 then
    vartype = strcat(vartype);
  end
  if isempty(vartype)
    vartype = asciimat(ascii('C')*ones(nx, 1));
  elseif (type(vartype)~=10 | length(vartype)~=nx) then
    error(sprintf('VARTYPE must be a char valued vector of length %d', nx));
  else
    for i=1:nx
      if ~or(part(vartype,i)==['c','C','i','I','b','B']) then
        error(sprintf('VARTYPE must contain only C, I or B'));
      end
    end
  end
else
  // As default we consider continuous vars
  vartype = ascii(ascii('C')*ones(nx, 1));
end

// 9) Parameters vector
if (nargin>=9) then
  if (typeof(param)~='plist') then
    error('PARAM must be a parameters list (plist)');
  end
else
  param = init_param();
end

c   = full(c);
lhs = full(lhs);
rhs = full(rhs);
ub  = full(ub);
lb  = full(lb);

[xmin, fmin, status, extra] = sciglpk(c, a, lhs, rhs, lb, ub, ctype, vartype, param);
endfunction
