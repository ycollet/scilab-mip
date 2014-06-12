function [xmin,fmin,status,extra] = lpsolve(c,A,lhs,rhs,lb,ub,btype,vartype,options,special)
// LPSOLVE Interface to MILP solver LPSOLVE
//
// min c'*x
// s.t   lhs < A*x < rhs,
//       lb  < x < ub
//
// Options structure (see LPSOLVE user manual for details)
// Author Yann Collette

[nargout,nargin] = argn();

//
// Check input
//

if ~isdef('options','local') then option  = init_param(); end
if ~isdef('vartype','local') then vartype = []; end
if ~isdef('ub','local')      then ub      = []; end
if ~isdef('lb','local')      then lb      = []; end
if ~isdef('btype','local')   then btype   = []; end
if ~isdef('rhs','local')     then rhs     = []; end
if ~isdef('lhs','local')     then lhs     = []; end
if ~isdef('A','local')       then A       = []; end
if ~isdef('c','local')       then c       = []; end
if ~isdef('special','local') then special = init_constraint(); end

xmin   = [];
lambda = [];
status = -1;
fmin   = [];

if nargin < 1 then
  printf('warning: some parameters must be defined\n');
  printf('function [xmin,fmin,lambda,status] = lpsolve(c,A,lhs,rhs,lb,ub,btype,vartype,options,special)\n');
  return;
end

if isempty(c) then
  printf('error: the c matrix must not be empty\n');
  return;
end

if length(lb)~=length(ub) | length(lb)~=length(c) then
  printf('error: dimension mismatch between lb, ub and / or c\n');
  return
end

if or(lb>ub) then
  printf('error: a lower variable boundary is higher than an upper variable boundary\n');
  return;
end

if or(lhs>rhs) then
   printf('error: a lower constraint boundary is higher than an upper constraint boundary\n');
   return;
end

if size(A,1)~=length(rhs) | size(A,2)~=length(ub) then
  printf('error: problem with the dimensions of the constraint matrix\n');
  return;
end

if isempty(btype) then
  // 'L' - smaller than - <=
  // 'E' - equality     - =
  // 'G' - greater than - >=
  // 'R' - Range        - <= + >=
  // 'N' - Free         - no constraints
  dimension = length(b);
  btype = ascii(ascii('L')*ones(dimension,1));
end

if isempty(vartype) then
  dimension = length(lb);
  vartype = ascii(ascii('C')*ones(dimension,1));
end

//
// LPSOLVE sparse format
//

c   = full(c);
rhs = full(rhs);
lhs = full(lhs);
lb  = full(lb);
ub  = full(ub);

//
// Call Scilab LPSOLVE interface
//

if (typeof(special)=='clist') then
  [xmin,fmin,status,extra] = scilpsolve(c,A,lhs,rhs,lb,ub,btype,vartype,options,special);
elseif (isempty(special)) then
  special = init_constraint();
  [xmin,fmin,status,extra] = scilpsolve(c,A,lhs,rhs,lb,ub,btype,vartype,options,special);   
else
  error('lpsolve: the special parameter must be a clist or must be empty');
end
endfunction

