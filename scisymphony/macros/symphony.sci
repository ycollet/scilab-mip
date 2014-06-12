function [xmin,fmin,status,extra] = symphony(c1,c2,A,lhs,rhs,lb,ub,btype,vartype,options)
// Symphony Interface to MILP solver Symphony
//
// [x,z,status] = symphony(c1,c2,A,lhs,rhs,lb,ub,btype,vartype,options)
//
// min c1'*x
// s.t   lhs < A*x < rhs,
//       lb  < x < ub
//
// or
//
// min [c1'*x c2'x]
// s.t   lhs < A*x < rhs,
//       lb  < x < ub

//
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
if ~isdef('c1','local')      then c1      = []; end
if ~isdef('c2','local')      then c2      = []; end

x      = [];
lambda = [];
status = -1;

if nargin < 1 then
  printf('warning: some parameters must be defined\n');
  printf('function [xmin,fmin,lambda,status] = symphony(c1,c2,A,lhs,rhs,lb,ub,btype,vartype,options)\n');
  return;
end

if isempty(c1) then
  printf('error: the c1 matrix must not be empty\n');
  return;
end

if length(lb)~=length(ub) | length(lb)~=length(c1) then
  printf('error: dimension mismatch between lb, ub and / or c\n');
  return
end

if or(lb>ub) then
  printf('error: a lower variable boundary is higher than an upper variable boundary\n');
  return;
end

if length(lhs)~=length(rhs) then
  printf('error: dimension mismatch between rhs and lhs\n');
  return;
end

if or(lhs>rhs) then
   printf('error: a lower constraint boundary is higher than an upper constraint boundary\n');
   return;
end

if size(A,1)~=length(rhs) | size(A,2) ~=length(ub) then
  printf('error: problem with the dimensions of the constraint matrix\n');
  return;
end

if isempty(btype) then
  // 'L' - smaller than - <=
  // 'E' - equality     - =
  // 'G' - greater than - >=
  // 'R' - Range        - <= + >=
  // 'N' - Free         - no constraints
  dimension = max(length(rhs),length(lhs));
  btype = ascii(ascii('L')*ones(dimension,1));
end

btype = convstr(btype,'u');

if isempty(vartype) then
  dimension = length(lb);
  vartype = ascii(ascii('C')*ones(dimension,1));
end

vartype = convstr(vartype,'u');

//
// deal with rhs and lhs
//

rhs2 = rhs;

// L part
Index = strindex(btype,'L');
rhs2(Index) = rhs(Index);
// G part
Index = strindex(btype,'G');
rhs2(Index) = lhs(Index);
// E part
Index = strindex(btype,'E');
rhs2(Index) = lhs(Index);

rhs = rhs2;

clear rhs2;

//
// Symphony sparse format
//

c1  = full(c1);
c2  = full(c2);
lhs = full(lhs);
rhs = full(rhs);
lb  = full(lb);
ub  = full(ub);

//
// Call Scilab Symphony interface
//

[xmin,fmin,status,extra] = scisymphony(A',c1,c2,lhs,rhs,ub,lb,btype,vartype,options);   
endfunction
