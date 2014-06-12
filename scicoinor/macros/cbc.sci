function [xmin,fmin,status,extra] = cbc(Q,c,A,lhs,rhs,lb,ub,btype,vartype,options,special)
// CBC Interface to the MIP solver CBC
//
// [x,z,status] = cbc(Q,c,A,lhs,rhs,lb,ub,btype,vartype,options)
//
// min c'*x
// s.t   lhs < A*x < rhs,
//       lb  < x < ub
//
// Options structure (see CBC user manual for details)
//    options.maxnumiterations [int>=0 (default 99999999)]
//    options.maxnumseconds    [int>=0 (default 3600)]
//    options.primaltolerance  [double>=0 (default 1e-7)]
//    options.dualtolerance    [double>=0 (default 1e-7)]
//    options.verbose          [0|1|... (default 0)]
//     - 0 - none
//     - 1 - just final
//     - 2 - just factorizations
//     - 3 - as 2 plus a bit more 
//     - 4 - verbose above that
//     - 8,16,32 etc just for selective debug. 
//    options.solver
//     - 1 - primal (default)
//     - 2 - dual
//     - 3 - barrier
//     - 4 - barrier no cross
//     - 5 - reduced gradient
//     - 6 - primal dual interior point
//     - 7 - pdco
//    options.optim_dir        [1 (minimize) | -1 (maximize) | 0 (ignore), (default 1)].
//    options.writemps         [filename].
//    options.perturb:
//     - 50 - switch on perturbation 
//     - 100 - auto perturb if takes too long (1.0e-6 largest nonzero) 
//     - 101 - we are perturbed 
//     - 102 - don't try perturbing again 
//     - default is 100
//     - others are for playing
//    options.fact_freq        [integer].
//    options.presolve         [0 (false) | 1 (true), (default 0)].
//    options.postsolve        [1 (false) | 2 (true), (default 0)].
//
// output
//  x      : primal
//  z      : dual
//  status : 0 - optimal, 1 - infeasible, 2- unbounded

// Author: Yann Collette

[nargout,nargin] = argn();

//
// Check input
//

if ~isdef('options','local')      then option       = init_param(); end
if ~isdef('vartype','local')      then vartype      = []; end
if ~isdef('ub','local')           then ub           = []; end
if ~isdef('lb','local')           then lb           = []; end
if ~isdef('btype','local')        then btype        = []; end
if ~isdef('rhs','local')          then rhs          = []; end
if ~isdef('lhs','local')          then lhs          = []; end
if ~isdef('A','local')            then A            = []; end
if ~isdef('c','local')            then c            = []; end
if ~isdef('Q','local')            then Q            = []; end
if ~isdef('special','local')      then special      = init_constraint(); end
  
x      = [];
lambda = [];
status = -1;

if nargin < 1 then
  printf('warning: some parameters must be defined\n');
  printf('function [xmin,fmin,lambda,status] = cbc(Q,c,A,lhs,rhs,lb,ub,btype,vartype,options)\n');
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

if isempty(vartype) then
  dimension = length(lb);
  vartype = ascii(ascii('C')*ones(dimension,1));
end

//
// CBC sparse format
//

c   = full(c);
lhs = full(lhs);
rhs = full(rhs);
lb  = full(lb);
ub  = full(ub);

if isempty(options) then
  options = init_param();
  options = add_param(options,'maxnumiterations',10000000);
  options = add_param(options,'maxnumseconds',1000);
  options = add_param(options,'primaltolerance',1e-7);
  options = add_param(options,'dualtolerance',1e-7);
  options = add_param(options,'verbose',0);
  options = add_param(options,'optim_dir', 1); // optimisation direction: 1 - minimize, -1 - maximize, 0 - ignore

  ///////////////////////////
  // Set the compare class //
  ///////////////////////////

  ////////////////////
  // cbccompareuser //
  ////////////////////
  options = add_param(options,'cbccompareuser',0.0); // from sample1.cpp
  
  ///////////////////////////////
  // Set the Cgl cut generator //
  ///////////////////////////////

  ///////////////////
  // cglpreprocess //
  ///////////////////
  options = add_param(options,'cglpreprocess',1);

  ////////////////
  // cglprobing //
  ////////////////
  options = add_param(options,'cglprobing',1);

  ///////////////
  // cglgomory //
  ///////////////
  // Description of the method:
  options = add_param(options,'cglgomory',1);

  //////////////////////
  // cglknapsackcover //
  //////////////////////
  options = add_param(options,'cglknapsackcover',1);

  /////////////////
  // cglredsplit //
  /////////////////
  options = add_param(options,'cglredsplit',1);

  ///////////////
  // cglclique //
  ///////////////
  options = add_param(options,'cglclique',1);

  /////////////////////////////
  // cglmixedintegerrounding //
  /////////////////////////////
  options = add_param(options,'cglmixedintegerrounding',1);

  //////////////////////////////
  // cglmixedintegerrounding2 //
  //////////////////////////////
  options = add_param(options,'cglmixedintegerrounding2',1);

  //////////////////
  // cglflowcover //
  //////////////////
  options = add_param(options,'cglflowcover',1);

  ////////////////
  // cgloddhole //
  ////////////////
  options = add_param(options,'cgloddhole',1);

  ///////////////
  // cgltwomir //
  ///////////////
  options = add_param(options,'cgltwomir',1);

  /////////////////////
  // cglalldifferent //
  /////////////////////
  options = add_param(options,'cglalldifferent',1);

  /////////////////////
  // cglduplicaterow //
  /////////////////////
  options = add_param(options,'cglduplicaterow',1);

  ///////////////////////
  // cglliftandproject //
  ///////////////////////
  options = add_param(options,'cglliftandproject',1);

  ///////////////////////
  // cglsimplerounding //
  ///////////////////////
  options = add_param(options,'cglresidualcapacity',1);

  /////////////////////////
  // Set some heuristics //
  /////////////////////////

  /////////////////
  // cbcrounding //
  /////////////////
  options = add_param(options,'cbcrounding',1);

  //////////////////////////////////
  // Set Branching decision class //
  //////////////////////////////////

  //////////////////////////////
  // cbcbranchdefaultdecision //
  options = add_param(options,'cbcbranchdefaultdecision_bestcriterion',1);

  ////////////////////////
  // Set Strategy class //
  ////////////////////////

  ////////////////////////
  // cbcstrategydefault //
  ////////////////////////
  options = add_param(options,'cbcstrategydefault',1);

  //////////////////////////////////////
  // Some Cbc parameters are settable //
  //////////////////////////////////////
  options = add_param(options,'cbc_numberstrong', 5);
  options = add_param(options,'cbc_printfrequency', 100);
  options = add_param(options,'cbc_maximumnodes', 10000);
  options = add_param(options,'cbc_maximumseconds', 10000); 
  options = add_param(options,'cbc_dobranchandbound',1);
end

//
// Call Scilab Cbc interface
//

if (typeof(special)=='clist') then
  [xmin,fmin,status,extra] = scicbc(A,c,lhs,rhs,ub,lb,btype,vartype,Q,options,special);
elseif (isempty(special)) then
  special = init_constraint();
  [xmin,fmin,status,extra] = scicbc(A,c,lhs,rhs,ub,lb,btype,vartype,Q,options,special);
else
  error('cbc: the special parameter must be a clist or must be empty');
end
endfunction
