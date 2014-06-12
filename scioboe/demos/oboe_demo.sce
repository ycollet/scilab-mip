// Initialization of some oboe parameters
params = init_param();
[params, err] = add_param(params, 'NumSubProblems', 1); // Be careful:
//if this parameter is not correctly set, oboe may hangs
[params, err] = add_param(params, 'ProblemName', 'LATest');
[params, err] = add_param(params, 'MaxOuterIterations', 50);
[params, err] = add_param(params, 'MaxInnerIterations', 10);
// [params, err] = add_param(params, 'ObjectiveLB', -0.416666666666666);
[params, err] = add_param(params, 'Verbosity', 3);
[params, err] = add_param(params, 'Proximal', 1);
[params, err] = add_param(params, 'Ball', 1);
[params, err] = add_param(params, 'RadiusBall', 5);
[params, err] = add_param(params, 'Tolerance', 0.00001);
[params, err] = add_param(params, 'WeightEpigraphCutInit', 1.0);
[params, err] = add_param(params, 'WeightEpigraphCutInc', 0.0);

// Definition of the objective function
// The prototype must be: [fobj_val, subgrad_val,info] = fobj(y)
function [fobj_val, subgrad_val, info] = fobj_oboe(y_in)
  global _A;
  global _b;
  fobj_val = 0.5*y_in'*_A*y_in - _b'*y_in;
  subgrad_val = _A*y_in - _b;
   
  info = 1;
endfunction

// Initialization parameters
addEqC      = %T;
n           = 5;
x_start     = zeros(n,1);
x_lower     = -100*ones(x_start);
x_upper     =  100*ones(x_start);
center_ball = 0.5*ones(n,1);

// Initialization of equality constraints
constraint = [];
rhs        = [];
if (addEqC) then
  constraint      = zeros(3,n);
  constraint(1,2) = 1;
  constraint(1,4) = 1;
  constraint(2,3) = 1;
  constraint(3,1) = 1;
  constraint(3,5) = 1;
  rhs = zeros(3,1);
  rhs(1) = 1;
  rhs(2) = 0.5;
  rhs(3) = 1;
end

// Initialization of problem matrixes
global _A;
global _b;
_A = zeros(n,n);
_b = zeros(n,1);

// Initialize A as the symmetric difference matrix and b as an identity vector
_b(1)     = 1;
_A(1,1)   = 2; 
_A(1,2)   = -1;
_A(n,n)   = 2;
_A(n,n-1) = -1;

for(i=2:n-1)
  _A(i,i-1) = -1;
  _A(i,i)   = 2;
  _A(i,i+1) = -1;
end

// Start optimization
[x_out, status_out] = oboe(x_start,fobj_oboe,x_lower,x_upper,constraint,rhs,center_ball,params);

printf('The solution found is:'); disp(x_out);

printf('Status code are:\n\n');
printf('-5: LOCSET_EMPTY\n');
printf('-4: CONVEXITY_FAILURE\n');
printf('-3: LA_ERROR\n');
printf('-2: CHOLESKY_FAILURE\n');
printf('-1: UNKNOWN\n');
printf(' 0: ITERATING\n');
printf(' 2: RELATIVE_GAP_REACHED\n');
printf(' 3: USER_STOP\n');
printf(' 4: MAX_OUTER_ITERATIONS\n\n');

printf('The status of the solver is %d\n', status_out);
