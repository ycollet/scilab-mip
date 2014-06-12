function _beta = lars(X, y, method, stop, useGram, Gram, _trace)
// LARS  The LARS algorithm for performing LAR or LASSO.
//    _beta = LARS(X, Y) performs least angle regression on the variables in
//    X to approximate the response Y. Variables X are assumed to be
//    normalized (zero mean, unit length), the response Y is assumed to be
//    centered.
//    _beta = LARS(X, Y, METHOD), where METHOD is either 'LARS' or 'LARS'  
//    determines whether least angle regression or lasso regression should
//    be performed. 
//    _beta = LARS(X, Y, METHOD, STOP) with nonzero STOP will perform least
//    angle or lasso regression with early stopping. If STOP is negative,
//    STOP is an integer that determines the desired number of variables. If
//    STOP is positive, it corresponds to an upper bound on the L1-norm of
//    the _beta coefficients.
//    _beta = LARS(X, Y, METHOD, STOP, USEGRAM) specifies whether the Gram
//    matrix X'X should be calculated (USEGRAM = 1) or not (USEGRAM = 0).
//    Calculation of the Gram matrix is suitable for low-dimensional
//    problems. By default, the Gram matrix is calculated.
//    _beta = LARS(X, Y, METHOD, STOP, USEGRAM, GRAM) makes it possible to
//    supply a pre-computed Gram matrix. Set USEGRAM to 1 to enable. If no
//    Gram matrix is available, exclude argument or set GRAM = [].
//    _beta = LARS(X, Y, METHOD, STOP, USEGRAM, GRAM, _trace) with nonzero
//    _trace will print the adding and subtracting of variables as all
//    LARS/lasso solutions are found.
//    Returns _beta where each row contains the predictor coefficients of
//    one iteration. A suitable row is chosen using e.g. cross-validation,
//    possibly including interpolation to achieve sub-iteration accuracy.
//
// Author: Karl Skoglund, IMM, DTU, kas@imm.dtu.dk
// Reference: 'Least Angle Regression' by Bradley Efron et al, 2003.

// Input checking
// Set default values.

[nargout, nargin] = argn();

if (nargin < 7) then
    _trace = 0;
end
if (nargin < 6) then
  Gram = [];
end
if (nargin < 5) then
  useGram = 1;
end
if (nargin < 4) then
  stop = 0;
end
if (nargin < 3) then
  method = 'lars';
end
if (method=='lasso') then
  lasso = 1;
else
  lasso = 0;
end

// LARS variable setup
[n p] = size(X);
nvars = min(n-1,p); 
maxk = 8*nvars; // Maximum number of iterations

if (stop == 0) then
  _beta = zeros(2*nvars, p);
elseif (stop < 0) then
  _beta = zeros(2*round(-stop), p);
else
  _beta = zeros(100, p);
end
mu = zeros(n, 1); // current "position" as LARS travels towards lsq solution
I = 1:p; // inactive set
A = [];  // active set

// Calculate Gram matrix if necessary
if (isempty(Gram) & useGram) then
  Gram = X'*X; // Precomputation of the Gram matrix. Fast but memory consuming.
end
if (~useGram) then
  R = []; // Cholesky factorization R'R = X'X where R is upper triangular
end

lassocond = 0; // LASSO condition boolean
stopcond  = 0; // Early stopping condition boolean
k    = 0; // Iteration count
vars = 0; // Current number of variables

if (_trace) then
  printf('Step\t\tAdded\t\tDropped\t\tActive set size\n');
end

// LARS main loop
while (vars < nvars) & (~stopcond) & (k < maxk)
  k = k + 1;
  c = X'*(y - mu);
  [C j] = max(abs(c(I)));
  j     = I(j);

  if (~lassocond) then // if a variable has been dropped, do one iteration with this configuration (don't add new one right away)
    if (~useGram) then
      R = _cholinsert(R,X(:,j),X(:,A));
    end
    A         = [A j];
    I(I == j) = [];
    vars      = vars + 1;
    if (_trace) then
      printf('%d\t\t%d\t\t\t\t\t%d\n', k, j, vars);
    end
  end

  s = sign(c(A)); // get the signs of the correlations

  if (useGram) then
    S   = s*ones(1,vars);
    GA1 = inv(Gram(A,A).*S'.*S)*ones(vars,1);
    AA  = 1/sqrt(sum(GA1));
    w   = AA*GA1.*s; // weights applied to each active variable to get equiangular direction
  else
    GA1 = R\(R'\s);
    AA  = 1/sqrt(sum(GA1.*s));
    w   = AA*GA1;
  end
  u = X(:,A)*w; // equiangular direction (unit vector)
  
  if (vars == nvars) then // if all variables active, go all the way to the lsq solution
    _gamma = C/AA;
  else
    a      = X'*u; // correlation between each variable and eqiangular vector
    temp   = [(C - c(I))./(AA - a(I)); (C + c(I))./(AA + a(I))];
    _gamma = min([temp(temp > 0); C/AA]);
  end

  // LASSO modification
  if (lasso) then
    lassocond = 0;
    temp      = -_beta(k,A)./w';
    [_gamma_tilde] = min([temp(temp > 0) _gamma]);
    j = find(temp == _gamma_tilde);
    if (_gamma_tilde < _gamma) then
      _gamma    = _gamma_tilde;
      lassocond = 1;
    end
  end

  mu = mu + _gamma*u;
  if (size(_beta,1) < k+1) then
    _beta = [_beta; zeros(size(_beta,1), p)];
  end
  _beta(k+1,A) = _beta(k,A) + _gamma*w';

  // Early stopping at specified bound on L1 norm of _beta
  if (stop > 0) then
    t2 = sum(abs(_beta(k+1,:)));
    if (t2 >= stop) then
      t1 = sum(abs(_beta(k,:)));
      s = (stop - t1)/(t2 - t1); // interpolation factor 0 < s < 1
      _beta(k+1,:) = _beta(k,:) + s*(_beta(k+1,:) - _beta(k,:));
      stopcond = 1;
    end
  end
  
  // If LASSO condition satisfied, drop variable from active set
  if (lassocond == 1) then
    if (~useGram) then
      R = _choldelete(R,j);
    end
    I = [I A(j)];
    A(j) = [];
    vars = vars - 1;
    if (_trace) then
      printf('%d\t\t\t\t%d\t\t\t%d\n', k, j, vars);
    end
  end
  
  // Early stopping at specified number of variables
  if (stop < 0) then
    stopcond = (vars >= -stop);
  end
end

// trim _beta
if (size(_beta,1) > k+1) then
  _beta(k+2:size(_beta,1), :) = [];
end

if (k == maxk) then
  printf('LARS warning: Forced exit. Maximum number of iteration reached.\n');
end
endfunction

