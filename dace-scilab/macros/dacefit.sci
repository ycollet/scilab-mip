function  [dmodel, perf] = dacefit(S, Y, regr, _corr, theta0, lob, upb)
//DACEFIT Constrained non-linear least-squares fit of a given correlation
// model to the provided data set and regression model
//
// Call
//   [dmodel, perf] = dacefit(S, Y, regr, _corr, theta0)
//   [dmodel, perf] = dacefit(S, Y, regr, _corr, theta0, lob, upb)
//
// Input
// S, Y    : Data points (S(i,:), Y(i,:)), i = 1,...,m
// regr    : Function handle to a regression model
// _corr   : Function handle to a correlation function
// theta0  : Initial guess on theta, the correlation function parameters
// lob,upb : If present, then lower and upper bounds on theta
//           Otherwise, theta0 is used for theta
//
// Output
// dmodel  : DACE model: a struct with the elements
//    regr   : function handle to the regression model
//    corr   : function handle to the correlation function
//    theta  : correlation function parameters
//    beta   : generalized least squares estimate
//    gamma  : correlation factors
//    sigma2 : maximum likelihood estimate of the process variance
//    S      : scaled design sites
//    Ssc    : scaling factors for design arguments
//    Ysc    : scaling factors for design ordinates
//    C      : Cholesky factor of correlation matrix
//    Ft     : Decorrelated regression matrix
//    G      : From QR factorization: Ft = Q*G' .
// perf    : struct with performance information. Elements
//    nv     : Number of evaluations of objective function
//    perf   : (q+2)*nv array, where q is the number of elements 
//             in theta, and the columns hold current values of
//                 [theta;  psi(theta);  type]
//             |type| = 1, 2 or 3, indicate 'start', 'explore' or 'move'
//             A negative value for type indicates an uphill step

// hbn@imm.dtu.dk  
// Last update September 3, 2002

[nargout,nargin] = argn();

// Check design points
[m,n] = size(S);  // number of design sites and their dimension
sY    = size(Y);
if min(sY) == 1 then  
  Y  = Y(:);
  lY = max(sY);
  sY = size(Y);
else 
  lY = sY(1);
end
if m ~= lY then
  error('S and Y must have the same number of rows');
end

// Check correlation parameters
lth = length(theta0);
if  nargin > 5 then // optimization case
  if  length(lob) ~= lth | length(upb) ~= lth then
    error('theta0, lob and upb must have the same length');
  end
  if  or(lob <= 0) | or(upb < lob) then
    error('The bounds must satisfy  0 < lob <= upb');
  end
else  // given theta
  if  or(theta0 <= 0) then
    error('theta0 must be strictly positive');
  end
end

// Normalize data
mS = mean(S,'r');   sS = stdev(S,'r');
mY = mean(Y,'r');   sY = stdev(Y,'r');
// 02.08.27: Check for 'missing dimension'
j = find(sS == 0);
if  ~isempty(j) then sS(j) = 1; end
j = find(sY == 0);
if  ~isempty(j) then sY(j) = 1; end

S = (S - (ones(m,1) .*. mS)) ./ (ones(m,1) .*. sS);
Y = (Y - (ones(m,1) .*. mY)) ./ (ones(m,1) .*. sY);

// Calculate distances D between points
mzmax = m*(m-1) / 2;        // number of non-zero distances
ij = zeros(mzmax, 2);       // initialize matrix with indices
D  = zeros(mzmax, n);        // initialize matrix with distances
ll = 0;
for k = 1:m-1
  ll  = ll($) + (1:m-k);
  ij(ll,:) = [(ones(m-k,1) .*. k) (k+1 : m)']; // indices for sparse matrix
  D(ll,:)  = (ones(m-k,1) .*. S(k,:)) - S(k+1:m,:); // differences between points
end

if min(sum(abs(D),2) ) == 0 then
  error('Multiple design sites are not allowed');
end

// Regression matrix
F = regr(S);  [mF,p] = size(F);
if mF ~= m then error('number of rows in  F  and  S  do not match'), end
if p > mF  then error('least squares problem is underdetermined'), end

// parameters for objective function
par = mlist(['dace_par','corr','regr','y','F','D','ij','scS'], ...
             _corr, regr, Y, F, D, ij, sS);
// Determine theta
if  nargin > 5 then
  // Bound constrained non-linear optimization
  [theta,f,fit,perf] = _boxmin(theta0, lob, upb, par);
  if isinf(f) then
    error('Bad parameter region.  Try increasing  upb');
  end
else
  // Given theta
  theta = theta0(:);   
  [f,fit] = objfunc(theta, par);
  perf = mlist(['dace_perf','perf','nv'],[theta;f;1],1);
  if isinf(f) then
    error('Bad point. Try increasing theta0');
  end
end

// Return values
dmodel = mlist(['dace_model','regr','corr','theta','beta','gamma','sigma2','S','Ssc','Ysc','C','Ft','G'], ...
                             regr, _corr, theta', fit('beta'),fit('gamma'), ...
                             sY.^2.*fit('sigma2'), S, [mS;sS],[mY;sY],fit('C'),fit('Ft'),fit('G'));
endfunction

// >>>>>>>>>>>>>>>>   Auxiliary functions  ====================

function  [obj, fit] = objfunc(theta, par)
// Initialize
obj = %inf; 
fit = mlist(['dace_fit','sigma2','beta','gamma','C','Ft','G'],%nan,%nan,%nan,%nan,%nan,%nan);
m   = size(par('F'),1);
// Set up  R
clear tmp; tmp = par('corr');
r   = tmp(theta, par('D'));
o   = (1:m)';
mu  = (10+m)*%eps;
idx = find(r>0);   
R   = sparse([[par('ij')(idx,1); o] [par('ij')(idx,2); o]], [r(idx); ones(m,1)+mu]);
// Cholesky factorization with check for pos. def.
//[C,rd] = spchol(R);
//if rd then return, end // not positive definite
[C] = chol(full(R));

// Get least squares solution
C = C';   Ft = C \ par('F');
//[Q G] = qr(Ft,0);
[Q,G] = qr(Ft,'e');
if  rcond(G) < 1e-10 then
  // Check   F  
  if  cond(par('F')) > 1e15 then
    T = sprintf('F is too ill conditioned\nPoor combination of regression model and design sites');
    error(T);
  else  // Matrix  Ft  is too ill conditioned
    return 
  end 
end
Yt   = C \ par('y'); _beta = G \ (Q'*Yt);
rho  = Yt - Ft*_beta;  sigma2 = sum(rho.^2)/m;
detR = prod(diag(C) .^ (2/m));
obj  = sum(sigma2) * detR;
if nargout > 1 then
  fit = mlist(['dace_fit','sigma2','beta','gamma','C','Ft','G'],sigma2,_beta,rho'/C,C,Ft,G');
end
endfunction

// --------------------------------------------------------

function  [t, f, fit, perf] = _boxmin(t0, lo, up, par)
//BOXMIN  Minimize with positive box constraints

// Initialize
[t,f,fit,itpar] = _start(t0, lo, up, par);

if  ~isinf(f) then
  // Iterate
  p = length(t);
  if  p <= 2 then  kmax = 2; else,  kmax = min(p,4); end
  for  k = 1 : kmax
    th = t;
    [t, f, fit, itpar] = _explore(t, f, fit, itpar, par);
    [t, f, fit, itpar] = _move(th, t, f, fit, itpar, par);
  end
end
perf = mlist(['dace_perf','nv','perf'],itpar('nv'),itpar('perf')(:,1:itpar('nv')));
endfunction

// --------------------------------------------------------

function  [t,f,fit,itpar] = _start(t0, lo, up, par)
// Get starting point and iteration parameters

// Initialize
t = t0(:);  lo = lo(:);   up = up(:);   p = length(t);
D = 2 .^ ([1:p]'/(p+2));
ee = find(up == lo);  // Equality constraints
if  ~isempty(ee) then
  D(ee) = ones(length(ee),1);   t(ee) = up(ee); 
end
ng = find(t < lo | up < t);  // Free starting values
if  ~isempty(ng) then
  t(ng) = (lo(ng) .* up(ng).^7).^(1/8);  // Starting point
end
ne = find(D ~= 1);

// Check starting point and initialize performance info
[f,fit] = objfunc(t,par); nv = 1;
itpar = mlist(['dace_par','D','ne','lo','up','perf','nv'],D,ne,lo,up,zeros(p+2,200*p),1);
itpar('perf')(:,1) = [t; f; 1];
if isinf(f) then   // Bad parameter region
  return
end

if  length(ng) > 1 then // Try to improve starting guess
  d0 = 16;  d1 = 2;   q = length(ng);
  th = t;   fh = f;   jdom = ng(1);  
  for  k = 1 : q
    j = ng(k);    fk = fh;  tk = th;
    DD = ones(p,1);  DD(ng) = ones(q,1) .*. (1/d1);  DD(j) = 1/d0;
    alpha = min(log(lo(ng) ./ th(ng)) ./ log(DD(ng))) / 5;
    v = DD .^ alpha;   tk = th;
    for  rept = 1 : 4
      tt = tk .* v; 
      [ff  fitt] = objfunc(tt,par);  nv = nv+1;
      itpar('perf')(:,nv) = [tt; ff; 1];
      if  ff <= fk then
        tk = tt;  fk = ff;
        if  ff <= f then
          t = tt;  f = ff;  fit = fitt; jdom = j;
        end
      else
        itpar('perf')($,nv) = -1;
        break;
      end
    end
  end // improve
  
  // Update Delta  
  if  jdom > 1 then
    D([1 jdom]) = D([jdom 1]); 
    itpar('D') = D;
  end
end // free variables

itpar('nv') = nv;
endfunction

// --------------------------------------------------------

function  [t, f, fit, itpar] = _explore(t, f, fit, itpar, par)
// Explore step

nv = itpar('nv');   ne = itpar('ne');
for  k = 1 : length(ne)
  j = ne(k);   tt = t;   DD = itpar('D')(j);
  if  t(j) == itpar('up')(j) then
    atbd = 1;   tt(j) = t(j) / sqrt(DD);
  elseif  t(j) == itpar('lo')(j) then
    atbd = 1;  tt(j) = t(j) * sqrt(DD);
  else
    atbd = 0;  tt(j) = min(itpar('up')(j), t(j)*DD);
  end
  [ff,fitt] = objfunc(tt,par);  nv = nv+1;
  itpar('perf')(:,nv) = [tt; ff; 2];
  if  ff < f then
    t = tt;  f = ff;  fit = fitt;
  else
    itpar('perf')($,nv) = -2;
    if ~atbd then // try decrease
      tt(j) = max(itpar('lo')(j), t(j)/DD);
      [ff  fitt] = objfunc(tt,par);  nv = nv+1;
      itpar('perf')(:,nv) = [tt; ff; 2];
      if  ff < f then
        t = tt;  f = ff;  fit = fitt;
      else
        itpar('perf')($,nv) = -2;
      end
    end
  end
end // k

itpar('nv') = nv;
endfunction

// --------------------------------------------------------

function  [t, f, fit, itpar] = _move(th, t, f, fit, itpar, par)
// Pattern move

nv = itpar('nv'); ne = itpar('ne'); p = length(t);
v = t ./ th;
if  and(v == 1) then
  itpar('D') = itpar('D')([2:p 1]).^.2;
  return
end

// Proper move
rept = 1;
while  rept
  tt = min(itpar('up'), max(itpar('lo'), t .* v));  
  [ff,fitt] = objfunc(tt,par);  nv = nv+1;
  itpar('perf')(:,nv) = [tt; ff; 3];
  if  ff < f then
    t = tt;  f = ff;  fit = fitt;
    v = v .^ 2;
  else
    itpar('perf')($,nv) = -3;
    rept = 0;
  end
  if  or(tt == itpar('lo') | tt == itpar('up')) then rept = 0; end
end

itpar('nv') = nv;
itpar('D') = itpar('D')([2:p 1]).^.25;
endfunction

