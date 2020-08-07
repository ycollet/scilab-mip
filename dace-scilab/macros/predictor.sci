function  [y, or1, or2, dmse] = predictor(x, dmodel)
//PREDICTOR  Predictor for y(x) using the given DACE model.
//
// Call:   y = predictor(x, dmodel)
//         [y, or] = predictor(x, dmodel)
//         [y, dy, mse] = predictor(x, dmodel) 
//         [y, dy, mse, dmse] = predictor(x, dmodel) 
//
// Input
// x      : trial design sites with n dimensions.  
//          For mx trial sites x:
//          If mx = 1, then both a row and a column vector is accepted,
//          otherwise, x must be an mx*n matrix with the sites stored
//          rowwise.
// dmodel : Struct with DACE model; see DACEFIT
//
// Output
// y    : predicted response at x.
// or   : If mx = 1, then or = gradient vector/Jacobian matrix of predictor
//        otherwise, or is an vector with mx rows containing the estimated
//                   mean squared error of the predictor
// Three or four results are allowed only when mx = 1,
// dy   : Gradient of predictor; column vector with  n elements
// mse  : Estimated mean squared error of the predictor;
// dmse : Gradient vector/Jacobian matrix of mse

[nargout,nargin] = argn();

// hbn@imm.dtu.dk
// Last update August 26, 2002
 
or1 = %nan;   or2 = %nan;  dmse = %nan;  // Default return values
if isnan(dmodel('beta')) then
  y = %nan;   
  error('DMODEL has not been found')
end

[m,n] = size(dmodel('S'));      // number of design sites and number of dimensions
sx = size(x);                // number of trial sites and their dimension
if min(sx) == 1 & n > 1 then // Single trial point 
  nx = max(sx);
  if nx == n then
    mx = 1; x = x(:)';
  end
else
  mx = sx(1); nx = sx(2);
end
if nx ~= n then
  error(sprintf('Dimension of trial sites should be %d',n))
end

// Normalize trial sites  
x = (x - (ones(mx,1) .*. dmodel('Ssc')(1,:))) ./ (ones(mx,1) .*. dmodel('Ssc')(2,:));
q = size(dmodel('Ysc'),2);  // number of response functions
y = zeros(mx,q);         // initialize result

f = []; df = []; r = []; dr = [];

if mx == 1 then // one site only
  dx = ones(m,1) .*. x - dmodel('S');  // distances to design sites
  if  nargout > 1 then                 // gradient/Jacobian wanted
    clear tmp; tmp = dmodel('regr');
    [f,df] = tmp(x);
    clear tmp; tmp = dmodel('corr');
    [r,dr] = tmp(dmodel('theta'), dx);
    // Scaled Jacobian
    dy = (df * dmodel('beta'))' + dmodel('gamma') * dr;
    // Unscaled Jacobian
    or1 = (dy .* (ones(1,nx) .*. dmodel('Ysc')(2, :)')) ./ (ones(q,1) .*. dmodel('Ssc')(2,:));
    if q == 1 then
      // Gradient as a column vector
      or1 = or1';
    end
    if  nargout > 2 then // MSE wanted
      rt = dmodel('C') \ r;
      u = dmodel('Ft')' * rt - f';
      v = dmodel('G') \ u;
      or2 = (ones(mx,1) .*. dmodel('sigma2')) .* (ones(1,q) .*. (1 + sum(v.^2) - sum(rt.^2))');
      
      if  nargout > 3 then // gradient/Jacobian of MSE wanted
        // Scaled gradient as a row vector
        Gv = dmodel('G')' \ v;
        g = (dmodel('Ft') * Gv - rt)' * (dmodel('C') \ dr) - (df * Gv)';
        // Unscaled Jacobian
        dmse = (ones(1,nx) .*. (2 * dmodel('sigma2')')) .* (ones(q,1) .*. (g ./ dmodel('Ssc')(2,:)));
        if q == 1 then
          // Gradient as a column vector
          dmse = dmse';
        end
      end
    end
  else  // predictor only
    clear tmp; tmp = dmodel('regr');
    f = tmp(x);
    clear tmp; tmp = dmodel('corr');
    r = tmp(dmodel('theta'), dx);
  end

  // Scaled predictor
  sy = f * dmodel('beta') + (dmodel('gamma')*r)';
  // Predictor
  y = (dmodel('Ysc')(1,:) + dmodel('Ysc')(2,:) .* sy)';
else  // several trial sites
  // Get distances to design sites  
  dx = zeros(mx*m,n);  kk = 1:m;
  for  k = 1 : mx
    dx(kk,:) = (ones(m,1) .*. x(k,:)) - dmodel('S');
    kk = kk + m;
  end
  // Get regression function and correlation
  clear tmp; tmp = dmodel('regr');
  f = tmp(x);
  clear tmp; tmp = dmodel('corr');
  r = tmp(dmodel('theta'), dx);
  r = matrix(r, m, mx);
  
  // Scaled predictor 
  sy = f * dmodel('beta') + (dmodel('gamma') * r)';
  // Predictor
  y = (ones(mx,1) .*. dmodel('Ysc')(1,:)) + (ones(mx,1) .*. dmodel('Ysc')(2,:)) .* sy;
  
  if  nargout > 1 then  // MSE wanted
    rt  = dmodel('C') \ r;
    u   = dmodel('G') \ (dmodel('Ft')' * rt - f');
    or1 = (ones(mx,1) .*. dmodel('sigma2')) .* (ones(1,q) .*. (1 + colsum(u.^2) - colsum(rt.^2))');
    if nargout > 2 then
      disp('WARNING from PREDICTOR.  Only  y  and  or1=mse  are computed')
    end
  end  
end // of several sites
endfunction

// >>>>>>>>>>>>>>>>   Auxiliary function  ====================

function  s = colsum(x)
// Columnwise sum of elements in  x
if size(x,1) == 1 then
  s = x;
else 
  s = sum(x,'r');
end
endfunction

