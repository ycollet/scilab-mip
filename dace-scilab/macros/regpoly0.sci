function  [f, df] = regpoly0(S)
//REGPOLY0  Zero order polynomial regression function
//
// Call:    f = regpoly0(S)
//          [f, df] = regpoly0(S)
//
// S  : m*n matrix with design sites
// f  : ones(m,1)
// df : Jacobian at the first point (first row in S) 

// hbn@imm.dtu.dk  
// Last update  April 12, 2002

[nargout,nargin] = argn();

[m,n] = size(S);
f = ones(m,1);
if nargout > 1 then
  df = zeros(n,1);
end
endfunction

