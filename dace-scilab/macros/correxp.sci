function  [r, dr] = correxp(theta, d)
//CORREXP  Exponential correlation function
//
//           n
//   r_i = prod exp(-theta_j * |d_ij|)
//          j=1
//
// If length(theta) = 1, then the model is isotropic: 
// theta_j = theta(1), j=1,...,n
//
// Call:    r = correxp(theta, d)
//          [r, dr] = correxp(theta, d)
//
// theta :  parameters in the correlation function
// d     :  m*n matrix with differences between given data points
// r     :  correlation
// dr    :  m*n matrix with the Jacobian of r at x. It is
//          assumed that x is given implicitly by d(i,:) = x - S(i,:), 
//          where S(i,:) is the i'th design site. 

// hbn@imm.dtu.dk  
// Last update April 12, 2002

[nargout, nargin] = argn();

[m,n] = size(d);  // number of differences and dimension of data
lt = length(theta);
if  lt == 1 then
  theta = ones(1,n) .*. theta;
elseif  lt ~= n then
  error(sprintf('Length of theta must be 1 or %d',n))
else
  theta = theta(:).';
end

td = abs(d) .* (ones(m,1) .*. (-theta));
r  = exp(sum(td,'c'));

if nargout > 1 then
  dr = (ones(m,1) .*. (-theta)) .* sign(d) .* (ones(1,n) .*. r);
end
endfunction
