function  [r, dr] = corrgauss(theta, d)
//CORRGAUSS  Gaussian correlation function,
//
//           n
//   r_i = prod exp(-theta_j * d_ij^2) ,  i = 1,...,m
//          j=1
//
// If length(theta) = 1, then the model is isotropic:
// all  theta_j = theta .
//
// Call:    r = corrgauss(theta, d)
//          [r, dr] = corrgauss(theta, d)
//
// theta :  parameters in the correlation function
// d     :  m*n matrix with differences between given data points
// r     :  correlation
// dr    :  m*n matrix with the Jacobian of r at x. It is
//          assumed that x is given implicitly by d(i,:) = x - S(i,:), 
//          where S(i,:) is the i'th design site. 

// hbn@imm.dtu.dk  
// Last update June 2, 2002

[nargout,nargin] = argn();

[m n] = size(d);  // number of differences and dimension of data
if  length(theta) == 1 then
  theta = ones(1,n) .*. theta;
elseif  length(theta) ~= n then
  error(sprintf('Length of theta must be 1 or %d',n))
end

td = d.^2 .* (ones(m,1) .*. (-theta(:).'));
r = exp(sum(td,'c'));

if  nargout > 1 then
  dr = (ones(m,1) .*. (-2*theta(:).')) .* d .* (ones(1,n) .*. r);
end
endfunction

