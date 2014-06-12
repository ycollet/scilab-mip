function [f,J] = trigexp(x)
// Trigexp function 
n = size(x,1);
f = zeros(n,1);

f(1) = 3*x(1)^3 + 2*x(2) - 5 + sin(x(1) - x(2))*sin(x(1) + x(2));

for i=2:n-1 
  f(i) = -x(i-1)*exp(x(i-1)-x(i)) + x(i)*(4+3*x(i)^2) + 2*x(i+1) ...
      + sin(x(i) - x(i+1))*sin(x(i) + x(i+1)) - 8;
end

f(n) = -x(n-1)*exp(x(n-1) - x(n)) + 4*x(n) - 3;

[nargout,nargin] = argn();
// evaluate the Jacobian if nargout > 1 
if nargout > 1 
  // The Jacobian is tridiagonal 
  // with the following number of nonzeros 
  nnzjac = 2 + 3*(n-2) + 2;
  
  irow = zeros(nnzjac,1);
  jcol = zeros(nnzjac,1);
  vals = zeros(nnzjac,1);
  
  k = 1
  
  // \partial f1/ \partial x1 
  irow(k) = 1;
  jcol(k) = 1;
  vals(k) = 9*x(1)^2 + cos(x(1) + x(2))*sin(x(1)-x(2)) ...
      + cos(x(1) - x(2))*sin(x(1) + x(2));
  k = k+1;
  
  // \partial f1 / \partial x2 
  irow(k) = 1;
  jcol(k) = 2;
  vals(k) = 2 + cos(x(1) + x(2))*sin(x(1) - x(2)) - cos(x(1) - x(2))*sin(x(1)+x(2));
  k = k+1;
  
  for i=2:n-1
    // \partial fi / \partial x_{i-1}
    irow(k) = i;
    jcol(k) = i-1;
    vals(k) = -exp(-x(i)+x(i-1)) - exp(-x(i) + x(i-1))*x(i-1)
    k = k+1
    
    // \partial fi / \partial xi 
    irow(k) = i;
    jcol(k) = i;
    vals(k) = 4 + 9*x(i)^2 + exp(-x(i) + x(i-1))*x(i-1) +  ...
	cos(x(i) + x(i+1))*sin(x(i) - x(i+1)) + ...
	cos(x(i) - x(i-1))*sin(x(i) + x(i+1));
    k = k+1;
    
    // \partial fi / \partial x_{i+1}
    irow(k) = i;
    jcol(k) = i+1;
    vals(k) = 2 + cos(x(i) + x(i+1))*sin(x(i)-x(i+1)) - ...
	cos(x(i) - x(i+1))*sin(x(i)+x(i+1));
    k=k+1;
  end
  
  // \partial fn/ \partial x_{n-1}
  irow(k) = n;
  jcol(k) = n-1;
  vals(k) = -exp(-x(n) + x(n-1)) - exp(-x(n)+x(n-1))*x(n-1);
  k = k + 1;
  
  // \partial fn/ \partial x_n 
  irow(k) = n;
  jcol(k) = n;
  vals(k) = 4 + exp(-x(n) + x(n-1))*x(n-1);
  
  // create a scilab sparse Jacobian 
  J = sparse([irow jcol], vals, [n,n]);
end 

endfunction 


function [f,J] = trigexpdense(x) 
   [nargout,nargin] = argn();
   
   if nargout == 1 
     f = trigexp(x);
   end
  
   if nargout > 1 
     [f,J] = trigexp(x);
     J = full(J);
   end
endfunction 