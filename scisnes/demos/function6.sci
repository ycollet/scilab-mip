function [f,J] = function6(x)
   n = size(x,1);
   f = zeros(n,1);

   // the orginal paper doesn't have a x(2) term in f(1) but it 
   // does have a nonzero term in the jacobian at (1,2) 
   // i've decided to add -2*x(2) here 
   f(1) = -2*x(1)^2 + 3*x(1) -2*x(2) + 3*x(n-4) - x(n-3) - x(n-2) + ...
       0.5*x(n-1) - x(n) + 1;
   
   for i=2:n-1
     f(i) = -2*x(i)^2 + 3*x(i) - x(i-1) - 2*x(i+1) + 3*x(n-4) ...
	 -x(n-3) - x(n-2) + 0.5*x(n-1) - x(n) + 1;
   end 

   f(n) = -2*x(n)^2 + 3*x(n) - x(n-1) + 3*x(n-4) - x(n-3) - x(n-2) ...
       + 0.5*x(n-1) - x(n) + 1
   
   [nargout,nargin] = argn();
   // if nargout > 1 supply the sparse Jacobian 
   if nargout > 1 
     // the sparsity structure of the jacobian is 
     // tridiagonal 
     nnzjac1 = 3*(n-2) + 2 + 2; 
     // with 5 dense columns at the end 
     nnzjac2 = 5*n;
     
     // make a jacobian for the tridiagonal part 
     irow = zeros(nnzjac1,1);
     jcol = zeros(nnzjac1,1);
     vals = zeros(nnzjac1,1); 
  
     k = 1; 
  
     for i=1:n
       // \partial f_i/ \partial x_{i-1}
       if (i > 1)
	 irow(k) = i;
	 jcol(k) = i-1;
	 vals(k) = -1;
	 k = k + 1;
       end
       
       // \partial f_i / \partial x_i
       irow(k) = i;
       jcol(k) = i; 
       vals(k) = -4*x(i) + 3;
       k = k+1;
     
       // \partial f_i / \partial x_{i+1} 
       if (i < n)
	 irow(k) = i;
	 jcol(k) = i+1;
	 vals(k) = -2;
	 k = k+1;
       end
    end 
    J1  = sparse([irow jcol],vals,[n,n]);
    
    // now make a jacobian for the dense columns 
    irow = zeros(nnzjac2,1);
    jcol = zeros(nnzjac2,1);
    vals = zeros(nnzjac2,1);
    k = 1;
    
    for i=1:n
      // \partial f_i / \partial x_{n-4}
       irow(k) = i;
       jcol(k) = n-4;
       vals(k) = 3;
       k = k+1;
       
       // \partial f_i / \partial x_{n-3}
       irow(k) = i;
       jcol(k) = n-3;
       vals(k) = -1;
       k = k+1;
       // \partial f_i / \partial x_{n-2}
       irow(k) = i;
       jcol(k) = n-2;
       vals(k) = -1;
       k = k+1;
       
       // \partial f_i / \partial x_{n-1}
       irow(k) = i;
       jcol(k) = n-1;
       vals(k) = 0.5;
       k = k+1;
       
       // \partial f_i / \partial x_{n}
       irow(k) = i;
       jcol(k) = n;
       vals(k) = -1;
       k = k + 1;
    end
    J2 = sparse([irow jcol],vals,[n,n]);
    
    J = J1 + J2;     
   end
endfunction 
 


function [f,J] = function6dense(x)
    [nargout,nargin] = argn();
    if nargout == 1
      f = function6(x);
    else 
      [f,J] = function6(x);
      J = full(J);
    end
endfunction 
 