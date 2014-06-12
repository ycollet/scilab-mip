function J=fdjac(userfun,x) 
// J = fdjac(userfun,x) Computes a finite-difference approximation of the Jacobian
//
// Note that fdjac computes a dense Jacobian, as such it requires n+1 function 
// evaluations. If the Jacobian is sparse, and the sparsity pattern is known 
// a finite-difference approximation of the Jacobian can be estimated in 
// far fewer function evaulations
//
// Input:  userfun      A user-defined function of the form 
//                      f = userfun(x)  
//                      which returns f(x), a vector of length m. 
//         x            A vector of length n, at which to evaluate J(x).
// Output: J            A m x n matrix containing the finite-difference 
//                      approximation of the Jacobian.

[n, l] = size(x); 
if l ~= 1
  error("Input x should be a n x 1 vector.\n");
end 

f = userfun(x);

[m, l] = size(f);
if l ~= 1 
  error("userfun should return a m x 1 vector.\n");
end

J = zeros(m,n);

epsilon = sqrt(%eps);
for j=1:n
  // We use the simple forward difference formula
  // J*e_j =~ (f(x + epsilon*e_j) - f(x))/epsilon
  // to compute each column of the Jacobian
  ej = zeros(n,1);
  ej(j) = 1;
  
  fp = userfun(x + epsilon*ej);
  J(:,j) = (fp - f)/epsilon;
end 

endfunction 