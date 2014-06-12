function J=sfdjac(userfun,x,Jpattern,coloring)
// J = sfdjac(userfun,x,Jpattern,coloring) Computes a finite-difference 
// approximation of a sparse Jacobian 
//
// Input:  userfun      A user-defined function of the form 
//                      f = userfun(x) 
//                      which returns f(x), a vector of length m
//
//         x            A vector of length n, at which to evaluate J(x) 
//
//         Jpattern     A m x n matrix that contains the sparsity pattern 
//                      of the Jacobian.  Jpattern = spones(J(x))
//
//         coloring     A coloring of the columns of the Jacobian such 
//                      that coloring(j) = k implies that the jth column 
//                      of J is colored k.  
//                       
//                      A coloring is a partition of the columns of 
//                      the Jacobian into distinct groups (or colors) such 
//                      that no column of the same color has a nonzero entry 
//                      in the same row.
// 
//                      All columns of the same color can be estimated via 
//                      a single finite-difference. Thus, the number of
//                      of function evaluations required is equal to the number
//                      of colors plus one.
//
// Output:  J           A m x n sparse matrix containing the finite-difference
//                      approximation of the Jacobian 

epsilon = sqrt(%eps);

[m,n] = size(Jpattern);

if length(coloring) ~= n 
    error('Must specify a color for each column');
end 

numcolors = max(coloring);
if numcolors > n 
  fprintf('Poor coloring, contains more than %d colors but only %d required\n',n,n);
  // If there are more than n colors use a dense coloring 
  coloring = 1:n;
end


f = userfun(x);
if size(f) ~= [m, 1]
  error('Inconsistent size of f and Jpattern\n');
end

// This ensures that J has the sparsity structure as Jpattern. 
J = spones(Jpattern);

// we use the english spelling of color to avoid clashes with 
// the built-in function color 

for colour=1:numcolors
  // p is a binary vector corresponding to those columns with a
  // particular color
  p = (coloring == colour);
  j = find(p);
   
  fp = userfun(x+epsilon*p);
  
  // w is a dense m x 1 vector 
  // w =~ J*p
  w = (fp-f)/epsilon;

  // This is probably inefficent, but there doesn't seem to 
  // be a way to adjust just the nonzero pattern of J
  for l=1:length(j)
      i = find(J(:,j(l)));
      for k=1:length(i)
          J(i(k),j(l)) = w(i(k));
      end
  end 

end


endfunction 