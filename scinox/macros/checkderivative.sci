function checkderivative(fd_deriv, user_deriv)
// checkderivative checks the conistency of the users analytic derivative 

tol = 1e-6; 

if size(fd_deriv) ~= size(user_deriv)
  error('Dimensions of input arguments do not match\n.');
end 

// Since Scilab can't handle max(1.0,A) where A is a sparse matrix 
// we need to check to see wheter the user passed in a sparse or dense matrix 
if ~issparse(fd_deriv) & ~issparse(user_deriv)
  relativeError = abs( user_deriv - fd_deriv )./(max(1.0,abs(user_deriv)));
elseif issparse(fd_deriv) & issparse(user_deriv)
    // the contortions below compute relativeError in the same 
    // manner as the dense case. They are necessary because 
    // max(1.0, abs(user_deriv)) 
    // would be a full matrix even if user_deriv is sparse. 
    // They are also necessary because the above command will 
    // not work in Scilab
    absdiff = abs(fd_deriv - user_deriv);
    pat  = spones(absdiff);
    
    [ij, valud, mn] = spget(user_deriv);
    valud = max(1.0, abs(valud));

    // This is inefficent because we are doing lots of alterations 
    // to the sparse matrix. Luckily there are no inserts. But 
    // I wish knew a better way to do this. 
    for k = 1:size(ij,1)
      if pat(ij(k,1), ij(k,2)) ~= 0 
	pat(ij(k,1),ij(k,2)) = valud(k);
      end
    end
    
    [ij, valdiff, mn] = spget(absdiff);
    [ij, valdiv,  mn] = spget(pat);
    
    relErrorval = valdiff./valdiv; 
    relativeError = sparse(ij,relErrorval,mn);
else
  error('Input arguments must both be dense or sparse');
end
  

// if we call max on a zero sparse matrix this causes a crash! 
// we work around this here.  
if norm(relativeError,%inf) == 0
  maxvalue = 0;
  ij = [1 1];
else
  [maxvalue,ij] = max(relativeError);
end 

// again we need to work around problems with max of a sparse matrix
if length(ij) == 2
  i = ij(1)
  j = ij(2)
elseif length(ij) == 1
  [i,j] = ind2sub(size(relativeError),ij)
else 
  error('Bad index value');
end

printf('Maximum relative discrepancy between derivatives = %g\n',maxvalue);
if or(relativeError > tol)
  printf('Caution: user-supplied and finite-difference derivatives\n');
  printf('do not agree to %g relative tolerance.\n',tol);
  
  printf('Maximum discrepancy occurs in element (%d,%d) of Jacobian:\n',i,j);

  // We need to wrap fulls here to deal with sparse matrices. If A is a sparse
  // matrix then the scalar A(i,j) is also a sparse matrix! Why Scilab? why? 
  printf('\tUser-supplied     derivative:   %g\n',full(user_deriv(i,j)));
  printf('\tFinite-difference derivative:   %g\n',full(fd_deriv(i,j)));
  
  printf('Type resume to continue or abort to stop the solve.\n');
  pause 
end
endfunction 