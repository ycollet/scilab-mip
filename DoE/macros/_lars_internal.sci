function [G,x] = _planerot(x)
// Givens plane rotation.
//   [G,Y] = _PLANEROT(X), where X is a 2-component column vector,
//   returns a 2-by-2 orthogonal matrix G so that Y = G*X has Y(2) = 0.
//
if (x(2) ~= 0) then
   r = norm(x);
   G = [x'; -x(2) x(1)]/r;
   x = [r; 0];
else
   G = eye(2);
end
endfunction

// Fast Cholesky insert and remove functions
// Updates R in a Cholesky factorization R'R = X'X of a data matrix X. R is
// the current R matrix to be updated. x is a column vector representing the
// variable to be added and X is the data matrix containing the currently
// active variables (not including x).
function R = _cholinsert(R, x, X)
diag_k = x'*x; // diagonal element k in X'X matrix
if (isempty(R)) then
  R = sqrt(diag_k);
else
  col_k = x'*X; // elements of column k in X'X matrix
  R_k   = R'\col_k'; // R'R_k = (X'X)_k, solve for R_k
  R_kk  = sqrt(diag_k - R_k'*R_k); // norm(x'x) = norm(R'*R), find last element by exclusion
  R     = [R R_k; [zeros(1,size(R,2)) R_kk]]; // update R
end
endfunction

// Deletes a variable from the X'X matrix in a Cholesky factorisation R'R =
// X'X. Returns the downdated R. This function is just a stripped version of
// Matlab's qrdelete.
function R = _choldelete(R,j)
R(:,j) = []; // remove column j
n      = size(R,2);
for k = j:n
  p = k:k+1;
  [G,R(p,k)] = _planerot(R(p,k)); // remove extra element in column
  if (k < n) then
    R(p,k+1:n) = G*R(p,k+1:n); // adjust rest of row
  end
end
R(size(R,1),:) = []; // remove zero'ed out row
endfunction

