function X = normalize(X, replace)
// NORMALIZE  Normalize the observations of a data matrix.
//    X = NORMALIZE(X) centers and scales the observations of a data
//    matrix such that each variable (column) has unit length.
//
// Author: Karl Skoglund, IMM, DTU, kas@imm.dtu.dk
[nargout, nargin] = argn();

if (nargin==1) then
  replace = %F;
end
[n p] = size(X);
X = center(X);
for i=1:p
  B = sum(X(:,1).^2);
  if (replace) then
    if (B<=%eps) then
      B = 1;
    end
  end
  X(:,i) = X(:,i)./sqrt(ones(n,1)*sum(X(:,i).^2));
end
endfunction
