function R = AddDOFsToVect(U,L,V)
// R = AddDOFsToVect(U,L,V)
// U : input vector
// R : resulting vector
// L : inserting position
// V : value to add (optional). If not present: add 0
// A. Seghir, le 01/08/04

if ~isdef('V','local') then
  V = zeros(length(L),1);
end

n = length(L);
L = matrix(L,n,1);
L = L(find(L));

L = gsort(-L); L = -L;
n = length(L);

R = U;
for i=1:n
  tmp     = R(L(i):$);
  R(L(i)) = V(i);
  R(L(i)+1:$+1) = tmp;
end
endfunction
