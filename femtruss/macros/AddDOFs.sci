function U = AddDOFs(U,L)
//
// U = K \ F : FE solution
// L : Table de localisation des noeuds d'appuis
//
// A. Seghir, le 05/07/04

[m,n] = size(L);
L = matrix(L,m*n,1);
L = L(find(L));

L = gsort(-L); L = -L;
n = length(L);
for i = 1:n
  U = AddToVect(U,L(i),0);
end
U = U';
endfunction
