function A = DelDOFs(A,L)
// A = DelDOFs(A,L)
//
// A : matrice globale apres assemblage
// L : table de localisation des noeuds d'appuis (les DDL a enlever)
//
// A. Seghir, le 01/08/04 modifie le 08/08/04

[m,n] = size(L);

L = matrix(L,m*n,1);
L = L(find(L)); // on supprime les indices nuls
L = gsort(-L); L = unique(-L); // on classe de facon unique les noeuds

n = length(L);

if (min(size(A)) == 1) then
  for i = n:-1:1
    A(L(i)) =[];
  end
else
  for i = n:-1:1
    A(L(i),:) =[];
    A(:,L(i)) =[];
  end
end
endfunction
