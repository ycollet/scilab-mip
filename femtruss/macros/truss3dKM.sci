function [K,M] = truss3dKM(t,p,A,E,rho)
// [K,M] = truss2dKM(t,p,A, rho)
// 
// K  : matrice de rigidite assemblee 
// M  : matrice masse assemblee 
// t  : table de connectivites des elements
// p  : table des coordonees des noeuds
// A  : sections des elements (des barres)
// E  : modules d'elasticite en vecteur des elements 
// rho: masse volumique
//
//A. Seghir, le 07/08/04, modifie : 26/10/04

net = size(t,1);
nnt = 3*size(p,1);

// YC  
// K = sparse(nnt,nnt);
// M = sparse(nnt,nnt);
K = zeros(nnt,nnt);
M = zeros(nnt,nnt);
for i = 1:net
  ti = t(i,:);
  Li = localise3d(ti);
  Ke = truss3dKe(p(ti,:),A(i),E(i));
  Me = truss3dMe(p(ti,:),A(i),rho(i));
  K(Li,Li) = K(Li,Li) + Ke;
  M(Li,Li) = M(Li,Li) + Me;
end
endfunction
