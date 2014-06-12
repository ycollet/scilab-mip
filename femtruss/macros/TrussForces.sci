function [F,R] = TrussForces(t,p,A,E,U)
// Forces = TrusAxialForces(t,p,A,E,U)
// t : table de connectivites des elements
// p : table des coordonees des noeuds
// A : sections des elements
// E : modules d'elasticite 
// U : solution en deplacelment
//
// A. Seghir, le 07/08/04 modifie : 27/08/04

net = size(t,1);
nnt = size(p,1);

_3D_problem = (size(p,2)==3);

if _3D_problem then R = zeros(3*nnt,1);
else                R = zeros(2*nnt,1);
end

for ie = 1:net
  if _3D_problem then
    L  = localise3d(t(ie,:));
    ke = truss3dKe(p(t(ie,:),:),A(ie),E(ie));
  else
    L  = localise2d(t(ie,:));
    ke = truss2dKe(p(t(ie,:),:),A(ie),E(ie));
  end
  ue = U(L);
  fe = ke*ue';
  R(L) = R(L) + fe;
  if _3D_problem then
    [L,cx,cy,cz] = EltLen(p(t(ie,:),:));
    F(ie,:) = -( cx*fe(1) + cy*fe(2) + cz*fe(3));
  else
    [L,cx,cy] = EltLen(p(t(ie,:),:));
    F(ie,:) = -( cx*fe(1) + cy*fe(2) );
  end
end
endfunction
