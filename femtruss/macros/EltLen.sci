function [ds,c,s,p] = EltLen(XY)
//
// [c,s] = EltLen(XY)  
// Calcul cos et sin de la matrice de rotation pour un  
// element barre a deux noeuds (x1,y1) (x2,y2)
// c : cosinus de l'angle
// s : sinus de l'angle 
// XY: cordonnees des noeuds  XY = [x1,y1; x2,y2]
// 
// A.Seghir, 06/08/04

[nargout,nargin] = argn();

_3D_problem = (nargout==4);

if _3D_problem then
  dx = XY(2,1) - XY(1,1);
  dy = XY(2,2) - XY(1,2);
  dz = XY(2,3) - XY(1,3);

  ds = sqrt(sum((XY(1,:) - XY(2,:)).^2));
  c = dx/ds;
  s = dy/ds;
  p = dz/ds;
else
  dx = XY(2,1) - XY(1,1);
  dy = XY(2,2) - XY(1,2);

  ds = sqrt(sum((XY(1,:) - XY(2,:)).^2));
  c = dx/ds;
  s = dy/ds;
end
endfunction
