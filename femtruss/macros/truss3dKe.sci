function ke = truss3dKe(XY,A,E)
// ke = truss3dKe(XY,A,E)
// Calcul de la matrice elementaire pour un  
// element barre a deux noeuds (x1,y1,z1) (x2,y2,z2)
// A : section de l'element
// E : module d'elasticite 
// XY: cordonnees des noeuds  XY = [x1,y1,z1; x2,y2,z2]
// A. Seghir, 06/08/04
[L,cx,cy,cz] = EltLen(XY);  

ke = (A*E/L) * [ cx^2,   cx*cy,  cx*cz, -cx^2,  -cx*cy, -cx*cz; ...
                 cx*cy,  cy^2,   cy*cz, -cx*cy, -cy^2,  -cy*cz; ...
                 cx*cz,  cy*cz,  cz^2,  -cx*cz, -cy*cz, -cz^2; ...
                -cx^2,  -cx*cy, -cx*cz,  cx^2,   cx*cy,  cx*cz; ...
                -cx*cy, -cy^2,  -cy*cz,  cx*cy,  cy^2,   cy*cz; ...
                -cx*cz, -cy*cz, -cz^2,   cx*cz,  cy*cz,  cz^2];
endfunction



