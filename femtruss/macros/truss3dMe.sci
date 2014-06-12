function Me = truss3dMe(XY,A,rho)
// Me = truss3dMe(XY,A,rho)
//      Calcul de la matrice masse elementaire Me pour un  
//      element barre a deux noeuds (x1,y1,z1) (x2,y2,z2)
// A   : section de l'element
// rho : masse volumique de l'element 
// XY  : cordonnees des noeuds  XY = [x1,y1,z1; x2,y2,z2]
//
// A. Seghir, 07/08/04

L = EltLen(XY);
Me = 0.5*rho*A*L *eye(6,6); // Lumped mass matrix

// Consistent mass matrix
// pour la masse repartie : 
//Me =(rho*A*L/6) * [ 2  0  1  0 ; ...
//                    0  2  0  1 ; ...
//                    1  0  2  0 ; ...
//                    0  1  0  2 ; ...
//                ];
endfunction
