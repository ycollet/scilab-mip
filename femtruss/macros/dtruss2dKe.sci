function [dKdx1,dKdx2,dKdy1,dKdy2] = dtruss2dKe(XY,A,E,Computation)
// [dKdx1,dKdx2,dKdy1,dKdy2] = dtruss2dKe(XY,A,E)
// Calcul de la matrice elementaire des derivees partielles pour un  
// element barre a deux noeuds (x1,y1) (x2,y2)
// A : section de l'element
// E : module d'elasticite 
// XY: cordonnees des noeuds  XY = [x1,y1; x2,y2]
// A. Seghir, 06/08/04
// Y. Collette 01/10/07

[L,c,s] = EltLen(XY);  

x1 = XY(1,1);
x2 = XY(2,1);
y1 = XY(1,2);
y2 = XY(2,2);

dKdx1 = [];
dKdx2 = [];
dKdy1 = [];
dKdy2 = [];

// Original matrix

cc = c*c;
cs = c*s;
ss = s*s;

ke = [ cc,  cs, -cc, -cs; ...
       cs,  ss, -cs, -ss; ...
      -cc, -cs,  cc,  cs; ...
      -cs, -ss,  cs,  ss
     ];

// Maxima commands to generate the stiffness partial derivative matrix
// display2d:false;
// C(x1,x2,y1,y2):= (x2-x1)/sqrt((x1-x2)^2+(y1-y2)^2);
// S(x1,x2,y1,y2):= (y2-y1)/sqrt((x1-x2)^2+(y1-y2)^2);
// K(x1,x2,y1,y2):=A*E/sqrt((x1-x2)^2+(y1-y2)^2)*matrix([C(x1,x2,y1,y2)^2,               C(x1,x2,y1,y2)*S(x1,x2,y1,y2), -C(x1,x2,y1,y2)^2,              -C(x1,x2,y1,y2)*S(x1,x2,y1,y2) ], ...
//                                                      [C(x1,x2,y1,y2)*S(x1,x2,y1,y2),  S(x1,x2,y1,y2)^2,              -C(x1,x2,y1,y2)*S(x1,x2,y1,y2), -S(x1,x2,y1,y2)^2              ], ...
//                                                      [-C(x1,x2,y1,y2)^2,             -C(x1,x2,y1,y2)*S(x1,x2,y1,y2),  C(x1,x2,y1,y2)^2,               C(x1,x2,y1,y2)*S(x1,x2,y1,y2) ], ...
//                                                      [-C(x1,x2,y1,y2)*S(x1,x2,y1,y2),-S(x1,x2,y1,y2)^2,               C(x1,x2,y1,y2)*S(x1,x2,y1,y2),  S(x1,x2,y1,y2)^2              ]);
//
// Result:
// * dK/dx1
// dKdx1:diff(K(x1,x2,y1,y2),x1);
// dKdx1:factor(dKdx1);

dX = x2-x1;
dY = y2-y1;
L  = dY^2+dX^2;

dC2dx1 = (2*dX^3-2*dX*L)/L^2;
dCSdx1 = (2*dX^2*dY-dY*L)/L^2;
dS2dx1 = (2*dX*dY^2)/L^2;
dC2dx2 = - dC2dx1;
dCSdx2 = - dCSdx1;
dS2dx2 = - dS2dx1;
dC2dy1 = (2*dX^2*dY)/L^2;
dCSdy1 = (2*dX*dY^2 - dX*L)/L^2;
dS2dy1 = (2*dY^3-2*dY*L)/L^2
dC2dy2 = - dC2dy1;
dCSdy2 = - dCSdy1;
dS2dy2 = - dS2dy1;

if Computation=='xy1' | Computation=='both' then
  dKdx1 = (A*E/L)*[dC2dx1,  dCSdx1, -dC2dx1, -dCSdx1; ...
                   dCSdx1,  dS2dx1, -dCSdx1, -dS2dx1; ...
                  -dC2dx1, -dCSdx1,  dC2dx1,  dCSdx1; ...
                  -dCSdx1, -dS2dx1,  dCSdx1,  dS2dx1] + (A*E)*(dX)/(L^2)^(3/2)*ke;
end

// * dK/dx2
// dKdx2:diff(K(x1,x2,y1,y2),x2);
// dKdx2:factor(dKdx2);

if Computation=='xy2' | Computation=='both' then
  dKdx2 = (A*E/L)*[dC2dx2,  dCSdx2, -dC2dx2, -dCSdx2; ...
                   dCSdx2,  dS2dx2, -dCSdx2, -dS2dx2; ...
                  -dC2dx2, -dCSdx2,  dC2dx2,  dCSdx2; ...
                  -dCSdx2, -dS2dx2,  dCSdx2,  dS2dx2] - (A*E)*(dX)/(L^2)^(3/2)*ke;
end

// * dK/dy1
// dKdy1:diff(K(x1,x2,y1,y2),y1);
// dKdy1:factor(dKdy1);

if Computation=='xy1' | Computation=='both' then
  dKdy1 = (A*E/L)*[dC2dy1,  dCSdy1, -dC2dy1, -dCSdy1; ...
                   dCSdy1,  dS2dy1, -dCSdy1, -dS2dy1; ...
                  -dC2dy1, -dCSdy1,  dC2dy1,  dCSdy1; ...
                  -dCSdy1, -dS2dy1,  dCSdy1,  dS2dy1] + (A*E)*(dY)/(L^2)^(3/2)*ke;
end

// * dK/dy2
// dKdy2:diff(K(x1,x2,y1,y2),y2);
// dKdy2:factor(dKdy2);

if Computation=='xy2' | Computation=='both' then
  dKdy2 = (A*E/L)*[dC2dy2,  dCSdy2, -dC2dy2, -dCSdy2; ...
                   dCSdy2,  dS2dy2, -dCSdy2, -dS2dy2; ...
                  -dC2dy2, -dCSdy2,  dC2dy2,  dCSdy2; ...
                  -dCSdy2, -dS2dy2,  dCSdy2,  dS2dy2] - (A*E)*(dY)/(L^2)^(3/2)*ke;
end

if Computation=='dx' then
// C(dX,dY):= dX/sqrt(dX^2+dY^2);
// S(dX,dY):= dY/sqrt(dX^2+dY^2);
// K(dX,dY):=A*E/sqrt(dX^2+dY^2)*matrix([C(dX,dY)^2,         C(dX,dY)*S(x1,x2,y1,y2), -C(dX,dY)^2,        -C(dX,dY)*S(dX,dY) ], ...
//                                      [C(dX,dY)*S(dX,dY),  S(dX,dY)^2,              -C(dX,dY)*S(dX,dY), -S(dX,dY)^2        ], ...
//                                      [-C(dX,dY)^2,       -C(dX,dY)*S(x1,x2,y1,y2),  C(dX,dY)^2,         C(dX,dY)*S(dX,dY) ], ...
//                                      [-C(dX,dY)*S(dX,dY),-S(dX,dY)^2,               C(dX,dY)*S(dX,dY),  S(dX,dY)^2        ]);

//  dX = abs(dX);
//  dY = abs(dY);
  
  // * dK/dx
  // dKdx:diff(K(dX,dY),dX);
  // dKdx:factor(dKdx);

  dKdx1 = [2*A*dX*E/(dY^2+dX^2)^(3/2)-3*A*dX^3*E/(dY^2+dX^2)^(5/2), A*dY*E/(dY^2+dX^2)^(3/2)-3*A*dX^2*dY*E/(dY^2+dX^2)^(5/2), 3*A*dX^3*E/(dY^2+dX^2)^(5/2)-2*A*dX*E/(dY^2+dX^2)^(3/2), 3*A*dX^2*dY*E/(dY^2+dX^2)^(5/2)-A*dY*E/(dY^2+dX^2)^(3/2); ...
           A*dY*E/(dY^2+dX^2)^(3/2)-3*A*dX^2*dY*E/(dY^2+dX^2)^(5/2), -3*A*dX*dY^2*E/(dY^2+dX^2)^(5/2), 3*A*dX^2*dY*E/(dY^2+dX^2)^(5/2)-A*dY*E/(dY^2+dX^2)^(3/2), 3*A*dX*dY^2*E/(dY^2+dX^2)^(5/2); ...
           3*A*dX^3*E/(dY^2+dX^2)^(5/2)-2*A*dX*E/(dY^2+dX^2)^(3/2), 3*A*dX^2*dY*E/(dY^2+dX^2)^(5/2)-A*dY*E/(dY^2+dX^2)^(3/2), 2*A*dX*E/(dY^2+dX^2)^(3/2)-3*A*dX^3*E/(dY^2+dX^2)^(5/2), A*dY*E/(dY^2+dX^2)^(3/2)-3*A*dX^2*dY*E/(dY^2+dX^2)^(5/2); ...
           3*A*dX^2*dY*E/(dY^2+dX^2)^(5/2)-A*dY*E/(dY^2+dX^2)^(3/2), 3*A*dX*dY^2*E/(dY^2+dX^2)^(5/2), A*dY*E/(dY^2+dX^2)^(3/2)-3*A*dX^2*dY*E/(dY^2+dX^2)^(5/2), -3*A*dX*dY^2*E/(dY^2+dX^2)^(5/2)];

  // * dK/dy
  // dKdy:diff(K(dX,dY),dY);
  // dKdy:factor(dKdy);

  dKdy1 = [-3*A*dX^2*dY*E/(dY^2+dX^2)^(5/2), A*dX*E/(dY^2+dX^2)^(3/2)-3*A*dX*dY^2*E/(dY^2+dX^2)^(5/2), 3*A*dX^2*dY*E/(dY^2+dX^2)^(5/2), 3*A*dX*dY^2*E/(dY^2+dX^2)^(5/2)-A*dX*E/(dY^2+dX^2)^(3/2); ...
           A*dX*E/(dY^2+dX^2)^(3/2)-3*A*dX*dY^2*E/(dY^2+dX^2)^(5/2),2*A*dY*E/(dY^2+dX^2)^(3/2)-3*A*dY^3*E/(dY^2+dX^2)^(5/2), 3*A*dX*dY^2*E/(dY^2+dX^2)^(5/2)-A*dX*E/(dY^2+dX^2)^(3/2), 3*A*dY^3*E/(dY^2+dX^2)^(5/2)-2*A*dY*E/(dY^2+dX^2)^(3/2); ...
           3*A*dX^2*dY*E/(dY^2+dX^2)^(5/2),3*A*dX*dY^2*E/(dY^2+dX^2)^(5/2)-A*dX*E/(dY^2+dX^2)^(3/2),-3*A*dX^2*dY*E/(dY^2+dX^2)^(5/2), A*dX*E/(dY^2+dX^2)^(3/2)-3*A*dX*dY^2*E/(dY^2+dX^2)^(5/2);...
           3*A*dX*dY^2*E/(dY^2+dX^2)^(5/2)-A*dX*E/(dY^2+dX^2)^(3/2),3*A*dY^3*E/(dY^2+dX^2)^(5/2)-2*A*dY*E/(dY^2+dX^2)^(3/2),A*dX*E/(dY^2+dX^2)^(3/2)-3*A*dX*dY^2*E/(dY^2+dX^2)^(5/2), 2*A*dY*E/(dY^2+dX^2)^(3/2)-3*A*dY^3*E/(dY^2+dX^2)^(5/2)];
end
endfunction
