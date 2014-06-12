UseAnalytical = %T;
Order  = 2; // 1,2,3 for derivative   
x0     = [1.0; 1.0];

function [t,p,e,A,E,rho,F] = bridge_optim(x)
// t : donnees en lignes regroupees par 2 (connection extremite 1 - connection extremite 2)
// p : autant de lignes que de noeuds. Sur chaque ligne n on retrouve les coordonnees 2d du noeud n
// e : liste des noeuds d'appui
// A : liste des sections des elements. Vecteur contenant autant d'entrees que d'elements
// E : liste des modules d'elasticite des elements. Vecteur contenant autant d'entrees que d'elements
// rho : liste des masses volumiques des elements. Vecteur contenant autant d'entrees que d'elements
// F   : liste des forces appliquees aux noeuds. Vecteur colonne comprenant la coordonnes X du noeud 1 en premiere ligne, la coordonnees Y du noeud 1 en seconde 
//       ligne, etc ...

t = [1, 2; // Connectivite element 1
     2, 3  // Connectivite element 2
    ];

//      2
//      +
//     / \
//    /   \
//   /     \
// -+       +-
//  |1      |3

p = [0.0,  0.0;  // Noeud 1
     x(1), x(2); // Noeud 2
     2.0,  0.0;  // Noeud 3
    ];


e = [
     [1, 1] .* localise2d(1) // On localise les positions du noeud 1 dans les matrices et on immobilise les deux degrees de liberte du noeud 1
     [1, 1] .* localise2d(3) // On localise les positions du noeud 3 dans les matrices et on immobilise les deux degrees de liberte du noeud 3
    ]; 

A   = ones(1,size(t,1)) * 25e-4; // sections des elements
E   = ones(1,size(t,1)) * 210e9; // module d'elasticite des elements 
rho = ones(1,size(t,1)) * 7.8e3; // masse volumique 

F = [0;  0     // noeud 1
     0; -3.4e5 // noeud 2
     0;  0     // noeud 3
    ];
endfunction

function y = fobj_truss(x)
[t,p,e,A,E,rho,F] = bridge_optim(x);
[U,P,R]= femtruss(bridge_optim, %F, x);
// First objective: minimize the deformation at nodes 2, 3, 4 with respect to y

// The deck of the bridge
Pos2 = localise2d(2);

y = sqrt(sum(U(Pos2).^2));
endfunction

function dy = dfobj_truss(x)
if UseAnalytical then
  // Here, we will compute analytical derivatives for the objective function.
  // We use dfemtruss to compute analytical derivatives for displacement
  // We then compute by hand the analytical derivative for the objective function using the result of dfemtruss
  [t,p,e,A,E,rho,F] = bridge_optim(x);
  [U,P,R]= femtruss(bridge_optim, %F, x);

  dU = dfemtruss_ana(bridge_optim,U,[],%F, x);
    
  // The deck and The degree of freedom for the optimization
  Pos2 = localise2d(2);
  // dU contains only partial derivatives for nodes 2 3 and 4. So, we get partial derivatives with respect to y by accessing
  // value 2 4 and 6. Values 1 3 and 5 are partial derivatives wrt x.
  // dU([I1 I2 I3],I4) = dU(I1) / dI4 dU(I2) / dI4 dU(I3) / dI4
//  dy(1) = sum(dU([Pos2],Pos2(1)).*U([Pos2])'); // For x2
//  dy(2) = sum(dU([Pos2],Pos2(2)).*U([Pos2])'); // For y2
  dy(1) = sum(dU(Pos2(1),Pos2).*U(Pos2)); // For x
  dy(2) = sum(dU(Pos2(2),Pos2).*U(Pos2)); // For y2
    
  dy = dy / sqrt(sum(U(Pos2).^2));
else
  dy = derivative(fobj_truss,x,order=Order)';
end
endfunction

printf('initial solution:'); disp(x0');
UseAnalytical = %T;
printf('analytical derivative - initial objective function value = %f\n',dfobj_truss(x0));
UseAnalytical = %F;
printf('finite difference derivative - initial objective function value = %f\n',dfobj_truss(x0));

scf();
x2     = 0.1:0.01:1.9;
x_pos  = 3.2;
dy_an = zeros(2,length(x2));
dy_fd = zeros(2,length(x2));
for i=1:length(x2);
  y(i) = fobj_truss([x_pos x2(i)]);
  UseAnalytical = %T;
  dy_an(:,i) = dfobj_truss([x_pos x2(i)]');
  UseAnalytical = %F;
  dy_fd(:,i) = dfobj_truss([x_pos x2(i)]');
end

subplot(4,1,1);
plot(x2,y,'k-');
xtitle('objective function','x2','fobj');

subplot(4,1,2);
plot(x2(1:$-1),diff(y)/0.01,'r-');
xtitle('derivative of objective function','x2','fobj');

subplot(4,1,3);
plot(x2,dy_an(2,:),'g-');
xtitle('Analytical derivative of objective function','x2','dfobj/dx2');

subplot(4,1,4);
plot(x2,dy_fd(2,:),'r-');
xtitle('Finite difference derivative of objective function','x2','dfobj/dx2');

