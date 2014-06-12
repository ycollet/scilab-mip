global EvalFObj;

Order = 1; // 1 2 or 4 for derivative

x0 = [25e-4; // Section element 1
      25e-4; // Section element 2
      25e-4; // Section element 3
      25e-4; // Section element 4
      25e-4; // Section element 5
      25e-4; // Section element 6
      25e-4; // Section element 7
      25e-4; // Section element 8
      25e-4; // Section element 9
      25e-4; // Section element 10
      25e-4; // Section element 11
      25e-4; // Section element 12
      25e-4; // Section element 13
      25e-4; // Section element 14
      25e-4  // Section element 15
     ];
//x0 = x0.*((1 - 0.1)*rand(size(x0,1),size(x0,2)) + 0.1);

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
     2, 3; // Connectivite element 2
     3, 4; // Connectivite element 3
     4, 5; // Connectivite element 4
     6, 7; // Connectivite element 5
     7, 8; // Connectivite element 6
     1, 6; // Connectivite element 7
     8, 5; // Connectivite element 8
     6, 2; // Connectivite element 9
     7, 3; // Connectivite element 10
     8, 4; // Connectivite element 11
     6, 3; // Connectivite element 12
     7, 2; // Connectivite element 13
     7, 4; // Connectivite element 14
     8, 3  // Connectivite element 15
    ];

//      6    7    8
//      +----+----+
//     /|\  /|\  /|\  
//    / | /\ | /\ | \ 
//   /  |/  \|/  \|  \
// -+---+----+----+---+-
//  |1  2    3    4   |5

p = [0.0, 0.0; // Noeud 1
     1.0, 0.0; // Noeud 2
     2.0, 0.0; // Noeud 3
     3.0, 0.0; // Noeud 4
     4.0, 0.0; // Noeud 5
     1.0, 1.0; // Noeud 6
     2.0, 1.0; // Noeud 7
     3.0, 1.0  // Noeud 8
    ];


e = [
     [1, 1] .* localise2d(5) // On localise les positions du noeud 6 dans les matrices et on immobilise les deux degrees de liberte du noeud 5
     [1, 1] .* localise2d(1) // On localise les positions du noeud 6 dans les matrices et on immobilise les deux degrees de liberte du noeud 2
    ]; 

A = [x(1); // Section element 1
     x(2); // Section element 2
     x(3); // Section element 3
     x(4); // Section element 4
     x(5); // Section element 5
     x(6); // Section element 6
     x(7); // Section element 7
     x(8); // Section element 8
     x(9); // Section element 9
     x(10); // Section element 10
     x(11); // Section element 11
     x(12); // Section element 12
     x(13); // Section element 13
     x(14); // Section element 14
     x(15)  // Section element 15
    ];

E   = ones(1,size(t,1)) * 210e9; // module d'elasticite des elements 
rho = ones(1,size(t,1)) * 7.8e3; // masse volumique 

F = [0;  0     // noeud 1
     0; -3.4e5 // noeud 2
     0; -3.4e5 // noeud 3
     0; -3.4e5 // noeud 4
     0;  0     // noeud 5
     0;  0     // noeud 6
     0;  0     // noeud 7
     0;  0     // noeud 8
    ];
endfunction

function y = fobj_truss(x)
global EvalFObj;
EvalFObj = EvalFObj + 1;
[t,p,e,A,E,rho,F] = bridge_optim(x);
[U,P,R]= femtruss(bridge_optim, %F, x);
// First objective: minimize the deformation at nodes 2, 3, 4 with respect to y

// The deck of the bridge
Pos2 = localise2d(2);
Pos3 = localise2d(3);
Pos4 = localise2d(4);

Deformation = sqrt(U(Pos2(2))^2 + U(Pos3(2))^2 + U(Pos4(2))^2);
y = Deformation;
endfunction

function dy = dfobj_truss(x)
dy = derivative(fobj_truss,x,order=Order)';
endfunction

function [y, dy, ind] = optim_fobj_truss(x, ind)
y  = fobj_truss(x);
dy = dfobj_truss(x)';
endfunction

// Parameters for optim
EvalFObj = 0;
Solutions = list();
Solutions_Cost = [];
Solutions_Deformation = [];
Algorithm    = 'gc';
upper_bound  = 2500*ones(size(x0,1),size(x0,2));
lower_bound  = 2.5e-5*ones(size(x0,1),size(x0,2));
TOL          = 1.0e-6; // accuracy for convergence test (minimum)
MaxEvalFunc  = 400;

printf('optimization starting, be patient ... \n\n');

[f_opt, x_opt] = optim(optim_fobj_truss, 'b', lower_bound, upper_bound, x0, algo=Algorithm,'ar',MaxEvalFunc,MaxEvalFunc,TOL,TOL);

printf('number of call to objective function: %d\n',EvalFObj);

printf('initial solution:'); disp(x0');
printf('initial objective function value = %f\n',fobj_truss(x0));

printf('Final solution:'); disp(x_opt');
printf('Final objective function value = %f\n',fobj_truss(x_opt));

scf();
[t,p,e,A,E,rho,F] = bridge_optim(x_opt);
plotsection(t,p,A);
title('sections of the optimal solution');

