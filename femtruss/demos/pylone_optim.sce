printf('Optimization of a pole structure\n\n');

global EvalFObj;
global Weight;

EvalFObj = 0;

// On ne touche pas aux noeuds 1, 2, 16, 17, 18, 19 - les 4 derniers supportents les cables
// les  deux premiers correspondent a la prise du pylone sur la terre

// 7                    + (15)
// 6               + (13)    + (14)
// 5        + (20) + (11)    + (12) + (21)
// 4 + (16) + (17) + (9)     + (10) + (18) + (19)
// 3               + (7)     + (8)
// 2               + (5)     + (6)
// 1               + (3)     + (4)
// 0               + (1)     + (2)
//   0      1      2         3      4      5

x0 = [2.0; 1.0; ... // OK
      4.0; 1.0; ... // OK
      2.0; 2.0; ... // OK
      4.0; 2.0; ... // OK
      2.0; 3.0; ... // OK
      4.0; 3.0; ... // OK
      2.0; 4.0; ... // OK
      4.0; 4.0; ... // OK
      2.0; 5.0; ... // OK
      4.0; 5.0; ... // OK
      2.0; 6.0; ... // OK
      4.0; 6.0; ... // OK
      3.0; 7.0; ... // OK
      1.0; 5.0; ... // OK
      5.0; 5.0];    // OK

function [t,p,e,A,E,rho,F] = pylone_optim(x)
// t   : donnees en lignes regroupees par 2 (connection extremite 1 - connection extremite 2)
// p   : autant de lignes que de noeuds. Sur chaque ligne n on retrouve les coordonnees 2d du noeud n
// e   : liste des noeuds d'appui
// A   : liste des sections des elements. Vecteur contenant autant d'entrees que d'elements
// E   : liste des modules d'elasticite des elements. Vecteur contenant autant d'entrees que d'elements
// rho : liste des masses volumiques des elements. Vecteur contenant autant d'entrees que d'elements
// F   : liste des forces appliquees aux noeuds. Vecteur colonne comprenant la coordonnes X du noeud 1 en premiere ligne, la coordonnees Y du noeud 1 en seconde 
//       ligne, etc ...

t = [1,  2; // Les barres horizontales
     3,  4;
     5,  6;
     7,  8;
     9,  10;
     11, 12;
     13, 14;
     16, 17;
     17, 9;
     10, 18;
     18, 19;
     20, 11;
     12, 21;
     1,  3; // Les barres verticales
     2,  4;
     3,  5;
     4,  6;
     5,  7;
     6,  8;
     7,  9;
     8,  10;
     9,  11;
     10, 12;
     11, 13;
     12, 14;
     17, 20;
     18, 21;
     1,  4; // Les barres diagonales droites
     3,  6;
     5,  8;
     7,  10;
     9,  12;
     11, 14;
     13, 15;
     17, 11;
     20, 13;
     16, 20;
     10, 21;
     2,  3; // Les barres de diagonales gauches
     4,  5;
     6,  7;
     8,  9;
     10, 11;
     12, 13;
     14, 15;
     19, 21;
     18, 12;
     21, 14;
     9,  20
     ];

p = [2.0,   0.0;   // Node 1
     4.0,   0.0;   // Node 2
     x(1),  x(2);  // Node 3
     x(3),  x(4);  // Node 4
     x(5),  x(6);  // Node 5
     x(7),  x(8);  // Node 6
     x(9),  x(10); // Node 7
     x(11), x(12); // Node 8
     x(13), x(14); // Node 9
     x(15), x(16); // Node 10
     x(17), x(18); // Node 11
     x(19), x(20); // Node 12
     x(21), x(22); // Node 13
     x(23), x(24); // Node 14
     x(25), x(26); // Node 15
     0.0,   4.0;   // Node 16
     1.0,   4.0;   // Node 17
     5.0,   4.0;   // Node 18
     6.0,   4.0;   // Node 19
     x(27), x(28); // Node 20
     x(29), x(30)  // Node 21
    ];

// 7                    + (15)
// 6               + (13)    + (14)
// 5        + (20) + (11)    + (12) + (21)
// 4 + (16) + (17) + (9)     + (10) + (18) + (19)
// 3               + (7)     + (8)
// 2               + (5)     + (6)
// 1               + (3)     + (4)
// 0               + (1)     + (2)
//   0      1      2         3      4      5


e = [
     [1, 1] .* localise2d(1) // On localise les positions du noeud 1 dans les matrices et on immobilise les deux degrees de liberte du noeud 1 et 2
     [1, 1] .* localise2d(2) // On localise les positions du noeud 2 dans les matrices et on immobilise les deux degrees de liberte du noeud 1 et 2
    ]; 

A   = ones(1,size(t,1)) * 25e-4; // sections des elements
E   = ones(1,size(t,1)) * 210e9; // module d'elasticite des elements 
rho = ones(1,size(t,1)) * 7.8e3; // masse volumique 

// Le pylone soutient 4 cables de 100 Kg chacun
F = [0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; -1e5;
     0.0; -1e5;
     0.0; -1e5;
     0.0; -1e5;
     0.0; 0.0;
     0.0; 0.0
    ];
endfunction

function y = fobj_truss(x)
global EvalFObj;
global Weight;
k1 = 1.0;
k2 = 0.005;
EvalFObj = EvalFObj + 1;
[t,p,e,A,E,rho,F] = pylone_optim(x);
[U,P,R]= femtruss(pylone_optim, Log, x);
// First objective: minimize the deformation at nodes 16, 17, 18, 19
n = length(U);
U_aux = matrix([U(1:2:n),U(2:2:n)],size(p,1),size(p,2));
Deformation = sqrt(sum(U_aux(16,:).^2 + U_aux(17,:).^2 + U_aux(18,:).^2 + U_aux(19,:).^2));
// Second objective: minimize the total length of the bars
Cost = 0;
for i=1:size(t,1)
  Cost = Cost + sum((p(t(i,1),:) - p(t(i,2),:)).^2);
end
y = (1-Weight)*Deformation*k1 + Weight*Cost*k2;
endfunction

function dy = dfobj_truss(x)
dy = derivative(fobj_truss,x,h=1e-3)';
endfunction

function [y, dy, ind] = optim_fobj_truss(x, ind)
y  = fobj_truss(x);
dy = derivative(fobj_truss,x,h=1e-3);
endfunction

function plot_fobj_truss(x)
[t,p,e,A,E,rho,F] = pylone_optim(x);
[U,P,R]= femtruss(pylone_optim, Log, x);
plotdeforme(U,p,t,10);
endfunction

/////////////////////////////////
// Symetric objective function //
/////////////////////////////////

x0_sym = [2.0; 1.0; ...
          2.0; 2.0; ...
          2.0; 3.0; ...
          2.0; 4.0; ...
          2.0; 5.0; ...
          2.0; 6.0; ...
          2.5; 7.0; ...
          1.0; 5.0];

function [t,p,e,A,E,rho,F] = pylone_optim_sym(x)
// t   : donnees en lignes regroupees par 2 (connection extremite 1 - connection extremite 2)
// p   : autant de lignes que de noeuds. Sur chaque ligne n on retrouve les coordonnees 2d du noeud n
// e   : liste des noeuds d'appui
// A   : liste des sections des elements. Vecteur contenant autant d'entrees que d'elements
// E   : liste des modules d'elasticite des elements. Vecteur contenant autant d'entrees que d'elements
// rho : liste des masses volumiques des elements. Vecteur contenant autant d'entrees que d'elements
// F   : liste des forces appliquees aux noeuds. Vecteur colonne comprenant la coordonnes X du noeud 1 en premiere ligne, la coordonnees Y du noeud 1 en seconde 
//       ligne, etc ...

t = [1,  2; // Les barres horizontales
     3,  4;
     5,  6;
     7,  8;
     9,  10;
     11, 12;
     13, 14;
     16, 17;
     17, 9;
     10, 18;
     18, 19;
     20, 11;
     12, 21;
     1,  3; // Les barres verticales
     2,  4;
     3,  5;
     4,  6;
     5,  7;
     6,  8;
     7,  9;
     8,  10;
     9,  11;
     10, 12;
     11, 13;
     12, 14;
     17, 20;
     18, 21;
     1,  4; // Les barres diagonales droites
     3,  6;
     5,  8;
     7,  10;
     9,  12;
     11, 14;
     13, 15;
     17, 11;
     20, 13;
     16, 20;
     10, 21;
     2,  3; // Les barres de diagonales gauches
     4,  5;
     6,  7;
     8,  9;
     10, 11;
     12, 13;
     14, 15;
     19, 21;
     18, 12;
     21, 14;
     9,  20
     ];

p = [2.0,       0.0;   // Node 1
     3.0,       0.0;   // Node 2
     x(1),      x(2);  // Node 3
     5 - x(1),  x(2);  // Node 4
     x(3),      x(4);  // Node 5
     5 - x(3),  x(4);  // Node 6
     x(5),      x(6);  // Node 7
     5 - x(5),  x(6);  // Node 8
     x(7),      x(8);  // Node 9
     5 - x(7),  x(8);  // Node 10
     x(9),      x(10); // Node 11
     5 - x(9),  x(10); // Node 12
     x(11),     x(12); // Node 13
     5 - x(11), x(12); // Node 14
     x(13),     x(14); // Node 15
     0.0,       4.0;   // Node 16
     1.0,       4.0;   // Node 17
     4.0,       4.0;   // Node 18
     5.0,       4.0;   // Node 19
     x(15),     x(16); // Node 20
     5 - x(15), x(16)  // Node 21
    ];

// 7                    + (15)
// 6               + (13)    + (14)
// 5        + (20) + (11)    + (12) + (21)
// 4 + (16) + (17) + (9)     + (10) + (18) + (19)
// 3               + (7)     + (8)
// 2               + (5)     + (6)
// 1               + (3)     + (4)
// 0               + (1)     + (2)
//   0      1      2         3      4      5


e = [
     [1, 1] .* localise2d(1) // On localise les positions du Node 1 dans les matrices et on immobilise les deux degrees de liberte du noeud 1 et 2
     [1, 1] .* localise2d(2) // On localise les positions du noeud 2 dans les matrices et on immobilise les deux degrees de liberte du noeud 1 et 2
    ]; 

A   = ones(1,size(t,1)) * 25e-4; // sections des elements
E   = ones(1,size(t,1)) * 210e9; // module d'elasticite des elements 
rho = ones(1,size(t,1)) * 7.8e3; // masse volumique 

// Le pylone soutient 4 cables de 100 Kg chacun
F = [0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; -1e5;
     0.0; -1e5;
     0.0; -1e5;
     0.0; -1e5;
     0.0; 0.0;
     0.0; 0.0
    ];
endfunction

function y = fobj_truss_sym(x)
global EvalFObj;
global Weight;
k1 = 1.0;
k2 = 0.005;
EvalFObj = EvalFObj + 1;
[t,p,e,A,E,rho,F] = pylone_optim_sym(x);
[U,P,R]= femtruss(pylone_optim_sym, Log, x);
// First objective: minimize the deformation at nodes 16, 17, 18, 19
n = length(U);
U_aux = matrix([U(1:2:n),U(2:2:n)],size(p,1),size(p,2));
Deformation = sqrt(sum(U_aux(16,:).^2 + U_aux(17,:).^2 + U_aux(18,:).^2 + U_aux(19,:).^2));
// Second objective: minimize the total length of the bars
Cost = 0;
for i=1:size(t,1)
  Cost = Cost + sum((p(t(i,1),:) - p(t(i,2),:)).^2);
end
y = (1-Weight)*Deformation*k1 + Weight*Cost*k2;
endfunction

function dy = dfobj_truss_sym(x)
dy = derivative(fobj_truss_sym,x,h=1e-3)';
endfunction

function [y, dy, ind] = optim_fobj_truss_sym(x, ind)
y  = fobj_truss_sym(x);
dy = derivative(fobj_truss_sym,x,h=1e-3);
endfunction

function plot_fobj_truss_sym(x)
[t,p,e,A,E,rho,F] = pylone_optim_sym(x);
[U,P,R]= femtruss(pylone_optim_sym, Log, x);
plotdeforme(U,p,t,10);
endfunction

// 7                    + (15)
// 6               + (13)    + (14)
// 5        + (20) + (11)    + (12) + (21)
// 4 + (16) + (17) + (9)     + (10) + (18) + (19)
// 3               + (7)     + (8)
// 2               + (5)     + (6)
// 1               + (3)     + (4)
// 0               + (1)     + (2)
//   0      1      2         3      4      5


//p = [2.0,       0.0;   // 2.0,   0.0;   // Node 1
//     4.0,       0.0;   // 4.0,   0.0;   // Node 2
//     x(1),      x(2);  // x(1),  x(2);  // Node 3
//     5 - x(1),  x(2);  // x(3),  x(4);  // Node 4
//     x(3),      x(4);  // x(5),  x(6);  // Node 5
//     5 - x(3),  x(4);  // x(7),  x(8);  // Node 6
//     x(5),      x(6);  // x(9),  x(10); // Node 7
//     5 - x(5),  x(6);  // x(11), x(12); // Node 8
//     x(7),      x(8);  // x(13), x(14); // Node 9
//     5 - x(7),  x(8);  // x(15), x(16); // Node 10
//     x(9),      x(10); // x(17), x(18); // Node 11
//     5 - x(9),  x(10); // x(19), x(20); // Node 12
//     x(11),     x(12); // x(21), x(22); // Node 13
//     5 - x(11), x(12); // x(23), x(24); // Node 14
//     x(13),     x(14); // x(25), x(26); // Node 15
//     0.0,       4.0;   // 0.0,   4.0;   // Node 16
//     1.0,       4.0;   // 1.0,   4.0;   // Node 17
//     5.0,       4.0;   // 5.0,   4.0;   // Node 18
//     6.0,       4.0;   // 6.0,   4.0;   // Node 19
//     x(15),     x(16); // x(27), x(28); // Node 20
//     5 - x(15), x(16)  // x(29), x(30)  // Node 21
//    ];

MaxEvalFunc     = 1000;
Algorithm       = 'gc';    // 'qn', 'gc', 'nd' -> Ne marche qu'avec 'qn' (quasi-newton). Pour les autres, on obtient rapidement une structure mal
                          // conditionnée
TOL             = 1.0e-12; // accuracy for convergence test (minimum)
Log             = %F;

NbPtsToCompute  = 3;
MaxWeight       = 0.5;
MaxMinStepCount = 3;

UseSymetric = %F; // Use the symetric formulation of the problem
rand_ampl   = 0.001; // Amplitude of the value value we add to x0 -> to try ta have several solutions

///////////////////////////
// Starting optimization //
///////////////////////////

EvalFObj = 0;
Solutions = list();
Solutions_Cost = [];
Solutions_Deformation = [];

printf('optimization starting, be patient ... \n\n');

if (UseSymetric) then
  x0_sym = x0_sym + rand_ampl*rand(size(x0_sym,1),size(x0_sym,2));
else
  x0 = x0 + rand_ampl*rand(size(x0,1),size(x0,2));
end

if (UseSymetric) then
  printf('using symetric parametrization\n');
  Weight = 1; Solutions_Cost($+1)        = fobj_truss_sym(x0_sym);
  Weight = 0; Solutions_Deformation($+1) = fobj_truss_sym(x0_sym);
  printf('initial solution:'); disp(x0_sym');
else
  printf('using non-symetric parametrization\n');
  Weight = 1; Solutions_Cost($+1)        = fobj_truss(x0);
  Weight = 0; Solutions_Deformation($+1) = fobj_truss(x0);
  printf('initial solution:'); disp(x0');
end

printf('initial objective function value = %f\n',Solutions_Deformation($));

for i=0:NbPtsToCompute-1
  // Weigh will vary from 0 to MaxWeight
  printf('Optimisation of point %d / %d\n', i+1, NbPtsToCompute+1);
  if (NbPtsToCompute==1) then
    Weight = 0;
  else
    Weight = MaxWeight*i/(NbPtsToCompute-1);
  end
  
  if (UseSymetric) then
    [f_opt, x_opt] = optim(optim_fobj_truss_sym, x0_sym, algo=Algorithm, iter=MaxEvalFunc);
  else
    [f_opt, x_opt] = optim(optim_fobj_truss, x0, algo=Algorithm, iter=MaxEvalFunc);
  end

  Solutions($+1) = x_opt;

  printf('Final solution:'); disp(x_opt');
  
  if (UseSymetric) then
    printf('Final objective function value = %f\n',fobj_truss_sym(x_opt));

    Weight = 1; Solutions_Cost($+1)        = fobj_truss_sym(x_opt);
    Weight = 0; Solutions_Deformation($+1) = fobj_truss_sym(x_opt);
  else
    printf('Final objective function value = %f\n',fobj_truss(x_opt));

    Weight = 1; Solutions_Cost($+1)        = fobj_truss(x_opt);
    Weight = 0; Solutions_Deformation($+1) = fobj_truss(x_opt);
  end
  printf('number of call to objective function: %d\n',EvalFObj);
end

scf();
plot(Solutions_Cost, Solutions_Deformation, 'ro');
xtitle('Pareto front','Length of bars','Deformation');

scf();
if (UseSymetric) then
  plot_fobj_truss_sym(x0_sym);
else
  plot_fobj_truss(x0);
end
xtitle('Before optimization','x','y');

scf();
if (NbPtsToCompute==1) then
  if (UseSymetric) then
    plot_fobj_truss_sym(Solutions(1));
  else
    plot_fobj_truss(Solutions(1));
  end
  xtitle(sprintf('After optimization - Weight = %f', Weight),'x','y');
else
  Index = 1;
  NbCol = ceil(sqrt(NbPtsToCompute));
  for i=1:NbPtsToCompute
    Weight = MaxWeight*i/(NbPtsToCompute-1);
    subplot(NbCol,NbCol,i);
    if (UseSymetric) then
      plot_fobj_truss_sym(Solutions(i));
    else
      plot_fobj_truss(Solutions(i));
    end
    xtitle(sprintf('After optimization - Weight = %f', Weight),'x','y');
  end
end

