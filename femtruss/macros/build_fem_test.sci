function [t, p, e, A, E, rho, F] = build_fem_test(test_name)
// Most of these structures are extracted from
// L. Lamberti, C. Pappalettere, "Move limits definition in structural optimization with 
// sequential linear programming. Part II: Numerical examples", computer and structures 81 (2003), pp. 215-238

// t : donnees en lignes regroupees par 2 (connection extremite 1 - connection extremite 2)
// p : autant de lignes que de nodes. Sur chaque ligne n on retrouve les coordonnees 2d du node n
// e : liste des nodes d'appui
// A : liste des sections des elements. Vecteur contenant autant d'entrees que d'elements
// E : liste des modules d'elasticite des elements. Vecteur contenant autant d'entrees que d'elements
// rho : liste des masses volumiques des elements. Vecteur contenant autant d'entrées que d'elements
// F   : liste des forces appliquees aux nodes. Vecteur colonne comprenant la coordonnes X du node 1 en premiere ligne, la coordonnees Y du node 1 en seconde 
//       ligne, etc ...

select (test_name)
  case '2bars' then
    t = [3, 1; ...
         2, 3  ...
        ];

    // 1+         +2
    //   \       /
    //  1  \   /  2
    //       +
    //       3

    p = [0.0, 0.6; // node 1
         1.6, 0.6; // node 2
         0.8, 0.0  // node 3
        ];

    e = [
          [1 1].*localise2d(1)
          [1 1].*localise2d(2)
        ];

    A   = [1 1] * 25e-4; // sections of the elements
    E   = [1 1] * 210e9; // elasticity module of the elements
    rho = [1 1] * 7.8e3; // volumic mass

    F = [ 0;  0     // node 1
          0;  0     // node 2
          0; -3.4e6 // node 3
        ];    
  case 'truss' then
    t = [7,  6; ...
         8,  7; ...
         9,  8; ...
         2,  1; ...
         3,  2; ...
         4,  3; ...
         5,  4; ...
         6,  2; ...
         7,  1; ...
         6,  2; ...
         8,  2; ...
         7,  3; ...
         9,  3; ...
         8,  4; ...
         9,  5  ...
        ];

    //  6    7    8    9
    // -+----+----+----+
    //   \  / \  / \  / \    
    //    /\   /\   /\   \  
    //   /  \ /  \ /  \   \ 
    // -+----+----+----+---+
    //  1    2    3    4   5

    p = [0.0, 0.0; // Node 1
         1.0, 0.0; // Node 2
         2.0, 0.0; // Node 3
         3.0, 0.0; // Node 4
         4.0, 0.0; // Node 5
         0.0, 1.0; // Node 6
         1.0, 1.0; // Node 7
         2.0, 1.0; // Node 8
         3.0, 1.0  // Node 9
        ];


    e = [
         [1, 1] .* localise2d(6) // We localise the positions of node 6 in the matrix and we fix both degrees of freedom of the node
         [1, 1] .* localise2d(1) // We localise the positions of node 1 in the matrix and we fix both degrees of freedom of the node
        ]; 

    A   = ones(1,size(t,1)) * 25e-4;
    E   = ones(1,size(t,1)) * 210e9;
    rho = ones(1,size(t,1)) * 7.8e3;

    F = [0;  0     // node 1
         0;  0     // node 2
         0;  0     // node 3
         0;  0     // node 4
         0; -3.4e5 // node 5
         0;  0     // node 6
         0;  0     // node 7
         0;  0     // node 8
         0;  0     // node 9
        ];
  case 'dev5' then
    t = [1, 2;
         1, 3;
         2, 3 
        ];
    //
    //      + 3
    //     / \
    //  2 /   \ 3
    //   /     \
    //  /       \
    // +---------+
    // 1    1    2

    p = [0.0, 0.0; // node 1
         5.0, 0.0; // node 2
         3.0, 1.5  // node 3
        ];

    e = [
         [0, 1] .* localise2d(1) // On recupere la liste des nodes de l'elements 1 et on immobilise le second node
         [1, 1] .* localise2d(2) // On recupere la liste des nodes de l'elements 2 et on immobilise les deux nodes
        ];

    A   = [25, 49, 49]*1e-4 ; 
    E   = [1, 1, 1] * 210e9;
    rho = [1, 1, 1] * 7.8e3;

    F = [ 0;  0     // node 1
          0;  0     // node 2
          0; -950e3 // node 3
        ];
  case 'bridge2d' then
    //////////////////////
    // Bridge structure //
    //////////////////////

    t = [1, 2;
         2, 3;
         3, 4;
         4, 5;
         6, 7;
         7, 8;
         1, 6;
         8, 5;
         6, 2;
         7, 3;
         8, 4;
         6, 3;
         7, 2;
         7, 4;
         8, 3
        ];

    //      6    7    8
    //      +----+----+
    //     /|\  /|\  /|\  
    //    / | /\ | /\ | \ 
    //   /  |/  \|/  \|  \
    // -+---+----+----+---+-
    //  |1  2    3    4   |5

    p = [0.0, 0.0; // node 1
         1.0, 0.0; // node 2
         2.0, 0.0; // node 3
         3.0, 0.0; // node 4
         4.0, 0.0; // node 5
         1.0, 1.0; // node 6
         2.0, 1.0; // node 7
         3.0, 1.0  // node 8
        ];


    e = [
         [1, 1] .* localise2d(5)
         [1, 1] .* localise2d(1)
        ]; 

    A   = ones(1,size(t,1)) * 25e-4;
    E   = ones(1,size(t,1)) * 210e9;
    rho = ones(1,size(t,1)) * 7.8e3;

    F = [0;  0     // node 1
         0; -3.4e5 // node 2
         0; -3.4e5 // node 3
         0; -3.4e5 // node 4
         0;  0     // node 5
         0;  0     // node 6
         0;  0     // node 7
         0;  0     // node 8
        ];
  case 'bridge3d' then
    //////////////////////
    // Bridge structure //
    //////////////////////

    t = [1+0, 2+0; // Side 1
         2+0, 3+0;
         3+0, 4+0;
         4+0, 5+0;
         6+0, 7+0;
         7+0, 8+0;
         1+0, 6+0;
         8+0, 5+0;
         6+0, 2+0;
         7+0, 3+0;
         8+0, 4+0;
         6+0, 3+0;
         7+0, 2+0;
         7+0, 4+0;
         8+0, 3+0;
         1+8, 2+8; // Side 2
         2+8, 3+8;
         3+8, 4+8;
         4+8, 5+8;
         6+8, 7+8;
         7+8, 8+8;
         1+8, 6+8;
         8+8, 5+8;
         6+8, 2+8;
         7+8, 3+8;
         8+8, 4+8;
         6+8, 3+8;
         7+8, 2+8;
         7+8, 4+8;
         8+8, 3+8;
         1,   9; // Connect side 1 and 2
         2,   10;
         3,   11;
         4,   12;
         5,   13;
         6,   14;
         7,   15;
         8,   16;
         9,   2; // Crosses between side 1 and 2
         10,  1;
         2,   11; // Cross
         10,  3;
         3,   12; // Cross
         11,  4;
         4,   13; // Cross
         12,  5;
         9,   6; // Cross
         1,   14;
         6,   15; // Cross
         7,   14;
         7,   16; // Cross
         8,   15;
         8,   13; // Cross
         5,   16;
         6,   10; // Cross
         2,   14;
         7,   11; // Cross
         3,   15;
         8,   12; // Cross
         4,   16
         ];

    // Side 1
    //      6    7    8
    //      +----+----+
    //     /|\  /|\  /|\  
    //    / | /\ | /\ | \ 
    //   /  |/  \|/  \|  \
    // -+---+----+----+---+-
    //  |1  2    3    4   |5

    // Side 2
    //      14   15   16
    //      +----+----+
    //     /|\  /|\  /|\  
    //    / | /\ | /\ | \ 
    //   /  |/  \|/  \|  \
    // -+---+----+----+---+-
    //  |9  10   11   12  |13


    p = [0.0, 0.0, 0.0; // node 1
         1.0, 0.0, 0.0; // node 2
         2.0, 0.0, 0.0; // node 3
         3.0, 0.0, 0.0; // node 4
         4.0, 0.0, 0.0; // node 5
         1.0, 1.0, 0.0; // node 6
         2.0, 1.0, 0.0; // node 7
         3.0, 1.0, 0.0  // node 8
         0.0, 0.0, 1.0; // node 9
         1.0, 0.0, 1.0; // node 10
         2.0, 0.0, 1.0; // node 11
         3.0, 0.0, 1.0; // node 12
         4.0, 0.0, 1.0; // node 13
         1.0, 1.0, 1.0; // node 14
         2.0, 1.0, 1.0; // node 15
         3.0, 1.0, 1.0  // node 16
        ];


    e = [
         [1, 1, 1] .* localise3d(5)  // We localise the positions of node 6 in the matrix and we fix the degrees of fredom of the node
         [1, 1, 1] .* localise3d(1)
         [1, 1, 1] .* localise3d(13)
         [1, 1, 1] .* localise3d(9)
        ]; 

    A   = ones(1,size(t,1)) * 25e-4;
    E   = ones(1,size(t,1)) * 210e9;
    rho = ones(1,size(t,1)) * 7.8e3;

    F = [0;  0; 0      // node 1
         0;  -3.4e5; 0 // node 2
         0;  -3.4e5; 0 // node 3
         0;  -3.4e5; 0 // node 4
         0;  0; 0      // node 5
         0;  0; 0      // node 6
         0;  0; 0      // node 7
         0;  0; 0      // node 8
         0;  0; 0      // node 9
         0;  -3.4e5; 0 // node 10
         0;  -3.4e5; 0 // node 11
         0;  -3.4e5; 0 // node 12
         0;  0; 0      // node 13
         0;  0; 0      // node 14
         0;  0; 0      // node 15
         0;  0; 0      // node 16
         ];
  case 'bar' then
    ///////////////////
    // Bar structure //
    ///////////////////

    //      6    5    4
    //      +----+----+
    //       \  /|\  /|
    //        /\ | /\ |
    //       /  \|/  \|
    //      +----+----+
    //      1    2    3

    L = 360; // inches

    p = [0,   0;   // Node 1
         1*L, 0;   // Node 2
         2*L, 0;   // Node 3
         2*L, 1*L; // Node 4
         1*L, 1*L; // Node 5
         0,   1*L  // Node 6
        ];

    t = [1, 2;
         2, 3;
         3, 4;
         4, 5;
         5, 6;
         1, 5;
         6, 2;
         5, 2;
         2, 4;
         5, 3;
         4, 3
        ];

    e = [
         [1, 1] .* localise2d(6) // We localise the position of node 6 in the matrix and we fix the dergee of freedom of this node.
         [1, 1] .* localise2d(1) // We localise the position of node 1 in the matrix and we fix the dergee of freedom of this node.
        ]; 

    A   = ones(1,size(t,1)) * 25e-4; // sections of the elements
    E   = ones(1,size(t,1)) * 210e9; // elasticity module of the elements
    rho = ones(1,size(t,1)) * 7.8e3; // volumic mass

    F = [0;  0     // node 1
         0;  0     // node 2
         0; -3.4e5 // node 3
         0;  0     // node 4
         0;  0     // node 5
         0;  0     // node 6
        ];
  case 'pylon2d' then

    t = [1,  2; // Horizontal bars
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
         1,  3; // Vertical bars
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
         1,  4; // The right diagonal bars
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
         2,  3; // The left diagonal bars
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

    p = [2.0, 0.0; // node 1
         4.0, 0.0; // node 2
         2.0, 1.0; // node 3
         4.0, 1.0; // node 4
         2.0, 2.0; // node 5
         4.0, 2.0; // node 6
         2.0, 3.0; // node 7
         4.0, 3.0; // node 8
         2.0, 4.0; // node 9
         4.0, 4.0; // node 10
         2.0, 5.0; // node 11
         4.0, 5.0; // node 12
         2.0, 6.0; // node 13
         4.0, 6.0; // node 14
         3.0, 7.0; // node 15
         0.0, 4.0; // node 16
         1.0, 4.0; // node 17
         5.0, 4.0; // node 18
         6.0, 4.0; // node 19
         1.0, 5.0; // node 20
         5.0, 5.0  // node 21
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
         [1, 1] .* localise2d(1)
         [1, 1] .* localise2d(2)
        ]; 

    A   = ones(1,size(t,1)) * 25e-4;
    E   = ones(1,size(t,1)) * 210e9;
    rho = ones(1,size(t,1)) * 7.8e3;

    // The pylon handles 4 wire weighting 100 Kg each
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
  case 'pylon3d' then
    ////////////////////////////////////////
    // 3D pylone structure (with 25 bars) //
    ////////////////////////////////////////

    a = 25; // inches

    p = [-1.5*a,  0,     8*a;   // Node 1
          1.5*a,  0,     8*a;   // Node 2
         -1.5*a,  1.5*a, 4*a;   // Node 3
          1.5*a,  1.5*a, 4*a;   // Node 4
          1.5*a, -1.5*a, 4*a;   // Node 5
         -1.5*a, -1.5*a, 4*a;   // Node 6
           -4*a,    4*a,   0;   // Node 7
            4*a,    4*a,   0;   // Node 8
            4*a,   -4*a,   0;   // Node 9
           -4*a,   -4*a,   0    // Node 10
        ];

    t = [1, 2;
         1, 4;
         2, 3;
         1, 5;
         2, 6;
         2, 4;
         2, 5;
         1, 3;
         1, 6;
         3, 6;
         5, 4;
         6, 5;
         3, 4;
         6, 9;
         5, 9;
         4, 9;
         5, 8;
         4, 8;
         3, 8;
         4, 7;
         3, 7;
         6, 7;
         3, 10;
         6, 10;
         5, 10
        ];

    e = [
         [1, 1, 1] .* localise3d(7)  // We localise the position of node 7  in the matrix and we fix the dergee of freedom of this node.
         [1, 1, 1] .* localise3d(8)  // We localise the position of node 8  in the matrix and we fix the dergee of freedom of this node.
         [1, 1, 1] .* localise3d(9)  // We localise the position of node 9  in the matrix and we fix the dergee of freedom of this node.
         [1, 1, 1] .* localise3d(10) // We localise the position of node 10 in the matrix and we fix the dergee of freedom of this node.
        ]; 

    A   = ones(1,size(t,1)) * 25e-4; // sections of the elements
    E   = ones(1,size(t,1)) * 210e9; // elasticity module of the elements
    rho = ones(1,size(t,1)) * 7.8e3; // volumic mass

    F = [0; 0; -3.4e5; // node 1
         0; 0; -3.4e5; // node 2
         0; 0;  0;     // node 3
         0; 0;  0;     // node 4
         0; 0;  0;     // node 5
         0; 0;  0;     // node 6
         0; 0;  0;     // node 7
         0; 0;  0;     // node 8
         0; 0;  0;     // node 9
         0; 0;  0     // node 10
        ];
  case 'building3d' then
    //////////////////////////////////////////
    // 3D building structure (with 72 bars) //
    //////////////////////////////////////////

    b = 60; // inches

    p = [  0,   0, 4*b;   // Node 1
         2*b,   0, 4*b;   // Node 2
         2*b, 2*b, 4*b;   // Node 3
           0, 2*b, 4*b;   // Node 4
           0,   0, 3*b;   // Node 5
         2*b,   0, 3*b;   // Node 6
         2*b, 2*b, 3*b;   // Node 7
           0, 2*b, 3*b;   // Node 8
           0,   0, 2*b;   // Node 9
         2*b,   0, 2*b;   // Node 10
         2*b, 2*b, 2*b;   // Node 11
           0, 2*b, 2*b;   // Node 12
           0,   0, 1*b;   // Node 13
         2*b,   0, 1*b;   // Node 14
         2*b, 2*b, 1*b;   // Node 15
           0, 2*b, 1*b;   // Node 16
           0,   0, 0*b;   // Node 17
         2*b,   0, 0*b;   // Node 18
         2*b, 2*b, 0*b;   // Node 19
           0, 2*b, 0*b    // Node 20
           ];

    t = [1+0,  2+0; // Floor 4
         2+0,  3+0;
         3+0,  4+0;
         4+0,  1+0;
         4+0,  2+0;
         3+0,  1+0;
         2+0,  6+0;
         3+0,  7+0;
         4+0,  8+0;
         1+0,  5+0;
         1+0,  6+0;
         2+0,  5+0;
         2+0,  7+0;
         3+0,  6+0;
         4+0,  7+0;
         3+0,  8+0;
         4+0,  5+0;
         1+0,  8+0;
         1+4,  2+4; // Floor 3
         2+4,  3+4;
         3+4,  4+4;
         4+4,  1+4;
         4+4,  2+4;
         3+4,  1+4;
         2+4,  6+4;
         3+4,  7+4;
         4+4,  8+4;
         1+4,  5+4;
         1+4,  6+4;
         2+4,  5+4;
         2+4,  7+4;
         3+4,  6+4;
         4+4,  7+4;
         3+4,  8+4;
         4+4,  5+4;
         1+4,  8+4;
         1+8,  2+8; // Floor 2
         2+8,  3+8;
         3+8,  4+8;
         4+8,  1+8;
         4+8,  2+8;
         3+8,  1+8;
         2+8,  6+8;
         3+8,  7+8;
         4+8,  8+8;
         1+8,  5+8;
         1+8,  6+8;
         2+8,  5+8;
         2+8,  7+8;
         3+8,  6+8;
         4+8,  7+8;
         3+8,  8+8;
         4+8,  5+8;
         1+8,  8+8;
         1+12, 2+12; // Floor 2
         2+12, 3+12;
         3+12, 4+12;
         4+12, 1+12;
         4+12, 2+12;
         3+12, 1+12;
         2+12, 6+12;
         3+12, 7+12;
         4+12, 8+12;
         1+12, 5+12;
         1+12, 6+12;
         2+12, 5+12;
         2+12, 7+12;
         3+12, 6+12;
         4+12, 7+12;
         3+12, 8+12;
         4+12, 5+12;
         1+12, 8+12
        ];

    e = [
         [1, 1, 1] .* localise3d(17) // We localise the position of node 17 in the matrix and we fix the dergee of freedom of this node.
         [1, 1, 1] .* localise3d(18) // We localise the position of node 18 in the matrix and we fix the dergee of freedom of this node.
         [1, 1, 1] .* localise3d(19) // We localise the position of node 19 in the matrix and we fix the dergee of freedom of this node.
         [1, 1, 1] .* localise3d(20) // We localise the position of node 20 in the matrix and we fix the dergee of freedom of this node.
        ]; 

    A   = ones(1,size(t,1)) * 25e-4; // sections of the elements
    E   = ones(1,size(t,1)) * 210e9; // elasticity module of the elements
    rho = ones(1,size(t,1)) * 7.8e3; // volumic mass

    F = [0; 0; -3.4e5; // node 1
         0; 0; -3.4e5; // node 2
         0; 0; -3.4e5; // node 3
         0; 0; -3.4e5; // node 4
         0; 0; -3.4e5; // node 5
         0; 0; -3.4e5; // node 6
         0; 0; -3.4e5; // node 7
         0; 0; -3.4e5; // node 8
         0; 0; -3.4e5; // node 9
         0; 0; -3.4e5; // node 10
         0; 0; -3.4e5; // node 11
         0; 0; -3.4e5; // node 12
         0; 0; -3.4e5; // node 13
         0; 0; -3.4e5; // node 14
         0; 0; -3.4e5; // node 15
         0; 0; -3.4e5; // node 16
         0; 0; 0; // node 17
         0; 0; 0; // node 18
         0; 0; 0; // node 19
         0; 0; 0 // node 20
        ];
  case 'building2d' then
  ///////////////////////////////////////////
  // 2D building structure (with 200 bars) //
  ///////////////////////////////////////////

  b = 60; // inches

  p = [-480, 1800; // Floor 6
       -240, 1800;
          0, 1800;
        240, 1800;
        480, 1800;
       -480, 1656;
       -360, 1656;
       -240, 1656;
       -120, 1656;
          0, 1656;
        120, 1656;
        240, 1656;
        360, 1656;
        480, 1656;
       -480, 1512; // Floor 5
       -240, 1512;
          0, 1512;
        240, 1512;
        480, 1512;
       -480, 1368;
       -360, 1368;
       -240, 1368;
       -120, 1368;
          0, 1368;
        120, 1368;
        240, 1368;
        360, 1368;
        480, 1368;
       -480, 1224; // Floor 4
       -240, 1224;
          0, 1224;
        240, 1224;
        480, 1224;
       -480, 1080;
       -360, 1080;
       -240, 1080;
       -120, 1080;
          0, 1080;
        120, 1080;
        240, 1080;
        360, 1080;
        480, 1080;
       -480,  936; // Floor 3
       -240,  936;
          0,  936;
        240,  936;
        480,  936;
       -480,  792;
       -360,  792;
       -240,  792;
       -120,  792;
          0,  792;
        120,  792;
        240,  792;
        360,  792;
        480,  792;
       -480,  648; // Floor 2
       -240,  648;
          0,  648;
        240,  648;
        480,  648;
       -480,  504;
       -360,  504;
       -240,  504;
       -120,  504;
          0,  504;
        120,  504;
        240,  504;
        360,  504;
        480,  504;
       -480,  360; // Floor 1
       -240,  360;
          0,  360;
        240,  360;
        480,  360;
       -240,    0;
        240,    0
         ];

  t = [ 1+0,   2+0; // Floor 6
        2+0,   3+0;
        3+0,   4+0;
        4+0,   5+0;
        1+0,   6+0;
        1+0,   7+0;
        2+0,   7+0;
        2+0,   8+0;
        2+0,   9+0;
        3+0,   9+0;
        3+0,  10+0;
        3+0,  11+0;
        4+0,  11+0;
        4+0,  12+0;
        4+0,  13+0;
        5+0,  13+0;
        5+0,  14+0;
        6+0,   7+0;
        7+0,   8+0;
        8+0,   9+0;
        9+0,  10+0;
       10+0,  11+0;
       11+0,  12+0;
       12+0,  13+0;
       13+0,  14+0;
        6+0,  15+0;
        7+0,  15+0;
        7+0,  16+0;
        8+0,  16+0;
        9+0,  16+0;
        9+0,  17+0;
       10+0,  17+0;
       11+0,  17+0;
       11+0,  18+0;
       12+0,  18+0;
       13+0,  18+0;
       13+0,  19+0;
       14+0,  19+0;
        1+14,  2+14; // Floor 5
        2+14,  3+14;
        3+14,  4+14;
        4+14,  5+14;
        1+14,  6+14;
        1+14,  7+14;
        2+14,  7+14;
        2+14,  8+14;
        2+14,  9+14;
        3+14,  9+14;
        3+14, 10+14;
        3+14, 11+14;
        4+14, 11+14;
        4+14, 12+14;
        4+14, 13+14;
        5+14, 13+14;
        5+14, 14+14;
        6+14,  7+14;
        7+14,  8+14;
        8+14,  9+14;
        9+14, 10+14;
       10+14, 11+14;
       11+14, 12+14;
       12+14, 13+14;
       13+14, 14+14;
        6+14, 15+14;
        7+14, 15+14;
        7+14, 16+14;
        8+14, 16+14;
        9+14, 16+14;
        9+14, 17+14;
       10+14, 17+14;
       11+14, 17+14;
       11+14, 18+14;
       12+14, 18+14;
       13+14, 18+14;
       13+14, 19+14;
       14+14, 19+14;
        1+28,  2+28; // Floor 4
        2+28,  3+28;
        3+28,  4+28;
        4+28,  5+28;
        1+28,  6+28;
        1+28,  7+28;
        2+28,  7+28;
        2+28,  8+28;
        2+28,  9+28;
        3+28,  9+28;
        3+28, 10+28;
        3+28, 11+28;
        4+28, 11+28;
        4+28, 12+28;
        4+28, 13+28;
        5+28, 13+28;
        5+28, 14+28;
        6+28,  7+28;
        7+28,  8+28;
        8+28,  9+28;
        9+28, 10+28;
       10+28, 11+28;
       11+28, 12+28;
       12+28, 13+28;
       13+28, 14+28;
        6+28, 15+28;
        7+28, 15+28;
        7+28, 16+28;
        8+28, 16+28;
        9+28, 16+28;
        9+28, 17+28;
       10+28, 17+28;
       11+28, 17+28;
       11+28, 18+28;
       12+28, 18+28;
       13+28, 18+28;
       13+28, 19+28;
       14+28, 19+28;
        1+42,  2+42; // Floor 3
        2+42,  3+42;
        3+42,  4+42;
        4+42,  5+42;
        1+42,  6+42;
        1+42,  7+42;
        2+42,  7+42;
        2+42,  8+42;
        2+42,  9+42;
        3+42,  9+42;
        3+42, 10+42;
        3+42, 11+42;
        4+42, 11+42;
        4+42, 12+42;
        4+42, 13+42;
        5+42, 13+42;
        5+42, 14+42;
        6+42,  7+42;
        7+42,  8+42;
        8+42,  9+42;
        9+42, 10+42;
       10+42, 11+42;
       11+42, 12+42;
       12+42, 13+42;
       13+42, 14+42;
        6+42, 15+42;
        7+42, 15+42;
        7+42, 16+42;
        8+42, 16+42;
        9+42, 16+42;
        9+42, 17+42;
       10+42, 17+42;
       11+42, 17+42;
       11+42, 18+42;
       12+42, 18+42;
       13+42, 18+42;
       13+42, 19+42;
       14+42, 19+42;
        1+56,  2+56; // Floor 2
        2+56,  3+56;
        3+56,  4+56;
        4+56,  5+56;
        1+56,  6+56;
        1+56,  7+56;
        2+56,  7+56;
        2+56,  8+56;
        2+56,  9+56;
        3+56,  9+56;
        3+56, 10+56;
        3+56, 11+56;
        4+56, 11+56;
        4+56, 12+56;
        4+56, 13+56;
        5+56, 13+56;
        5+56, 14+56;
        6+56,  7+56;
        7+56,  8+56;
        8+56,  9+56;
        9+56, 10+56;
       10+56, 11+56;
       11+56, 12+56;
       12+56, 13+56;
       13+56, 14+56;
        6+56, 15+56;
        7+56, 15+56;
        7+56, 16+56;
        8+56, 16+56;
        9+56, 16+56;
        9+56, 17+56;
       10+56, 17+56;
       11+56, 17+56;
       11+56, 18+56;
       12+56, 18+56;
       13+56, 18+56;
       13+56, 19+56;
       14+56, 19+56;
          71,    72; // Floor 1
          72,    73;
          73,    74;
          74,    75;
          71,    76;
          72,    76;
          73,    76;
          73,    77;
          74,    77;
          75,    77
  ];

  e = [
       [1, 1] .* localise2d(76) // We localise the position of node 76 in the matrix and we fix the dergee of freedom of this node.
       [1, 1] .* localise2d(77) // We localise the position of node 77 in the matrix and we fix the dergee of freedom of this node.
      ]; 

  A   = ones(1,size(t,1)) * 25e-4; // sections of the elements
  E   = ones(1,size(t,1)) * 210e9; // elasticity module of the elements
  rho = ones(1,size(t,1)) * 7.8e3; // volumic mass

  F = [0; -3.4e5 // node 1
       0; -3.4e5 // node 2
       0; -3.4e5 // node 3
       0; -3.4e5 // node 4
       0; -3.4e5 // node 5
       0; -3.4e5 // node 6
       0; -3.4e5 // node 7
       0; -3.4e5 // node 8
       0; -3.4e5 // node 9
       0; -3.4e5 // node 10
       0; -3.4e5 // node 11
       0; -3.4e5 // node 12
       0; -3.4e5 // node 13
       0; -3.4e5 // node 14
       0; -3.4e5 // node 15
       0; -3.4e5 // node 16
       0; -3.4e5 // node 17
       0; -3.4e5 // node 18
       0; -3.4e5 // node 19
       0; -3.4e5 // node 20
       0; -3.4e5 // node 21
       0; -3.4e5 // node 22
       0; -3.4e5 // node 23
       0; -3.4e5 // node 24
       0; -3.4e5 // node 25
       0; -3.4e5 // node 26
       0; -3.4e5 // node 27
       0; -3.4e5 // node 28
       0; -3.4e5 // node 29
       0; -3.4e5 // node 30
       0; -3.4e5 // node 31
       0; -3.4e5 // node 32
       0; -3.4e5 // node 33
       0; -3.4e5 // node 34
       0; -3.4e5 // node 35
       0; -3.4e5 // node 36
       0; -3.4e5 // node 37
       0; -3.4e5 // node 38
       0; -3.4e5 // node 39
       0; -3.4e5 // node 40
       0; -3.4e5 // node 41
       0; -3.4e5 // node 42
       0; -3.4e5 // node 43
       0; -3.4e5 // node 44
       0; -3.4e5 // node 45
       0; -3.4e5 // node 46
       0; -3.4e5 // node 47
       0; -3.4e5 // node 48
       0; -3.4e5 // node 49
       0; -3.4e5 // node 50
       0; -3.4e5 // node 51
       0; -3.4e5 // node 52
       0; -3.4e5 // node 53
       0; -3.4e5 // node 54
       0; -3.4e5 // node 55
       0; -3.4e5 // node 56
       0; -3.4e5 // node 57
       0; -3.4e5 // node 58
       0; -3.4e5 // node 59
       0; -3.4e5 // node 60
       0; -3.4e5 // node 61
       0; -3.4e5 // node 62
       0; -3.4e5 // node 63
       0; -3.4e5 // node 64
       0; -3.4e5 // node 65
       0; -3.4e5 // node 66
       0; -3.4e5 // node 67
       0; -3.4e5 // node 68
       0; -3.4e5 // node 69
       0; -3.4e5 // node 70
       0; -3.4e5 // node 71
       0; -3.4e5 // node 72
       0; -3.4e5 // node 73
       0; -3.4e5 // node 74
       0; -3.4e5 // node 75
       0; 0      // node 76
       0; 0      // node 77
      ];
  case 'dome3d' then
    ///////////////////////
    // 3D dome structure //
    ///////////////////////

    deff('y=_f_dome(z)','y = 10*sqrt(20 - z)/sqrt(20);');
    r = 10;

    p = [0,            0,            20;                // Floor 11 - Node 1
         _f_dome(18)*cos(0*%pi/6),  _f_dome(18)*sin(0*%pi/6),  18; // Floor 10 - Node 2
         _f_dome(18)*cos(1*%pi/6),  _f_dome(18)*sin(1*%pi/6),  18;            // Node 3
         _f_dome(18)*cos(2*%pi/6),  _f_dome(18)*sin(2*%pi/6),  18;            // Node 4
         _f_dome(18)*cos(3*%pi/6),  _f_dome(18)*sin(3*%pi/6),  18;            // Node 5
         _f_dome(18)*cos(4*%pi/6),  _f_dome(18)*sin(4*%pi/6),  18;            // Node 6
         _f_dome(18)*cos(5*%pi/6),  _f_dome(18)*sin(5*%pi/6),  18;            // Node 7
         _f_dome(18)*cos(6*%pi/6),  _f_dome(18)*sin(6*%pi/6),  18;            // Node 8
         _f_dome(18)*cos(7*%pi/6),  _f_dome(18)*sin(7*%pi/6),  18;            // Node 9
         _f_dome(18)*cos(8*%pi/6),  _f_dome(18)*sin(8*%pi/6),  18;            // Node 10
         _f_dome(18)*cos(9*%pi/6),  _f_dome(18)*sin(9*%pi/6),  18;            // Node 11
         _f_dome(18)*cos(10*%pi/6), _f_dome(18)*sin(10*%pi/6), 18;            // Node 12
         _f_dome(18)*cos(11*%pi/6), _f_dome(18)*sin(11*%pi/6), 18;            // Node 13
         _f_dome(16)*cos(0*%pi/6),  _f_dome(16)*sin(0*%pi/6),  16; // Floor 9  - Node 14
         _f_dome(16)*cos(1*%pi/6),  _f_dome(16)*sin(1*%pi/6),  16;            // Node 15
         _f_dome(16)*cos(2*%pi/6),  _f_dome(16)*sin(2*%pi/6),  16;            // Node 16
         _f_dome(16)*cos(3*%pi/6),  _f_dome(16)*sin(3*%pi/6),  16;            // Node 17
         _f_dome(16)*cos(4*%pi/6),  _f_dome(16)*sin(4*%pi/6),  16;            // Node 18
         _f_dome(16)*cos(5*%pi/6),  _f_dome(16)*sin(5*%pi/6),  16;            // Node 19
         _f_dome(16)*cos(6*%pi/6),  _f_dome(16)*sin(6*%pi/6),  16;            // Node 20
         _f_dome(16)*cos(7*%pi/6),  _f_dome(16)*sin(7*%pi/6),  16;            // Node 21
         _f_dome(16)*cos(8*%pi/6),  _f_dome(16)*sin(8*%pi/6),  16;            // Node 22
         _f_dome(16)*cos(9*%pi/6),  _f_dome(16)*sin(9*%pi/6),  16;            // Node 23
         _f_dome(16)*cos(10*%pi/6), _f_dome(16)*sin(10*%pi/6), 16;            // Node 24
         _f_dome(16)*cos(11*%pi/6), _f_dome(16)*sin(11*%pi/6), 16;            // Node 25
         _f_dome(14)*cos(0*%pi/6),  _f_dome(14)*sin(0*%pi/6),  14; // Floor 8  - Node 26
         _f_dome(14)*cos(1*%pi/6),  _f_dome(14)*sin(1*%pi/6),  14;            // Node 27
         _f_dome(14)*cos(2*%pi/6),  _f_dome(14)*sin(2*%pi/6),  14;            // Node 28
         _f_dome(14)*cos(3*%pi/6),  _f_dome(14)*sin(3*%pi/6),  14;            // Node 29
         _f_dome(14)*cos(4*%pi/6),  _f_dome(14)*sin(4*%pi/6),  14;            // Node 30
         _f_dome(14)*cos(5*%pi/6),  _f_dome(14)*sin(5*%pi/6),  14;            // Node 31
         _f_dome(14)*cos(6*%pi/6),  _f_dome(14)*sin(6*%pi/6),  14;            // Node 32
         _f_dome(14)*cos(7*%pi/6),  _f_dome(14)*sin(7*%pi/6),  14;            // Node 33
         _f_dome(14)*cos(8*%pi/6),  _f_dome(14)*sin(8*%pi/6),  14;            // Node 34
         _f_dome(14)*cos(9*%pi/6),  _f_dome(14)*sin(9*%pi/6),  14;            // Node 35
         _f_dome(14)*cos(10*%pi/6), _f_dome(14)*sin(10*%pi/6), 14;            // Node 36
         _f_dome(14)*cos(11*%pi/6), _f_dome(14)*sin(11*%pi/6), 14;            // Node 37
         _f_dome(12)*cos(0*%pi/6),  _f_dome(12)*sin(0*%pi/6),  12; // Floor 7  - Node 38
         _f_dome(12)*cos(1*%pi/6),  _f_dome(12)*sin(1*%pi/6),  12;            // Node 39
         _f_dome(12)*cos(2*%pi/6),  _f_dome(12)*sin(2*%pi/6),  12;            // Node 40
         _f_dome(12)*cos(3*%pi/6),  _f_dome(12)*sin(3*%pi/6),  12;            // Node 41
         _f_dome(12)*cos(4*%pi/6),  _f_dome(12)*sin(4*%pi/6),  12;            // Node 42
         _f_dome(12)*cos(5*%pi/6),  _f_dome(12)*sin(5*%pi/6),  12;            // Node 43
         _f_dome(12)*cos(6*%pi/6),  _f_dome(12)*sin(6*%pi/6),  12;            // Node 44
         _f_dome(12)*cos(7*%pi/6),  _f_dome(12)*sin(7*%pi/6),  12;            // Node 45
         _f_dome(12)*cos(8*%pi/6),  _f_dome(12)*sin(8*%pi/6),  12;            // Node 46
         _f_dome(12)*cos(9*%pi/6),  _f_dome(12)*sin(9*%pi/6),  12;            // Node 47
         _f_dome(12)*cos(10*%pi/6), _f_dome(12)*sin(10*%pi/6), 12;            // Node 48
         _f_dome(12)*cos(11*%pi/6), _f_dome(12)*sin(11*%pi/6), 12;            // Node 49
         _f_dome(10)*cos(0*%pi/6),  _f_dome(10)*sin(0*%pi/6),  10; // Floor 6  - Node 50
         _f_dome(10)*cos(1*%pi/6),  _f_dome(10)*sin(1*%pi/6),  10;            // Node 51
         _f_dome(10)*cos(2*%pi/6),  _f_dome(10)*sin(2*%pi/6),  10;            // Node 52
         _f_dome(10)*cos(3*%pi/6),  _f_dome(10)*sin(3*%pi/6),  10;            // Node 53
         _f_dome(10)*cos(4*%pi/6),  _f_dome(10)*sin(4*%pi/6),  10;            // Node 54
         _f_dome(10)*cos(5*%pi/6),  _f_dome(10)*sin(5*%pi/6),  10;            // Node 55
         _f_dome(10)*cos(6*%pi/6),  _f_dome(10)*sin(6*%pi/6),  10;            // Node 56
         _f_dome(10)*cos(7*%pi/6),  _f_dome(10)*sin(7*%pi/6),  10;            // Node 57
         _f_dome(10)*cos(8*%pi/6),  _f_dome(10)*sin(8*%pi/6),  10;            // Node 58
         _f_dome(10)*cos(9*%pi/6),  _f_dome(10)*sin(9*%pi/6),  10;            // Node 59
         _f_dome(10)*cos(10*%pi/6), _f_dome(10)*sin(10*%pi/6), 10;            // Node 60
         _f_dome(10)*cos(11*%pi/6), _f_dome(10)*sin(11*%pi/6), 10;            // Node 61
         _f_dome(8)*cos(0*%pi/6),   _f_dome(8)*sin(0*%pi/6),   8; // Floor 5   - Node 62
         _f_dome(8)*cos(1*%pi/6),   _f_dome(8)*sin(1*%pi/6),   8;             // Node 63
         _f_dome(8)*cos(2*%pi/6),   _f_dome(8)*sin(2*%pi/6),   8;             // Node 64
         _f_dome(8)*cos(3*%pi/6),   _f_dome(8)*sin(3*%pi/6),   8;             // Node 65
         _f_dome(8)*cos(4*%pi/6),   _f_dome(8)*sin(4*%pi/6),   8;             // Node 66
         _f_dome(8)*cos(5*%pi/6),   _f_dome(8)*sin(5*%pi/6),   8;             // Node 67
         _f_dome(8)*cos(6*%pi/6),   _f_dome(8)*sin(6*%pi/6),   8;             // Node 68
         _f_dome(8)*cos(7*%pi/6),   _f_dome(8)*sin(7*%pi/6),   8;             // Node 69
         _f_dome(8)*cos(8*%pi/6),   _f_dome(8)*sin(8*%pi/6),   8;             // Node 70
         _f_dome(8)*cos(9*%pi/6),   _f_dome(8)*sin(9*%pi/6),   8;             // Node 71
         _f_dome(8)*cos(10*%pi/6),  _f_dome(8)*sin(10*%pi/6),  8;             // Node 72
         _f_dome(8)*cos(11*%pi/6),  _f_dome(8)*sin(11*%pi/6),  8;             // Node 73
         _f_dome(6)*cos(0*%pi/6),   _f_dome(6)*sin(0*%pi/6),   6; // Floor 4   - Node 74
         _f_dome(6)*cos(1*%pi/6),   _f_dome(6)*sin(1*%pi/6),   6;             // Node 75
         _f_dome(6)*cos(2*%pi/6),   _f_dome(6)*sin(2*%pi/6),   6;             // Node 76
         _f_dome(6)*cos(3*%pi/6),   _f_dome(6)*sin(3*%pi/6),   6;             // Node 77
         _f_dome(6)*cos(4*%pi/6),   _f_dome(6)*sin(4*%pi/6),   6;             // Node 78
         _f_dome(6)*cos(5*%pi/6),   _f_dome(6)*sin(5*%pi/6),   6;             // Node 79
         _f_dome(6)*cos(6*%pi/6),   _f_dome(6)*sin(6*%pi/6),   6;             // Node 80
         _f_dome(6)*cos(7*%pi/6),   _f_dome(6)*sin(7*%pi/6),   6;             // Node 81
         _f_dome(6)*cos(8*%pi/6),   _f_dome(6)*sin(8*%pi/6),   6;             // Node 82
         _f_dome(6)*cos(9*%pi/6),   _f_dome(6)*sin(9*%pi/6),   6;             // Node 83
         _f_dome(6)*cos(10*%pi/6),  _f_dome(6)*sin(10*%pi/6),  6;             // Node 84
         _f_dome(6)*cos(11*%pi/6),  _f_dome(6)*sin(11*%pi/6),  6;             // Node 85
         _f_dome(4)*cos(0*%pi/6),   _f_dome(4)*sin(0*%pi/6),   4; // Floor 3   - Node 86
         _f_dome(4)*cos(1*%pi/6),   _f_dome(4)*sin(1*%pi/6),   4;             // Node 87
         _f_dome(4)*cos(2*%pi/6),   _f_dome(4)*sin(2*%pi/6),   4;             // Node 88
         _f_dome(4)*cos(3*%pi/6),   _f_dome(4)*sin(3*%pi/6),   4;             // Node 89
         _f_dome(4)*cos(4*%pi/6),   _f_dome(4)*sin(4*%pi/6),   4;             // Node 90
         _f_dome(4)*cos(5*%pi/6),   _f_dome(4)*sin(5*%pi/6),   4;             // Node 91
         _f_dome(4)*cos(6*%pi/6),   _f_dome(4)*sin(6*%pi/6),   4;             // Node 92
         _f_dome(4)*cos(7*%pi/6),   _f_dome(4)*sin(7*%pi/6),   4;             // Node 93
         _f_dome(4)*cos(8*%pi/6),   _f_dome(4)*sin(8*%pi/6),   4;             // Node 94
         _f_dome(4)*cos(9*%pi/6),   _f_dome(4)*sin(9*%pi/6),   4;             // Node 95
         _f_dome(4)*cos(10*%pi/6),  _f_dome(4)*sin(10*%pi/6),  4;             // Node 96
         _f_dome(4)*cos(11*%pi/6),  _f_dome(4)*sin(11*%pi/6),  4;             // Node 97
         _f_dome(2)*cos(0*%pi/6),   _f_dome(2)*sin(0*%pi/6),   2; // Floor 2   - Node 98
         _f_dome(2)*cos(1*%pi/6),   _f_dome(2)*sin(1*%pi/6),   2;             // Node 99
         _f_dome(2)*cos(2*%pi/6),   _f_dome(2)*sin(2*%pi/6),   2;             // Node 100
         _f_dome(2)*cos(3*%pi/6),   _f_dome(2)*sin(3*%pi/6),   2;             // Node 101
         _f_dome(2)*cos(4*%pi/6),   _f_dome(2)*sin(4*%pi/6),   2;             // Node 102
         _f_dome(2)*cos(5*%pi/6),   _f_dome(2)*sin(5*%pi/6),   2;             // Node 103
         _f_dome(2)*cos(6*%pi/6),   _f_dome(2)*sin(6*%pi/6),   2;             // Node 104
         _f_dome(2)*cos(7*%pi/6),   _f_dome(2)*sin(7*%pi/6),   2;             // Node 105
         _f_dome(2)*cos(8*%pi/6),   _f_dome(2)*sin(8*%pi/6),   2;             // Node 106
         _f_dome(2)*cos(9*%pi/6),   _f_dome(2)*sin(9*%pi/6),   2;             // Node 107
         _f_dome(2)*cos(10*%pi/6),  _f_dome(2)*sin(10*%pi/6),  2;             // Node 108
         _f_dome(2)*cos(11*%pi/6),  _f_dome(2)*sin(11*%pi/6),  2;             // Node 109
         _f_dome(0)*cos(0*%pi/6),   _f_dome(0)*sin(0*%pi/6),   0; // Floor 1   - Node 110
         _f_dome(0)*cos(1*%pi/6),   _f_dome(0)*sin(1*%pi/6),   0;             // Node 111
         _f_dome(0)*cos(2*%pi/6),   _f_dome(0)*sin(2*%pi/6),   0;             // Node 112
         _f_dome(0)*cos(3*%pi/6),   _f_dome(0)*sin(3*%pi/6),   0;             // Node 113
         _f_dome(0)*cos(4*%pi/6),   _f_dome(0)*sin(4*%pi/6),   0;             // Node 114
         _f_dome(0)*cos(5*%pi/6),   _f_dome(0)*sin(5*%pi/6),   0;             // Node 115
         _f_dome(0)*cos(6*%pi/6),   _f_dome(0)*sin(6*%pi/6),   0;             // Node 116
         _f_dome(0)*cos(7*%pi/6),   _f_dome(0)*sin(7*%pi/6),   0;             // Node 117
         _f_dome(0)*cos(8*%pi/6),   _f_dome(0)*sin(8*%pi/6),   0;             // Node 118
         _f_dome(0)*cos(9*%pi/6),   _f_dome(0)*sin(9*%pi/6),   0;             // Node 119
         _f_dome(0)*cos(10*%pi/6),  _f_dome(0)*sin(10*%pi/6),  0;             // Node 120
         _f_dome(0)*cos(11*%pi/6),  _f_dome(0)*sin(11*%pi/6),  0              // Node 121
        ];

    t = [ 1,     2;   // Floor 10
          1,     3;
          1,     4;
          1,     5;
          1,     6;
          1,     7;
          1,     8;
          1,     9;
          1,    10;
          1,    11;
          1,    12;
          1,    13;
          2+0,   3+0; // Floor 9
          3+0,   4+0;
          4+0,   5+0;
          5+0,   6+0;
          6+0,   7+0;
          7+0,   8+0;
          8+0,   9+0;
          9+0,  10+0;
         10+0,  11+0;
         11+0,  12+0;
         12+0,  13+0;
          2+0,  14+0;
          3+0,  15+0;
          4+0,  16+0;
          5+0,  17+0;
          6+0,  18+0;
          7+0,  19+0;
          8+0,  20+0;
          9+0,  21+0;
         10+0,  22+0;
         11+0,  23+0;
         12+0,  24+0;
         13+0,  25+0;
          2+0,  15+0;
          3+0,  14+0;
          3+0,  16+0;
          4+0,  15+0;
          4+0,  17+0;
          5+0,  16+0;
          5+0,  18+0;
          6+0,  17+0;
          6+0,  19+0;
          7+0,  18+0;
          7+0,  20+0;
          8+0,  19+0;
          8+0,  21+0;
          9+0,  20+0;
          9+0,  22+0;
         10+0,  21+0;
         10+0,  23+0;
         11+0,  22+0;
         11+0,  24+0;
         12+0,  23+0;
         12+0,  25+0;
         13+0,  24+0;
         13+0,  14+0;
          2+0,  25+0;
          2+12,  3+12; // Floor 8
          3+12,  4+12;
          4+12,  5+12;
          5+12,  6+12;
          6+12,  7+12;
          7+12,  8+12;
          8+12,  9+12;
          9+12, 10+12;
         10+12, 11+12;
         11+12, 12+12;
         12+12, 13+12;
          2+12, 14+12;
          3+12, 15+12;
          4+12, 16+12;
          5+12, 17+12;
          6+12, 18+12;
          7+12, 19+12;
          8+12, 20+12;
          9+12, 21+12;
         10+12, 22+12;
         11+12, 23+12;
         12+12, 24+12;
         13+12, 25+12;
          2+12, 15+12;
          3+12, 14+12;
          3+12, 16+12;
          4+12, 15+12;
          4+12, 17+12;
          5+12, 16+12;
          5+12, 18+12;
          6+12, 17+12;
          6+12, 19+12;
          7+12, 18+12;
          7+12, 20+12;
          8+12, 19+12;
          8+12, 21+12;
          9+12, 20+12;
          9+12, 22+12;
         10+12, 21+12;
         10+12, 23+12;
         11+12, 22+12;
         11+12, 24+12;
         12+12, 23+12;
         12+12, 25+12;
         13+12, 24+12;
         13+12, 14+12;
          2+12, 25+12;
          2+24,  3+24; // Floor 7
          3+24,  4+24;
          4+24,  5+24;
          5+24,  6+24;
          6+24,  7+24;
          7+24,  8+24;
          8+24,  9+24;
          9+24, 10+24;
         10+24, 11+24;
         11+24, 12+24;
         12+24, 13+24;
          2+24, 14+24;
          3+24, 15+24;
          4+24, 16+24;
          5+24, 17+24;
          6+24, 18+24;
          7+24, 19+24;
          8+24, 20+24;
          9+24, 21+24;
         10+24, 22+24;
         11+24, 23+24;
         12+24, 24+24;
         13+24, 25+24;
          2+24, 15+24;
          3+24, 14+24;
          3+24, 16+24;
          4+24, 15+24;
          4+24, 17+24;
          5+24, 16+24;
          5+24, 18+24;
          6+24, 17+24;
          6+24, 19+24;
          7+24, 18+24;
          7+24, 20+24;
          8+24, 19+24;
          8+24, 21+24;
          9+24, 20+24;
          9+24, 22+24;
         10+24, 21+24;
         10+24, 23+24;
         11+24, 22+24;
         11+24, 24+24;
         12+24, 23+24;
         12+24, 25+24;
         13+24, 24+24;
         13+24, 14+24;
          2+24, 25+24;
          2+36,  3+36; // Floor 6
          3+36,  4+36;
          4+36,  5+36;
          5+36,  6+36;
          6+36,  7+36;
          7+36,  8+36;
          8+36,  9+36;
          9+36, 10+36;
         10+36, 11+36;
         11+36, 12+36;
         12+36, 13+36;
          2+36, 14+36;
          3+36, 15+36;
          4+36, 16+36;
          5+36, 17+36;
          6+36, 18+36;
          7+36, 19+36;
          8+36, 20+36;
          9+36, 21+36;
         10+36, 22+36;
         11+36, 23+36;
         12+36, 24+36;
         13+36, 25+36;
          2+36, 15+36;
          3+36, 14+36;
          3+36, 16+36;
          4+36, 15+36;
          4+36, 17+36;
          5+36, 16+36;
          5+36, 18+36;
          6+36, 17+36;
          6+36, 19+36;
          7+36, 18+36;
          7+36, 20+36;
          8+36, 19+36;
          8+36, 21+36;
          9+36, 20+36;
          9+36, 22+36;
         10+36, 21+36;
         10+36, 23+36;
         11+36, 22+36;
         11+36, 24+36;
         12+36, 23+36;
         12+36, 25+36;
         13+36, 24+36;
         13+36, 14+36;
          2+36, 25+36;
          2+48,  3+48; // Floor 5
          3+48,  4+48;
          4+48,  5+48;
          5+48,  6+48;
          6+48,  7+48;
          7+48,  8+48;
          8+48,  9+48;
          9+48, 10+48;
         10+48, 11+48;
         11+48, 12+48;
         12+48, 13+48;
          2+48, 14+48;
          3+48, 15+48;
          4+48, 16+48;
          5+48, 17+48;
          6+48, 18+48;
          7+48, 19+48;
          8+48, 20+48;
          9+48, 21+48;
         10+48, 22+48;
         11+48, 23+48;
         12+48, 24+48;
         13+48, 25+48;
          2+48, 15+48;
          3+48, 14+48;
          3+48, 16+48;
          4+48, 15+48;
          4+48, 17+48;
          5+48, 16+48;
          5+48, 18+48;
          6+48, 17+48;
          6+48, 19+48;
          7+48, 18+48;
          7+48, 20+48;
          8+48, 19+48;
          8+48, 21+48;
          9+48, 20+48;
          9+48, 22+48;
         10+48, 21+48;
         10+48, 23+48;
         11+48, 22+48;
         11+48, 24+48;
         12+48, 23+48;
         12+48, 25+48;
         13+48, 24+48;
         13+48, 14+48;
          2+48, 25+48;
          2+60,  3+60; // Floor 3
          3+60,  4+60;
          4+60,  5+60;
          5+60,  6+60;
          6+60,  7+60;
          7+60,  8+60;
          8+60,  9+60;
          9+60, 10+60;
         10+60, 11+60;
         11+60, 12+60;
         12+60, 13+60;
          2+60, 14+60;
          3+60, 15+60;
          4+60, 16+60;
          5+60, 17+60;
          6+60, 18+60;
          7+60, 19+60;
          8+60, 20+60;
          9+60, 21+60;
         10+60, 22+60;
         11+60, 23+60;
         12+60, 24+60;
         13+60, 25+60;
          2+60, 15+60;
          3+60, 14+60;
          3+60, 16+60;
          4+60, 15+60;
          4+60, 17+60;
          5+60, 16+60;
          5+60, 18+60;
          6+60, 17+60;
          6+60, 19+60;
          7+60, 18+60;
          7+60, 20+60;
          8+60, 19+60;
          8+60, 21+60;
          9+60, 20+60;
          9+60, 22+60;
         10+60, 21+60;
         10+60, 23+60;
         11+60, 22+60;
         11+60, 24+60;
         12+60, 23+60;
         12+60, 25+60;
         13+60, 24+60;
         13+60, 14+60;
          2+60, 25+60;
          2+72,  3+72; // Floor 2
          3+72,  4+72;
          4+72,  5+72;
          5+72,  6+72;
          6+72,  7+72;
          7+72,  8+72;
          8+72,  9+72;
          9+72, 10+72;
         10+72, 11+72;
         11+72, 12+72;
         12+72, 13+72;
          2+72, 14+72;
          3+72, 15+72;
          4+72, 16+72;
          5+72, 17+72;
          6+72, 18+72;
          7+72, 19+72;
          8+72, 20+72;
          9+72, 21+72;
         10+72, 22+72;
         11+72, 23+72;
         12+72, 24+72;
         13+72, 25+72;
          2+72, 15+72;
          3+72, 14+72;
          3+72, 16+72;
          4+72, 15+72;
          4+72, 17+72;
          5+72, 16+72;
          5+72, 18+72;
          6+72, 17+72;
          6+72, 19+72;
          7+72, 18+72;
          7+72, 20+72;
          8+72, 19+72;
          8+72, 21+72;
          9+72, 20+72;
          9+72, 22+72;
         10+72, 21+72;
         10+72, 23+72;
         11+72, 22+72;
         11+72, 24+72;
         12+72, 23+72;
         12+72, 25+72;
         13+72, 24+72;
         13+72, 14+72;
          2+72, 25+72;
          2+84,  3+84; // Floor 1
          3+84,  4+84;
          4+84,  5+84;
          5+84,  6+84;
          6+84,  7+84;
          7+84,  8+84;
          8+84,  9+84;
          9+84, 10+84;
         10+84, 11+84;
         11+84, 12+84;
         12+84, 13+84;
          2+84, 14+84;
          3+84, 15+84;
          4+84, 16+84;
          5+84, 17+84;
          6+84, 18+84;
          7+84, 19+84;
          8+84, 20+84;
          9+84, 21+84;
         10+84, 22+84;
         11+84, 23+84;
         12+84, 24+84;
         13+84, 25+84;
          2+84, 15+84;
          3+84, 14+84;
          3+84, 16+84;
          4+84, 15+84;
          4+84, 17+84;
          5+84, 16+84;
          5+84, 18+84;
          6+84, 17+84;
          6+84, 19+84;
          7+84, 18+84;
          7+84, 20+84;
          8+84, 19+84;
          8+84, 21+84;
          9+84, 20+84;
          9+84, 22+84;
         10+84, 21+84;
         10+84, 23+84;
         11+84, 22+84;
         11+84, 24+84;
         12+84, 23+84;
         12+84, 25+84;
         13+84, 24+84;
         13+84, 14+84;
          2+84, 25+84;
        ];

    e = [
         [1, 1, 1] .* localise3d(110)  // We localise the position of node 110 in the matrix and we fix the dergee of freedom of this node.
         [1, 1, 1] .* localise3d(111)  // We localise the position of node 111 in the matrix and we fix the dergee of freedom of this node.
         [1, 1, 1] .* localise3d(112)  // We localise the position of node 112 in the matrix and we fix the dergee of freedom of this node.
         [1, 1, 1] .* localise3d(113)  // We localise the position of node 113 in the matrix and we fix the dergee of freedom of this node.
         [1, 1, 1] .* localise3d(114)  // We localise the position of node 114 in the matrix and we fix the dergee of freedom of this node.
         [1, 1, 1] .* localise3d(115)  // We localise the position of node 115 in the matrix and we fix the dergee of freedom of this node.
         [1, 1, 1] .* localise3d(116)  // We localise the position of node 116 in the matrix and we fix the dergee of freedom of this node.
         [1, 1, 1] .* localise3d(117)  // We localise the position of node 117 in the matrix and we fix the dergee of freedom of this node.
         [1, 1, 1] .* localise3d(118)  // We localise the position of node 118 in the matrix and we fix the dergee of freedom of this node.
         [1, 1, 1] .* localise3d(119)  // We localise the position of node 119 in the matrix and we fix the dergee of freedom of this node.
         [1, 1, 1] .* localise3d(120)  // We localise the position of node 120 in the matrix and we fix the dergee of freedom of this node.
         [1, 1, 1] .* localise3d(121)  // We localise the position of node 121 in the matrix and we fix the dergee of freedom of this node.
        ]; 

    A   = ones(1,size(t,1)) * 25e-4; // sections of the elements
    E   = ones(1,size(t,1)) * 210e9; // elasticity module of the elements
    rho = ones(1,size(t,1)) * 7.8e3; // volumic mass

    F = [0; 0; -3.4e5 // node 1
         0; 0; -3.4e5 // node 2
         0; 0; -3.4e5 // node 3
         0; 0; -3.4e5 // node 4
         0; 0; -3.4e5 // node 5
         0; 0; -3.4e5 // node 6
         0; 0; -3.4e5 // node 7
         0; 0; -3.4e5 // node 8
         0; 0; -3.4e5 // node 9
         0; 0; -3.4e5 // node 10
         0; 0; -3.4e5 // node 11
         0; 0; -3.4e5 // node 12
         0; 0; -3.4e5 // node 13
         0; 0; -3.4e5 // node 14
         0; 0; -3.4e5 // node 15
         0; 0; -3.4e5 // node 16
         0; 0; -3.4e5 // node 17
         0; 0; -3.4e5 // node 18
         0; 0; -3.4e5 // node 19
         0; 0; -3.4e5 // node 20
         0; 0; -3.4e5 // node 21
         0; 0; -3.4e5 // node 22
         0; 0; -3.4e5 // node 23
         0; 0; -3.4e5 // node 24
         0; 0; -3.4e5 // node 25
         0; 0; -3.4e5 // node 26
         0; 0; -3.4e5 // node 27
         0; 0; -3.4e5 // node 28
         0; 0; -3.4e5 // node 29
         0; 0; -3.4e5 // node 30
         0; 0; -3.4e5 // node 31
         0; 0; -3.4e5 // node 32
         0; 0; -3.4e5 // node 33
         0; 0; -3.4e5 // node 34
         0; 0; -3.4e5 // node 35
         0; 0; -3.4e5 // node 36
         0; 0; -3.4e5 // node 37
         0; 0; -3.4e5 // node 38
         0; 0; -3.4e5 // node 39
         0; 0; -3.4e5 // node 40
         0; 0; -3.4e5 // node 41
         0; 0; -3.4e5 // node 42
         0; 0; -3.4e5 // node 43
         0; 0; -3.4e5 // node 44
         0; 0; -3.4e5 // node 45
         0; 0; -3.4e5 // node 46
         0; 0; -3.4e5 // node 47
         0; 0; -3.4e5 // node 48
         0; 0; -3.4e5 // node 49
         0; 0; -3.4e5 // node 50
         0; 0; -3.4e5 // node 51
         0; 0; -3.4e5 // node 52
         0; 0; -3.4e5 // node 53
         0; 0; -3.4e5 // node 54
         0; 0; -3.4e5 // node 55
         0; 0; -3.4e5 // node 56
         0; 0; -3.4e5 // node 57
         0; 0; -3.4e5 // node 58
         0; 0; -3.4e5 // node 59
         0; 0; -3.4e5 // node 60
         0; 0; -3.4e5 // node 61
         0; 0; -3.4e5 // node 62
         0; 0; -3.4e5 // node 63
         0; 0; -3.4e5 // node 64
         0; 0; -3.4e5 // node 65
         0; 0; -3.4e5 // node 66
         0; 0; -3.4e5 // node 67
         0; 0; -3.4e5 // node 68
         0; 0; -3.4e5 // node 69
         0; 0; -3.4e5 // node 70
         0; 0; -3.4e5 // node 71
         0; 0; -3.4e5 // node 72
         0; 0; -3.4e5 // node 73
         0; 0; -3.4e5 // node 74
         0; 0; -3.4e5 // node 75
         0; 0; -3.4e5 // node 76
         0; 0; -3.4e5 // node 77
         0; 0; -3.4e5 // node 78
         0; 0; -3.4e5 // node 79
         0; 0; -3.4e5 // node 80
         0; 0; -3.4e5 // node 81
         0; 0; -3.4e5 // node 82
         0; 0; -3.4e5 // node 83
         0; 0; -3.4e5 // node 84
         0; 0; -3.4e5 // node 85
         0; 0; -3.4e5 // node 86
         0; 0; -3.4e5 // node 87
         0; 0; -3.4e5 // node 88
         0; 0; -3.4e5 // node 89
         0; 0; -3.4e5 // node 90
         0; 0; -3.4e5 // node 91
         0; 0; -3.4e5 // node 92
         0; 0; -3.4e5 // node 93
         0; 0; -3.4e5 // node 94
         0; 0; -3.4e5 // node 95
         0; 0; -3.4e5 // node 96
         0; 0; -3.4e5 // node 97
         0; 0; -3.4e5 // node 98
         0; 0; -3.4e5 // node 99
         0; 0; -3.4e5 // node 100
         0; 0; -3.4e5 // node 101
         0; 0; -3.4e5 // node 102
         0; 0; -3.4e5 // node 103
         0; 0; -3.4e5 // node 104
         0; 0; -3.4e5 // node 105
         0; 0; -3.4e5 // node 106
         0; 0; -3.4e5 // node 107
         0; 0; -3.4e5 // node 108
         0; 0; -3.4e5 // node 109
         0; 0; 0 // node 110
         0; 0; 0 // node 111
         0; 0; 0 // node 112
         0; 0; 0 // node 113
         0; 0; 0 // node 114
         0; 0; 0 // node 115
         0; 0; 0 // node 116
         0; 0; 0 // node 117
         0; 0; 0 // node 118
         0; 0; 0 // node 119
         0; 0; 0 // node 120
         0; 0; 0 // node 121
        ];
  else
    error('build_fem_test: error - wrong instance name\n');
end
endfunction  
