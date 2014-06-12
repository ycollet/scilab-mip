function dU = dfemtruss_ana(ffd, U, listofnodes, Log, varargin)
// dU = dfemtruss(ffd,U,listofnodes, Log)
// resolution des systemes d'assemblage de barres bidimensionnels 
// dU          : derive partielle de la solution en deplacements nodaux 
// U           : solution en deplacements nodaux
// ffd         : fichier fonction de donnees du probleme
// listofnodes : a list of nodes where to compute the partial derivative
// Log         : affiche quelques messages informatifs
// A. Seghir, le 06/08/04 modifie le 27/10/04

//
// Il faut que dU = [Ndof x Ndof]
// Les colonnes contiennent dUi et les lignes dxj
// Ensuite, il faut utiliser ces morceaux dans l'expression de la fonction objectif derivee
//

if ~isdef('Log','local') then
  Log = %F;
end
if ~isdef('listofnodes','local') then
  listofnodes = 1:size(p,1);
elseif isempty(listofnodes) then
  listofnodes = 1:size(p,1);
end

if (size(varargin)~=0) then 
  [t,p,e,A,E,rho,F] = ffd(varargin(1));
else
  [t,p,e,A,E,rho,F] = ffd();
end
  
_3D_problem = (size(p,2)==3);

if _3D_problem then
  dU = zeros(3*length(listofnodes),length(U));
  net = size(t,1);
  nnt = 3*size(p,1);
  // Computation of the stiffness matrix
  K = sparse([1 1;nnt nnt],[0 0],[nnt,nnt]);
  for i=1:net
    ti = t(i,:);
    Li = localise3d(ti);
    Ke = truss3dKe(p(ti,:),A(i),E(i));
    K(Li,Li) = K(Li,Li) + Ke;
  end
  K     = DelDOFs(K,e);
  Uprim = DelDOFs(U,e);
  KInv  = inv(K);

  // Computation of the stiffness partial derivatives matrix
  Index = 1;
  for i=1:length(listofnodes)
    dKdxG = sparse([1 1;nnt nnt],[0 0],[nnt,nnt]);
    dKdyG = sparse([1 1;nnt nnt],[0 0],[nnt,nnt]);
    dKdzG = sparse([1 1;nnt nnt],[0 0],[nnt,nnt]);

    // Pour chaque point:
    // * suivant x : il faut rechercher les barres qui commencent par ce point (on l'appelle x1) et construire la matrice dKdx1
    // * suivant x : il faut rechercher les barres qui finissent  par ce point (on l'appelle x2) et construire la matrice dKdx2
    // * suivant y : il faut rechercher les barres qui commencent par ce point (on l'appelle y1) et construire la matrice dKdy1
    // * suivant y : il faut rechercher les barres qui finissent  par ce point (on l'appelle y2) et construire la matrice dKdy2
    // * suivant z : il faut rechercher les barres qui commencent par ce point (on l'appelle z1) et construire la matrice dKdz1
    // * suivant z : il faut rechercher les barres qui finissent  par ce point (on l'appelle z2) et construire la matrice dKdz2
    Index_Begin = find(t(:,1)==listofnodes(i));
    for j=1:length(Index_Begin)
      ti = t(Index_Begin(j),:);
      Li = localise3d(ti);
      [dKdx1,dKdx2,dKdy1,dKdy2,dKdz1,dKdz2] = dtruss3dKe(p(ti,:),A(Index_Begin(j)),E(Index_Begin(j)),'xyz1');
      dKdxG(Li,Li) = dKdxG(Li,Li) + dKdx1;
      dKdyG(Li,Li) = dKdyG(Li,Li) + dKdy1;
      dKdzG(Li,Li) = dKdzG(Li,Li) + dKdz1;
    end

    Index_End = find(t(:,2)==listofnodes(i));
    for j=1:length(Index_End)
      ti = t(Index_End(j),:);
      Li = localise3d(ti);
      [dKdx1,dKdx2,dKdy1,dKdy2,dKdz1,dKdz2] = dtruss3dKe(p(ti,:),A(Index_End(j)),E(Index_End(j)),'xyz2');
      dKdxG(Li,Li) = dKdxG(Li,Li) + dKdx2;
      dKdyG(Li,Li) = dKdyG(Li,Li) + dKdy2;
      dKdzG(Li,Li) = dKdzG(Li,Li) + dKdz2;
    end
  
    // Computation of the partial derivative for the current point
    dKdxGprim = DelDOFs(dKdxG,e);
    dKdyGprim = DelDOFs(dKdyG,e);
    dKdzGprim = DelDOFs(dKdzG,e);

    dU(Index,:) = full(AddDOFsToVect((- KInv*dKdxGprim*Uprim'),e)');
    Index = Index + 1;
    dU(Index,:) = full(AddDOFsToVect((- KInv*dKdyGprim*Uprim'),e)');
    Index = Index + 1;
    dU(Index,:) = full(AddDOFsToVect((- KInv*dKdzGprim*Uprim'),e)');
    Index = Index + 1;
  end
else
  dU = zeros(2*length(listofnodes),length(U));
  net = size(t,1);
  nnt = 2*size(p,1);
  // Computation of the stiffness matrix
  K = sparse([1 1;nnt nnt],[0 0],[nnt,nnt]);
  for i=1:net
    ti = t(i,:);
    Li = localise2d(ti);
    Ke = truss2dKe(p(ti,:),A(i),E(i));
    K(Li,Li) = K(Li,Li) + Ke;
  end

  K     = DelDOFs(K,e);
  Uprim = DelDOFs(U,e);
  KInv  = inv(K);

  // Computation of the stiffness partial derivatives matrix
  Index = 1;
  for i=1:length(listofnodes)
    dKdxG = sparse([1 1;nnt nnt],[0 0],[nnt,nnt]);
    dKdyG = sparse([1 1;nnt nnt],[0 0],[nnt,nnt]);

    // Pour chaque point:
    // * suivant x : il faut rechercher les barres qui commencent par ce point (on l'appelle x1) et construire la matrice dKdx1
    // * suivant x : il faut rechercher les barres qui finissent  par ce point (on l'appelle x2) et construire la matrice dKdx2
    // * suivant y : il faut rechercher les barres qui commencent par ce point (on l'appelle y1) et construire la matrice dKdy1
    // * suivant y : il faut rechercher les barres qui finissent  par ce point (on l'appelle y2) et construire la matrice dKdy2
    Index_Begin = find(t(:,1)==listofnodes(i));
    for j=1:length(Index_Begin)
      ti = t(Index_Begin(j),:);
      Li = localise2d(ti);
      [dKdx1,dKdx2,dKdy1,dKdy2] = dtruss2dKe(p(ti,:),A(Index_Begin(j)),E(Index_Begin(j)),'xy1');
      dKdxG(Li,Li) = dKdxG(Li,Li) + dKdx1;
      dKdyG(Li,Li) = dKdyG(Li,Li) + dKdy1;
    end

    Index_End = find(t(:,2)==listofnodes(i));
    for j=1:length(Index_End)
      ti = t(Index_End(j),:);
      Li = localise2d(ti);
      [dKdx1,dKdx2,dKdy1,dKdy2] = dtruss2dKe(p(ti,:),A(Index_End(j)),E(Index_End(j)),'xy2');
      dKdxG(Li,Li) = dKdxG(Li,Li) + dKdx2;
      dKdyG(Li,Li) = dKdyG(Li,Li) + dKdy2;
    end
  
    // Computation of the partial derivative for the current point
    dKdxGprim = DelDOFs(dKdxG,e);
    dKdyGprim = DelDOFs(dKdyG,e);

    dU(Index,:) = full(AddDOFsToVect((- KInv*dKdxGprim*Uprim'),e)');
    Index = Index + 1;
    dU(Index,:) = full(AddDOFsToVect((- KInv*dKdyGprim*Uprim'),e)');
    Index = Index + 1;
  end
end
endfunction
