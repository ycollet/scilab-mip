function [U,P,R,Kref,Mref]= femtruss(ffd, Log, varargin)
// [U,P,T,phi] = trussfem(ffd)
// resolution des systemes d'assemblage de barres bidimensionnels 
// U          : solution en deplacements nodaux 
// P          : forces axiales dans les barres 
// R          : reaction aux noeuds immobilises
// Kref       : la matrice de rigidite (optionnel)
// Mref       : la matrice de masse (optionnel)
// ffd        : fichier fonction de donnees du probleme
// Log        : affiche quelques messages informatifs
// A. Seghir, le 06/08/04 modifie le 27/10/04

if ~isdef('Log','local') then
  Log = %F;
end

if (size(varargin)~=0) then 
  [t,p,e,A,E,rho,F] = ffd(varargin(1));
else
  [t,p,e,A,E,rho,F] = ffd();
end

_3D_problem = (size(p,2)==3);

if _3D_problem then
  [Kref,Mref] = truss3dKM(t,p,A,E,rho);
else
  [Kref,Mref] = truss2dKM(t,p,A,E,rho);
end

K = DelDOFs(Kref,e);
M = DelDOFs(Mref,e);
F = DelDOFs(F,e);

U = full(K \ F); 
U = AddDOFs(U,e); // Ajouter les DDL des noeuds encastre

[P,R] = TrussForces(t,p,A,E,U);

net = size(t,1);
nnt = size(p,1);

if Log then
  // impression des resultats
  printf(' Result of the computation on the truss structure\n');
  printf(' Number of elements : %d\n',net);
  printf(' Number of nodes    : %d\n',nnt);

  printf('Displacements of the nodes :\n');
  printf(' Node\t\t Ux\t\t\t\t Uy \n');
  for i=1:nnt
    if _3D_problem then printf(' %d\t\t\t%+5.5f\t\t\t%+5.5f\t\t\t%+5.5f\n',i,U(localise3d(i)));
    else                printf(' %d\t\t\t%+5.5f\t\t\t%+5.5f\n',i,U(localise2d(i))); 
    end
  end

  printf('\n Contraints in the elements :\n');
  printf(' Element \t\t P\n')
  for i=1:net
    printf(' %d\t\t\t%+1.5f\n',i,P(i));
  end

  printf('Force at the support nodes :\n');
  printf(' Support\t\t Rx\t\t\t\t\t Ry \n');
  for i=1:nnt
    if _3D_problem then
      L = localise3d(i);
      if find(e==L(1))
        printf(' %d\t\t\t%+1.5f\t\t%+1.5f\t\t%+1.5f\n',i,R(L(1)), R(L(2)), R(L(3)));
      end
    else
      L = localise2d(i);
      if find(e==L(1))
        printf(' %d\t\t\t%+1.5f\t\t%+1.5f\n',i,R(L(1)), R(L(2)));
      end
    end
  end
end
endfunction
