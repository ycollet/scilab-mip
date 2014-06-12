function [Umod,T_period,Phi]= femmode(ffd, K, M, Log, nbmodes,varargin)
// [U,P,T,phi] = trussfem(ffd)
// resolution des systemes d'assemblage de barres bidimensionnels 
// Umod       : solutions modales aux noeuds
// T_period   : periodes propres de la structure
// phi        : modes propres de la structure
// ffd        : fichier fonction de donnees du probleme
// K          : la matrice de rigidite (optionnel)
// M          : la matrice de masse (optionnel)
// Log        : affiche quelques messages informatifs
// nbmodes    : number of modes to compute
// A. Seghir, le 06/08/04 modifie le 27/10/04

if ~isdef('Log','local') then
  Log = %F;
end

lines(0);

if (size(varargin)~=0) then 
  [t,p,e,A,E,rho,F] = ffd(varargin(1));
else
  [t,p,e,A,E,rho,F] = ffd();
end

_3D_problem = (size(p,2)==3);

K = DelDOFs(K,e);
M = DelDOFs(M,e);
F = DelDOFs(F,e);

net = size(t,1);
nnt = size(p,1);

if ~isdef('nbmodes','local') then
  [m,m]=size(K);
  nbmodes = m;
end

[T_period,Phi] = EigenModes(K,M,nbmodes);
n=size(T_period);
// We solve M.x'' + lambda.M.x' + M.x = F in the fourier domain
// x_sol = (K - M.w^2 + j.w.lambda.M) / F

for i=1:n(1),
  Umod(:,i) = full(((K-M*(2*%pi/T_period(i))^2 + %i*0.02*M*(2*%pi/T_period(i))) \ F)); 
end

if Log then
  printf('\n Eigen period of the structure - number of modes computed: %d\n',nbmodes);
  printf(' mode\t\tT\n')
  for i=1:length(T_period)
    printf(' %d\t\t\t%5.7f\n',i,T_period(i));
  end
end
for i=1:size(Phi,2)
  Phi2(:,i) = AddDOFsToVect(Phi(:,i),e);
end

Phi = Phi2; clear Phi2;
endfunction

