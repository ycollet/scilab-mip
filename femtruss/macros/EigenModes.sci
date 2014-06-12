function [T,Phi] = EigenModes(K,M,nmodes)
// [T,Phi] = EigenModes(K,M,nmodes)
// T      : periodes propres 
// Phi    : modes propes
// nmodes : nombre de modes
  
[n,n]  = size(M);
nmodes = min([nmodes,n]);

if (typeof(M)=='sparse') & (typeof(K)=='sparse') then
  [al, be, Phi] = spec(full(M), full(K));
elseif (typeof(M)=='sparse') & (typeof(K)~='sparse') then
  [al, be, Phi] = spec(full(M), K);
elseif (typeof(M)~='sparse') & (typeof(K)=='sparse') then
  [al, be, Phi] = spec(M, full(K));
else
  [al, be, Phi] = spec(M, K);
end

Index = find(be==0);
if ~isempty(Index) then
  if Index(1)>1 then
    nmodes = Index(1)-1;
    al = al(1:nmodes);
    be = be(1:nmodes);
  else
    error('No eigen modes');
  end
end
Omega = al(1:nmodes) ./ be(1:nmodes);
Phi   = Phi(:,1:nmodes);
T     = 2 * %pi * sqrt(Omega); // YC diag(Omega) avant
Index = find(eval(string(Omega))<eval(string(1E-8))); // YC: bug avec cette comparaison - seul solution pour que ca marche
//Index = find(Omega<=1e-8); // YC: bug avec cette comparaison - seul solution pour que ca marche
T(Index) = [];
Phi(:,Index) =[];
endfunction
