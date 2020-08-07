function res = correl_materm(x,p);
// p(1) = sigma > 0
// p(2) = nu > 0
// p(3) = rho > 0

if p(1)<=0 then p(1) = %eps; end
if p(2)<=0 then p(2) = %eps; end
if p(3)<=0 then p(3) = %eps; end

Sigma2b = 0.01;
res = (Sigma2b + p(1));
Index = find(x<0.001);
res(Index) = p(1);

t = 2.0*sqrt(p(2)) * x / p(3);
// Limit behavior of Knu
aux = abs(t);
Index = find(aux==0);
aux(Index) = %eps;

bessellimit = 2.0^(p(2) - 1.0) * gamma(p(2)) * aux .^(-p(2));
Index = find((bessellimit<1e30) & (x>=0.01));
res(Index) = p(1) / (2.0^(p(2) - 1.0) * gamma(p(2))) * t(Index) .^p(2) .* besselk(p(2), t(Index));
endfunction
