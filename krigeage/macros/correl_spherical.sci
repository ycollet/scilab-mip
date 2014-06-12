function res=correl_spherical(x,p);
aux    = p(1)*abs(x);
Index = find(aux>=1);
aux(Index) = 1;
res   = p(2)*(1 - 1.5*aux + 0.5*aux.^3);
Index = find(res<0);
res(Index) = 0;
endfunction
