function res=correl_spline(x,p)
aux    = p(1)*abs(x);
aux1   = aux;
Index = find((aux>=0)&(aux<=0.2)); 
aux1(Index) = p(2)*(1 - 15*aux(Index)^2 + 30*aux(Index)^3);
Index = find((aux>0.2)&(aux<1)); 
aux1(Index) = p(2)*(1.25*(1-aux(Index))^3);
Index = find(aux>=1); 
aux1(Index) = 0;
res    = aux1;
endfunction
