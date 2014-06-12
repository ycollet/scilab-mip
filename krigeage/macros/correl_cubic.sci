function res=correl_cubic(x,p);
// p(1) > 0
aux = p(1)*abs(x);
Index = find(aux>=1)
aux(Index) = 1;
res = p(2)*(1 - 3*aux.^2 + 2*aux.^3);
endfunction
