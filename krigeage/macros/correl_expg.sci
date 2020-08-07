function res=correl_expg(x,p);
if p(2)==0 then p(2) = %eps; end
Index = find(x==0);
x(Index) = %eps;

res=p(1)*exp(p(2)*x^p(3));
endfunction
