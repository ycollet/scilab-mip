function res=correl_lin(x,p);
if p(2)==0 then p(2) = %eps; end
Index = find(x==0);
x(Index) = %eps;

res = p(2)*(1.0 - p(1)*x);
Index = find(res<0);
res(Index) = 0;
endfunction
