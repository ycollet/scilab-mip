function [parm_sol, err] = learn_correl_gauss(Z)
deff('e=G(p,z)','e=correl_gauss(z(2),p)-z(1)');
[parm_sol,err]=datafit(G,Z,'b',[-%inf, 0],[%inf, %inf],[2*rand(1,1)-1, rand(1,1)],algo='qn');
printf('Correlation - gauss : parameters a=%f b=%f\n',parm_sol(1),parm_sol(2));
endfunction
