function [parm_sol, err] = learn_correl_lin(Z)
deff('e=G(p,z)','e=correl_lin(z(2),p)-z(1)');
//[parm_sol,err]=datafit(G,Z,'b',[-%inf, -%inf],[%inf, %inf],[100*rand(1,1), 2*rand(1,1)-1],algo='gc','ar',10000,10000,1e-12,1e-12);
[parm_sol,err]=datafit(G,Z,'b',[-%inf, -%inf],[%inf, %inf],[2*rand(1,1)-1, 2*rand(1,1)-1],algo='gc');
printf('Correlation - lin = max(0, b*(1-a*x)) : parameters a=%f b =%f\n',parm_sol(1),parm_sol(2));
endfunction
