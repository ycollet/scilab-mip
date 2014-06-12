function [parm_sol, err] = learn_correl_exp(Z)
deff('e=G(p,z)','e=correl_exp(z(2),p)-z(1)');
[parm_sol,err]=datafit(G,Z,'b',[-%inf, -%inf],[%inf, %inf],[2*rand(1,1)-1, 2*rand(1,1)-1],algo='gc'); //qn avant
printf('Correlation - exp : parameters a=%f b=%f\n',parm_sol(1),parm_sol(2));
endfunction
