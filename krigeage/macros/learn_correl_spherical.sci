function [parm_sol, err] = learn_correl_spherical(Z)
deff('e=G(p,z)','e=correl_spherical(z(2),p)-z(1)');
[parm_sol,err]=datafit(G,Z,'b',[-%inf, -%inf],[%inf, %inf],[2*rand(1,1)-1, rand(1,1)],algo='qn');
printf('Correlation - spherical : parameters a=%f\n',parm_sol(1));
endfunction
