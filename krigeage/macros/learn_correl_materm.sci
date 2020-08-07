function [parm_sol, err] = learn_correl_materm(Z)
deff('e=G(p,z)','e=correl_materm(z(2),p)-z(1)');
[parm_sol,err]=datafit(G,Z,'b',[0,0,0],[%inf,%inf,%inf], [rand(1,1),rand(1,1),rand(1,1)],algo='qn');
printf('Correlation - exp : parameters sigma=%f nu=%f rho=%f\n', parm_sol(1), parm_sol(2), parm_sol(3));
endfunction
