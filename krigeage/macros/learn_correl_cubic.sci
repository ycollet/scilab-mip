function [parm_sol, err] = learn_correl_cubic(Z)
deff('e=G(p,z)','e = correl_cubic(z(2),p)-z(1);');
[parm_sol,err]=datafit(G,Z,'b',[-%inf, -%inf],[%inf, %inf],[rand(1,1),2*rand(1,1)-1],algo='qn');
printf('Correlation - cubic : parameters a=%f\n',parm_sol(1));
endfunction
