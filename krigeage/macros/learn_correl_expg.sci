function [parm_sol, err] = learn_correl_expg(Z)
deff('e=G(p,z)','e=correl_expg(z(2),p)-z(1)');
[parm_sol,err]=datafit(G,Z,'b',[-%inf,-%inf,0.001], [%inf %inf %inf], [2*rand(1,1)-1 , 2*rand(1,1)-1 , rand(1,1)],algo='qn');
printf('Correlation - expg : parameters a=%f b=%f c=%f\n',parm_sol(1),parm_sol(2),parm_sol(3));
endfunction
