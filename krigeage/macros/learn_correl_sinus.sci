function [parm_sol, err] = learn_correl_sinus(Z)
deff('e=G(p,z)','e=correl_sinus(z(2),p)-z(1)');
[parm_sol,err]=datafit(G,Z,2*rand(4,1)-1,algo='qn');
printf('Correlation - sinus : parameters a=%f b=%f c=%f d=%f\n',parm_sol(1),parm_sol(2),parm_sol(3),parm_sol(4));
endfunction

