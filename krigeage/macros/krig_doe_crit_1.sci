function Res = krig_doe_crit_1(Model,x,fmin)
// We compute Phi((fmin-y_estim)/y_var)
[Mean,Var] = computeKrig(Model,x);
[P,Q] = cdfnor('PQ',(fmin-Mean)/Var,0,1);
Res = P;
endfunction
