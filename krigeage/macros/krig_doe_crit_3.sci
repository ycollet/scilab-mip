function Res = krig_doe_crit_3(Model,x,fmin)
// We compute y_var
[Mean,Var] = computeKrig(Model,x);
Res = abs(Var);
endfunction
