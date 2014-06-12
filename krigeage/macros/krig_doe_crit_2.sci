function Res = krig_doe_crit_2(Model,x,fmin)
// We compute y_estim + (fmin-y_estim)*Phi(fmin') + y_var*phi(fmin') with fmin' = (fmin-y_estim)/y_var
[Mean,Var] = computeKrig(Model,x);
if Var>0 then
  [P,Q] = cdfnor('PQ',(fmin-Mean)/Var,0,1);
  pdfnor = 1/sqrt(2*%pi)*exp(-((fmin-Mean)/Var)^2);
  Res = Mean + (fmin - Mean)*P + Var*pdfnor;
else
  Res = 0;
end
endfunction
