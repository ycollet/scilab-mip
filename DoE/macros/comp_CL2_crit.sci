function Result = comp_CL2_crit(Data)
// Computation of the CL2 Discrepancy
Result = (13/12)^2;
n = size(Data,1);
m = size(Data,2);
for i=1:n
  Aux1 = 1;
  for k=1:m
    Aux1 = Aux1 * (1+0.5*abs(Data(i,k) - 0.5) - 0.5*abs(Data(i,k)-0.5)^2);
  end
  Result = Result - (2/n)*Aux1
  for j=1:n
    for k=1:m
      Aux1 = Aux1 * (1 + 0.5*abs(Data(i,k) - 0.5) + 0.5*abs(Data(j,k) - 0.5) - 0.5*abs(Data(i,k) - Data(j,k)));
    end
  Result = Result + 1/n^2*Aux1;
  end
end
Result = - Result;
endfunction
