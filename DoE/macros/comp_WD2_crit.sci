function Result = comp_WD2_crit(M_doe, Model)
// Computation of the WD2 Discrepancy (Wrap around)
n = size(M_doe,1);
m = size(M_doe,2);
Result = -(4/3)^m;
for i=1:n
  for j=1:n
    Aux1 = 1;
    for k=1:m
      Aux1 = Aux1 * (1.5 - abs(M_doe(i,k) - M_doe(j,k)) * (1 - abs(M_doe(i,k) - M_doe(j,k))));
    end
  Result = Result + 1/n^2*Aux1;
  end
end
Result = - Result;
endfunction
