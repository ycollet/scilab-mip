function Result = comp_ssd_max_khi2_crit(M_doe, Model)
n = size(M_doe,2);
AuxMat = ones(n,n)*(-%inf);
for i=1:n-1
  for j=i+1:n
    Aux = comp_ssd_khi2ij_crit(M_doe, i, j);
    AuxMat(i,j) = Aux;
    AuxMat(j,i) = Aux;
  end
end
Result = - max(AuxMat);
endfunction

