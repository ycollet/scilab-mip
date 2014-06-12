function Result = comp_ssd_r_value_crit(M_doe, Model)
  n    = size(M_doe,2);
  CMat = zeros(n,n);
  for i=1:n
    for j=1:n
      CMat(i,j) = corr(M_doe(:,i),M_doe(:,j),1);
    end
  end
  Result = -(norm(CMat));
endfunction
