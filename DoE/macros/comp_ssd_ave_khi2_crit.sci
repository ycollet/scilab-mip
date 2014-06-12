function Result = comp_ssd_ave_khi2_crit(M_doe, Model)
n = size(M_doe,2);
Result = 0;
for i=1:n-1
  for j=i+1:n
    Result = Result + comp_ssd_khi2ij_crit(M_doe, i, j);
  end
end
Result = - Result / (n*(n-1)/2);
endfunction

function Result = comp_ssd_khi2ij_crit(M_doe, i, j)
// Get a vector of number of levels
Min_Level_i = min(M_doe(:,i));
Max_Level_i = max(M_doe(:,i));
Min_Level_j = min(M_doe(:,j));
Max_Level_j = max(M_doe(:,j));

S_i = Max_Level_i - Min_Level_i;
S_j = Max_Level_j - Min_Level_j;

n = size(M_doe,1);
Result = 0;

for a=Min_Level_i:Max_Level_i
  for b=Min_Level_j:Max_Level_j
    n_ab = size(vectorfind([M_doe(:,i) M_doe(:,j)],[a b],'r'),2);
    Result = Result + (n_ab - n/(S_i * S_j)^2)/(n/(S_i*S_j));
  end
end
endfunction
