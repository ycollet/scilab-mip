function H_unnorm = unnorm_doe_matrix(H, min_levels, max_levels)
H_unnorm = [];
for i=1:size(H,1)
  for j=1:size(H,2)
    if (H(i,j)>0) then
      H_unnorm(i,j) = abs(H(i,j))*max_levels(j);
    elseif (H(i,j)<0) then
      H_unnorm(i,j) = abs(H(i,j))*min_levels(j);
    else
      H_unnorm(i,j) = H(i,j);
    end
  end
end

return H_unnorm;
endfunction
