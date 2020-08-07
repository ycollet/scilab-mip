function TSP_Matrix_Dist = compute_tsp_dist(TSP_List_Town, Log)
if (~isdef('Log','local')) then
  Log = %F;
end
n = size(TSP_List_Town,1);
TSP_Matrix_Dist = zeros(n,n);
Index = 0;
for i=1:n
  for j=1:n
    Index = Index + 1;
    if (modulo(Index,10)) then
      if (Log) then
        printf('compute_tsp_dist: step %d / %d\n', Index, n*n);
      end
    end
    if (i==j) then
      TSP_Matrix_Dist(i,j) = %inf;
    else
      TSP_Matrix_Dist(i,j) = sqrt((TSP_List_Town(i,2) - TSP_List_Town(j,2))^2 + (TSP_List_Town(i,3) - TSP_List_Town(j,3))^2);
    end
  end
end
endfunction
