function Res = compute_tsp_fobj(TSP_Matrix, TSP_Var)
Res = [];
for i=1:size(TSP_Var,1)-1
  Res = Res + TSP_Matrix(TSP_Var(i),TSP_Var(i+1));
end
Res = Res + TSP_Matrix(TSP_Var(i+1),TSP_Var(1));
endfunction
