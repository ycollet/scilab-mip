function wrong_struct = rank_efficient_struct(K,p,e)
if typeof(K)=='sparse' then
  rank_matr = rank(full(K));
else
  rank_matr = rank(K);
end

wrong_struct = (rank_matr==(prod(size(p))-length(e)));
endfunction
