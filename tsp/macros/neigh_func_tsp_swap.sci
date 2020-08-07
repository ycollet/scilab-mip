function x_neigh = neigh_func_tsp_swap(x_current, T, param)
if ~isdef('param','local') then
  param = [];
end

if is_param(param,'typeofmove') then
  sa_min_delta = get_param(param,'swap');
end

Index1 = ceil(rand(1,1)*(length(x_current)-1)) + 1;
Index2 = ceil(rand(1,1)*(length(x_current)-1)) + 1;

x_neigh = x_current;

x_neigh([Index1 Index2]) = x_neigh([Index2 Index1]);
endfunction
