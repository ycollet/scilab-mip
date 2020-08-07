function x_neigh = neigh_func_tsp_2opt(x_current, T, param)
if ~isdef('param','local') then
  param = [];
end

if is_param(param,'typeofmove') then
  sa_min_delta = get_param(param,'swap');
end

Index1a = ceil(rand(1,1)*(length(x_current)-1)) + 1;
Index1b = Index1a+1; 
if (Index1b>length(x_current)) then 
  Index1b = 1; 
end
Index2a = ceil(rand(1,1)*(length(x_current)-1)) + 1;
while Index2a==Index1a
  Index2a = ceil(rand(1,1)*(length(x_current)-1)) + 1;
end
Index2b = Index2a+1;
if (Index2b>length(x_current)) then 
  Index2b = 1; 
end
if Index2b==Index1a then 
  Index2b = Index2a-1; 
end
Index2b = max(Index2b,1);

x_neigh = x_current;

x_neigh([Index1b Index1a Index2b Index2a]) = x_neigh([Index1b Index2b Index1a Index2a]);
endfunction
