function x_init = nm_init(x0,simplex_relsize)
if ~isdef('simplex_relsize','local') then
  simplex_relsize = 0.1*ones(size(x0,1),size(x0,2));
end

// We set the initial simplex
for i=1:length(x0)+1
  x_init(:,i) = x0 + 0.5 * simplex_relsize .* grand(size(x0,1),size(x0,2),'unf',0,1);
end
endfunction

