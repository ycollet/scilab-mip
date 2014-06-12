function H = doe_composite(nb_var, alpha)
if ~isdef('alpha','local') then
  alpha = sqrt(2);
end
H1 = doe_factorial(nb_var);
H2 = doe_star(nb_var);
H = doe_union(H1, alpha*H2);
endfunction
