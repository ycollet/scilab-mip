function H = doe_star(nb_var)
H = eye(nb_var, nb_var);
H = [H -H]';
return H;
endfunction
