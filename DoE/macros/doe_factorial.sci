function H = doe_factorial(nb_var)
// We count in binary and replace each 0 by -1
H        = [];
nb_lines = 2^nb_var;
Aux      = zeros(1,nb_var);

for i=1:nb_lines-1
  Aux(1) = Aux(1) + 1;
  for j=2:nb_var
    if (Aux(j-1)>1) then
      Aux(j-1) = 0;
      Aux(j)   = Aux(j) + 1;
    end
  end
  H = [H' Aux']';
end

// We add a 0 line to the regression matrix
H = [H' zeros(1,nb_var)']'
// We replace each 0 by -1
I = find (H==0);
H(I) = -1;

return H;
endfunction
