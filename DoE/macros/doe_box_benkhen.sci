function H = doe_box_benkhen(nb_var, nb_center)
if (nb_var<3) then
  error('doe_box_benkhen: number of variables must be higher than 2');
end
if (~isdef('nb_center','local')) then
  nb_center = 1;
end

// First part, we compute a factorial doe with 2 parameters
// We count in binary and replace each 0 by -1
H_fact   = doe_factorial(2);
// Now, we populate the real doe with this doe

// On fait un plan factoriel sur chaque couple de dimensions.
// - Donc, on créé un plan factoriel à deux facteurs
// - on fait deux boucles

Index = 0;
nb_lines = nb_var*size(H_fact,1);
H        = zeros(nb_lines,nb_var);
for i=1:nb_var-1
  for j=i+1:nb_var
    Index = Index + 1;
    H(max([1 ((Index-1)*size(H_fact,1)+1)]):(Index)*size(H_fact,1),i) = H_fact(:,1);
    H(max([1 ((Index-1)*size(H_fact,1)+1)]):(Index)*size(H_fact,1),j) = H_fact(:,2);
  end
end

H = [H' doe_repeat_center(nb_var, nb_center)']';

return H;
endfunction
