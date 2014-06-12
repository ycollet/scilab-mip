function H = doe_diff(H1, H2)
  Aux_H1 = string(H1);
  Aux_H2 = string(H2);
  Aux    = intersect(Aux_H1, Aux_H2);
  H      = eval(Aux);
endfunction
