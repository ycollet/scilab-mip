function Result = comp_d_opti_crit(M_doe, model)
  H      = build_regression_matrix(M_doe,model);
  [m,e]  = det((H'*H)^(-1));
  Result = abs(m*10^e);
endfunction
