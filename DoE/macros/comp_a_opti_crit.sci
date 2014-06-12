function Result = comp_a_opti_crit(M_doe, model)
  H      = build_regression_matrix(M_doe,model);
  Result = - trace((H'*H)^(-1));
endfunction
