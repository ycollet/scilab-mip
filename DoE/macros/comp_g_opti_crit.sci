function Result = comp_g_opti_crit(M_doe, model)
  H      = build_regression_matrix(M_doe,model);
  Result = - max(diag((H'*H)^(-1)));
endfunction
