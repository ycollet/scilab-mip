function Result = comp_o_opti_crit(M_doe, model)
  H      = build_regression_matrix(M_doe,model);
  M_Aux  = (H'*H);
  Result = -(sum(abs(M_Aux)) - sum(abs(diag(M_Aux))));
endfunction
