function Result = comp_ssd_a_value_crit(M_doe, model)
  H      = build_regression_matrix(M_doe,model);
  n = size(H,2);
  IMat   = H'*H / n;
  Result = -(1 / (sum(spec(IMat).^(-1)) / n));
endfunction
