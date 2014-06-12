function var = var_regression_matrix(H, x, model, sigma)
  // This function computes the variance of the "regression error"
  if (~isdef('sigma','local')) then
    sigma = 1;
  end
  
  if (size(x,1)==1) then
    x = x';
  end
  
  x_mod = build_regression_matrix(x,model);
  var   = sigma^2 * x_mod'*(H'*H)^(-1)*x_mod;
endfunction
