function [f_opt, x_opt, status, meth_name] = nlopt_unc(x0, f, lower, upper, ItMX, params);
if is_param(params,'nb_ceq') then
  [params, err] = set_param(params,'nb_ceq',0);
else
  [params, err] = add_param(params,'nb_ceq',0);
end

if is_param(params,'nb_cineq') then
  [params, err] = set_param(params,'nb_cineq',0);
else
  [params, err] = add_param(params,'nb_cineq',0);
end

[f_opt, x_opt, status, meth_name] = nlopt(x0, f, [], [], lower, upper, ItMX, params);
endfunction
