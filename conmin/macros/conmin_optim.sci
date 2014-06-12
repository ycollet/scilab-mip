function [x_opt,f_opt,df_opt,g_opt,dg_opt,ic_res] = conmin_optim(x0,conmin_optim_f,conmin_optim_g,ncon,upper,lower,ItMX,Params)

[nargout,nargin] = argn();

if ~isdef('Params','local') then
  Params = init_param();
end

x_opt  = [];
f_opt  = [];
df_opt = [];
g_opt  = [];
dg_opt = [];
ic_res = [];

[x_opt,f_opt,df_opt,g_opt,dg_opt,ic_res] = conmin_internal(x0, conmin_optim_f, conmin_optim_g, ncon, upper, lower, ItMX, Params);
endfunction

