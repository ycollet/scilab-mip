function [U,dU]= dfemtruss_smw(ffd, x, ddl, params)
// computation of the gradient by finite differences + the sherman morrison woodruff formula
// Work done by Ms. Ibrahimcha, Nkamani, Hammad.

if ~isdef('params','local') then
  params = [];
end

U  = [];
dU = [];

[t,p,e,A,E,rho,F] = ffd(x);

_3D_problem = (size(p,2)==3);

if _3D_problem then
  [K,M] = truss3dKM(t,p,A,E,rho);
else
  [K,M] = truss2dKM(t,p,A,E,rho);
end

[delta_h,err] = get_param(params, 'delta_h', 1e-6);

K = DelDOFs(K,e);
M = DelDOFs(M,e);
F = DelDOFs(F,e);

Kinv = inv(K);
U = Kinv * F;
U = AddDOFs(U,e); // Add the degrees of freedom of the fixed nodes

// Perform the finite differences
delta = zeros(x);

dU(:,1) = zeros(U(ddl)');

for i=1:length(x)
  delta(i) = delta_h;
  x_new = x + delta;
  
  [t_new,p_new,e_new,A_new,E_new,rho_new,F_new] = ffd(x_new);
  
  if _3D_problem then
    [K_new,M_new] = truss3dKM(t_new,p_new,A_new,E_new,rho_new);
  else
    [K_new,M_new] = truss2dKM(t_new,p_new,A_new,E_new,rho_new);
  end

  K_new = DelDOFs(K_new,e_new);
  M_new = DelDOFs(M_new,e_new);
  F_new = DelDOFs(F_new,e_new);

  Kinv_new = sherwoodmor(Kinv,K_new - K);

  U_new = Kinv_new * F_new;
  U_new = AddDOFs(U_new,e_new); // Ajouter les DDL des noeuds encastre
    
  dU(:,i) = (U_new(ddl)' - U(ddl)') / delta(i);
  
  delta(i) = 0;
end
endfunction
