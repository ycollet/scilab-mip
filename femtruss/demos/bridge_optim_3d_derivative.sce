// Work done by Ms. Ibrahimcha, Nkamani, Hammad.

lines(0);
warning_old = warning('query');

warning('off');

Order = 4; // 1 2 or 4 for derivative
Debug = %T;
delta_h = 1e-6;

femtruss_path = get_absolute_file_path('bridge_optim_3d_derivative.sce');

//exec(femtruss_path + 'build_long_bridge_3d_bad.sce');
exec(femtruss_path + 'build_long_bridge_3d_good.sce');
  
function y = fobj_truss(x)
[t,p,e,A,E,rho,F] = bridge_optim_3d(x);
[U,P,R]= femtruss(bridge_optim_3d, %F, x);

// First objective: minimize the deformation at nodes of the deck with respect to y (coordinate 2)

// The deck of the bridge
Pos_deck  = localise3d(IndexNodeVarInf);

y = sqrt(sum(U(Pos_deck(2:3:$)).^2));
endfunction


function dy=dfobj_truss_smw(x)
Pos_deck  = localise3d(IndexNodeVarInf);
Pos_free  = localise3d(IndexNodeVarSup);

params = init_param();
[delta_h,err] = add_param(params, 'delta_h', delta_h);

[U,dU] = dfemtruss_smw(bridge_optim_3d,x,Pos_deck(2:3:$),params);

dy = [];
for i=1:length(Pos_free)
  dy(i) = 2*U(Pos_deck(2:3:$)) * dU(:,i);
end
endfunction

function dy=dfobj_truss_diff(x)
dy = derivative(fobj_truss,x,delta_h,Order)';
endfunction

function [y, dy, ind] = optim_fobj_truss(x, ind)
y  = fobj_truss(x);
dy1 = dfobj_truss_smw(x);
dy2 = dfobj_truss_diff(x);

if Debug then
  printf('y = %f\n', y);
  printf('norm(dy1 - dy2) = %f (percent)\n', 100*norm(dy1 - dy2)/max(%eps,norm(dy2)));
  printf('norm(dy1) = %f - abs(dy1): max = %f min = %f\n', norm(dy1), max(abs(dy1)), min(abs(dy1)));
  printf('norm(dy2) = %f - abs(dy2): max = %f min = %f\n', norm(dy2), max(abs(dy2)), min(abs(dy2)));
end

dy = dy1;
endfunction

MaxEvalFunc = 40;
Algorithm   = 'gc'; // 'qn', 'gc', 'nd' -> Ne marche qu'avec 'qn' (quasi-newton). Pour les autres, on obtient rapidement une structure mal
                    // conditionnee
Tol = 1e-12;

////////////////////////////
// Start the Optimization //
////////////////////////////

printf('optimization starting, be patient ... \n\n');

tic();
//[f_opt, x_opt] = optim(optim_fobj_truss, 'b',lower, upper, x0, 'gc','ar',1000,1000,1e-6,1e-6);
//[f_opt, x_opt] = optim(optim_fobj_truss, 'b',lower, upper, x0, 'gc','ar',1000,1000,1e-12,1e-12);

// Be careful: the 'gc' algorithm doesn't seems to accept bounds constraint

[f_opt, x_opt] = optim(optim_fobj_truss, x0, Algorithm,'ar',MaxEvalFunc,MaxEvalFunc,Tol,Tol);

printf('initial solution:'); disp(x0');
printf('initial objective function value = %f\n',fobj_truss(x0));
printf('Final solution:'); disp(x_opt');
printf('Final objective function value = %f\n',fobj_truss(x_opt));
  
t_end=toc();
printf('computation time = %f\n' , t_end);

////////////////////////
// Plot the solutions //
////////////////////////

scf();
plot_fobj_truss(x0);
xtitle('Before optimization','x','y','z');

scf();
plot_fobj_truss(x_opt);
xtitle('After optimization','x','y','z');

warning(warning_old);

