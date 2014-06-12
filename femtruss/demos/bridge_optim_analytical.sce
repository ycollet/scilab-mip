lines(0);
warning_old = warning('query');

Order = 2; // 1 2 or 4 for derivative
Debug = %T;

femtruss_path = get_absolute_file_path('bridge_optim_analytical.sce');
exec(femtruss_path + 'build_long_bridge_2d.sce');

function y = fobj_truss(x)
[t,p,e,A,E,rho,F] = bridge_optim(x);
[U,P,R]= femtruss(bridge_optim, %F, x);

y = sum(U(IndexNodeVarInf).^2);
endfunction

function dy = dfobj_truss(x)
// Here, we will compute analytical derivatives for the objective function.
// We use dfemtruss to compute analytical derivatives for displacement
// We then compute by hand the analytical derivative for the objective function using the result of dfemtruss
[t,p,e,A,E,rho,F] = bridge_optim(x);
[U,P,R]= femtruss(bridge_optim, %F, x);

// The deck - 3 4 5
Pos_deck = localise2d(IndexNodeVarInf);
// The degree of freedom for the optimization - 6 7 8
Pos_free = localise2d(IndexNodeVarSup);

dU = dfemtruss_ana(bridge_optim,U,[],%F, x);

dy = zeros(length(Pos_free),1);

// dU contains only partial derivatives for nodes 2 3 and 4. So, we get partial derivatives with respect to y by accessing
// value 2 4 and 6. Values 1 3 and 5 are partial derivatives wrt x.
// dU(I4,[I1 I2 I3]) = dU(I1) / dI4 dU(I2) / dI4 dU(I3) / dI4
for i=1:length(Pos_free)
  dy(i) = 2*sum(dU(Pos_free(i), Pos_deck(2:2:$)).*U(Pos_deck(2:2:$)));
end  
endfunction

function [y, dy, ind] = optim_fobj_truss(x, ind)
y  = fobj_truss(x);
dy = dfobj_truss(x)';
if Debug then
  printf('y = %f norm(dy) = %f\n',y, norm(dy));
end
endfunction

////////////////////////////////////////////
// Parameters for the optimization method //
////////////////////////////////////////////

MaxEvalFunc = 400;
Algorithm   = 'gc'; // 'qn', 'gc', 'nd' -> Ne marche qu'avec 'qn' (quasi-newton). Pour les autres, on obtient rapidement une structure mal
                    // conditionnee
Log         = %F;
TOL         = 1.0e-12; // accuracy for convergence test (minimum)

//////////////////////////////
// Main part of the problem //
//////////////////////////////

printf('optimization starting, be patient ... \n\n');

[f_opt, x_opt] = optim(optim_fobj_truss, x0, algo=Algorithm,'ar',MaxEvalFunc,MaxEvalFunc,TOL,TOL);

printf('initial solution:'); disp(x0');
printf('initial objective function value = %f\n',fobj_truss(x0));

printf('Final solution:'); disp(x_opt');
printf('Final objective function value = %f\n',fobj_truss(x_opt));
  
scf();
plot_fobj_truss(x0);
xtitle('Before optimization','x','y');

scf();
plot_fobj_truss(x_opt);
xtitle('After optimization','x','y');

warning(warning_old);

