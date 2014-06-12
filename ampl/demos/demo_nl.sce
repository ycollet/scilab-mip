// The ch3.nl problem (in Solve_asl + nl_ndex = 1) contains the following problem
// ===============================================================
// CHEBYQUAD FUNCTION (R. Fletcher, "Function Minimization Without
//		      Evaluating Derivatives -- A Review",
//		      Computer J. 8 (1965), pp. 163-168.)
// ===============================================================
//
// param n > 0;
//
// var x {j in 1..n} := j/(n+1);
//
// var T {i in 0..n, j in 1..n} =
//          if (i = 0) then 1 else
//          if (i = 1) then 2*x[j] - 1 else
//          2 * (2*x[j]-1) * T[i-1,j] - T[i-2,j];
//
// =======================
// Equation form (for nl2)
// =======================
//
// eqn {i in 1..n}:
//    (1/n) * sum {j in 1..n} T[i,j] = if (i mod 2 = 0) then 1/(1-i^2) else 0;
//
// ========================
// Objective form (for mng)
// ========================
//
// minimize ssq: 0.5*sum{i in 1..n}
//      ((1/n) * sum {j in 1..n} T[i,j] - if (i mod 2 = 0) then 1/(1-i^2))^2;
// data;
//  param n := 3;

lines(0);

stacksize('max');

path = get_absolute_file_path('demo_nl.sce');

exec(path + 'nl_data.sce');

Solve_macminlp = %F; // 92  files
Solve_coinor   = %F; // 127 files
Solve_asl      = %T; // 10  files
Solve_modnl    = %F; // 898 files

///////////////////////////
// Set global paremeters //
///////////////////////////

nl_index = 1; // 325CoinOR + Index 13  = hs100
              // ModNL  + Index 440 = rosenbr
                
if Solve_macminlp then nl_filename = MacMINLP(nl_index); end
if Solve_coinor   then nl_filename = CoinOR(nl_index);   end
if Solve_asl      then nl_filename = ASL(nl_index);      end
if Solve_modnl    then nl_filename = ModNL(nl_index);    end

///////////////////////////////
// Load and test the problem //
///////////////////////////////

printf('Optimization of the %s problem.\n\n',basename(nl_filename));

[asl, x0, bl, bu, v, cl, cu] = ampl_init(nl_filename);

[f, c] = ampl_evalf(asl, x0);
printf('\nObjective function informations:\n');
printf('objective function value: %d\n', f);

[g, jac] = ampl_evalg(asl, x0);
printf('constraint values:'); disp(g');
printf('\nGradient informations:\n');
printf('gradient of the objective function:'); disp(g');
printf('jacobian of the constraints:'); disp(jac');

ampl_write_sol('this is a message',x0,v);

W = ampl_evalw(asl,v);

printf('\nHessian informations:\n');
printf('- the Hessian is %d x %d\n', size(W,1), size(W,2));
printf('- there is %d non zeros values in the Hessian\n',length(find(W~=0)));

/////////////////////////////////////////////
// Validation of the sparse representation //
/////////////////////////////////////////////

[spg, spjac] = ampl_eval_sp_g(asl, x0);
printf('The difference between the dense Jacobian and the sparse one is: %f\n', norm(jac-full(spjac')));

spW = ampl_eval_sp_w(asl,v);
printf('The difference between the dense Hessian and the sparse one is: %f\n', norm(W-full(spW)));

//////////////////////////////////////////////////
// Get some informations related to the problem //
//////////////////////////////////////////////////

info = ampl_get_size(asl);

printf('\nProblem informations:\n');
printf('number of linear binary variables:     %d\n',info('nbv'));
printf('number of linear integer variables:    %d\n',info('niv'));
printf('total number of nonlinear constraints: %d\n',info('nlc'));
printf('number of equality constraints or -1 if unknown (ampl prior to 19970627) : %d\n',info('n_eqn'));
printf('total complementarity conditions:        %d\n',info('n_cc'));
printf('nonlinear complementarity conditions:    %d\n',info('nlcc'));
printf('number of nonlinear network constraints: %d\n',info('nlnc'));
printf('number of nonlinear objectives:          %d\n',info('nlo'));
printf('number of nonlinear variables in both constraints and objectives: %d\n',info('nlvb'));
printf('number of nonlinear variables in constraints:                     %d\n',info('nlvc'));
printf('number of nonlinear variables in objectives nlvc and nlvo include nlvb: %d\n',info('nlvo')); 
printf('integer nonlinear variables in both constraints and objectives :        %d\n',info('nlvbi'));
printf('integer nonlinear vars just in constraints :    %d\n',info('nlvci'));
printf('integer nonlinear vars just in objectives:      %d\n',info('nlvoi'));
printf('number of (linear) network variables (arcs):    %d\n',info('nwv'));
printf('number of nonzeros in constraints Jacobian:     %d\n',info('nzc'));
printf('number of nonzeros in all objective gradients : %d\n',info('nzo'));
printf('total number of variables:     %d\n',info('n_var'));
printf('total number of constraints:   %d\n',info('n_con'));
printf('total number of objectives:    %d\n',info('n_obj'));
printf('number of logical constraints: %d\n',info('n_lcon'));

///////////////////////////////////
// Get Complementarity condition //
///////////////////////////////////

cvar = ampl_get_compl(asl);

// if cvar == -1 then, there are no complementarity constraints
printf('Complementarity conditions. List of indexes:'); disp(cvar);

///////////////////////////
// Start an optimization //
///////////////////////////

Algorithm = 'qn';

// Some function seems to be unbounded. So, we replace the %inf bounds by finite values
Bound = 100;
bu(find(bu>=%inf))  =  Bound;
bu(find(bu<=-%inf)) = -Bound;
bl(find(bl>=%inf))  =  Bound;
bl(find(bl<=-%inf)) = -Bound;

xstart = (bu - bl) .* (2*rand(bu)-1) + bl;

xstart = x0;

deff('[f, g, ind] = fobj(x,ind)','[tmp_fobj, tmp_gobj] = ampl_evalf(asl,x); ...
                                  f = tmp_fobj + 100 * sum(tmp_gobj); ...
                                  [tmp_dfobj, tmp_dgobj] = ampl_evalg(asl,x); ...
                                  g = tmp_dfobj + 100 * sum(tmp_dgobj'',''c''); ...
                                  disp(g'');');
                                  
[fopt, xopt] = optim(fobj,'b',bl,bu,xstart,algo=Algorithm,imp=2);

printf('\nOptimization information:\n');

[tmp_fobj, tmp_gobj, ind] = fobj(xstart,2);

printf('initial point: fobj(xstart) = %f\n',tmp_fobj); 
printf('fopt = %f\n', fopt);

[g_opt, jac_opt] = ampl_evalg(asl, xopt);

if (~or(isnan(g_opt)) & ~or(isnan(jac_opt))) then
  printf('norm of the gradient: %f\n', norm(g_opt));
  printf('norm of the jacobian: %f\n', norm(jac_opt));
end
printf('the gradient: '); disp(g_opt');

////////////////////////////////////////////////////
// Check the derivative of the objective function //
////////////////////////////////////////////////////

deff('yvar = fobj_check_deriv(xvar)','[yvar,tmp] = ampl_evalf(asl,xvar);');

printf('Check performed at starting point:\n');
printf('derivative of the objective function given by ampl:');
[g, jac] = ampl_evalg(asl,x0);
disp(g');

printf('derivative of the objective function given by ''derivative'':');
g = derivative(fobj_check_deriv,x0)';
disp(g');

printf('Check performed at solution:\n');
printf('derivative of the objective function given by ampl:');
[g, jac] = ampl_evalg(asl,xopt);
disp(g');

printf('derivative of the objective function given by ''derivative'':');
g = derivative(fobj_check_deriv,xopt)';
disp(g');

/////////////////////////////////////////////
// Check the derivative of the constraints //
/////////////////////////////////////////////

deff('yvar = constr_check_deriv(xvar)','[tmp,yvar] = ampl_evalf(asl,xvar);');

printf('Check performed at starting point:\n');
printf('derivative of the constraints given by ampl:');
[g, jac] = ampl_evalg(asl,x0);
disp(jac');

printf('derivative of the constraints given by ''derivative'':');
jac = derivative(constr_check_deriv,x0);
disp(jac');

printf('Check performed at solution:\n');
printf('derivative of the constraints given by ampl:');
[g, jac] = ampl_evalg(asl,xopt);
disp(jac');

printf('derivative of the constraints given by ''derivative'':');
jac = derivative(constr_check_deriv,xopt);
disp(jac');

/////////////////////////////////////////////////////
// Check the sparsity structure of the constraints //
/////////////////////////////////////////////////////

deff('yvar = constr_check_deriv(xvar)','[tmp,yvar] = ampl_evalf(asl,xvar);');

printf('Check performed at starting point:\n');
printf('derivative of the constraints given by the sparsity structure functions ampl:');
[irow, jcol] = ampl_eval_spst_g_rc(asl,x0);
val          = ampl_eval_spst_g_val(asl,x0);
spjac = sparse([irow, jcol], val, [length(cl),length(x0)]);
disp(full(spjac)');

printf('derivative of the constraints given by ''derivative'':');
jac = derivative(constr_check_deriv,x0);
disp(jac');
///////////////////////////////////////
// Getting the type of the variables //
///////////////////////////////////////

printf('variables type: %s\n', ampl_get_type(asl));

//////////////////////////
// Writing the solution //
//////////////////////////

printf('\Writing the solution:\n');
printf('The solution of file %s will be written in %s.\n', nl_filename, strsubst(nl_filename,'.nl','.sol'));

ampl_write_sol(asl,'solution of the AMPL ' + nl_filename + ' problem', xopt, v);

///////////////////////////////////////////////
// Get the problem as a direct acyclic graph //
///////////////////////////////////////////////

// NOT YET WORKING: dag = ampl_get_dag(asl);

//////////////////////////////////////////////////
// Free the pointer handled in the asl variable //
//////////////////////////////////////////////////

ampl_free(asl);
asl = [];

