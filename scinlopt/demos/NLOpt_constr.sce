lines(0);

path = get_absolute_file_path('NLOpt_constr.sce');

exec(path + 'cont_constr_funcs.sci');

func_names =  get_constr_pb_name();

item = x_choose(func_names,'Constrained problems');

if item~=0 then
  func_name = func_names(item);
else
  abort();
end

clear f;
clear ineqconstraint;
clear eqconstraint;

f              = [];
ineqconstraint = [];
eqconstraint   = [];

execstr('res   = ' + func_name + '_has_df_obj();');
if res then
  deff('[y,dy]=f(x)','y  = ' + func_name + '_obj(x); ...
                      dy = df_' + func_name + '_obj,x);');
else
  deff('[y,dy]=f(x)','y  = ' + func_name + '_obj(x); ...
                      dy = derivative(' + func_name + '_obj,x);');
end

execstr('res   = ' + func_name + '_has_df_ineq();');
if res then
  deff('[y,dy]=ineqconstraint(x)','y  = ' + func_name + '_ineq(x); ...
                                   dy = df_' + func_name + '_ineq,x);');
else
  deff('[y,dy]=ineqconstraint(x)','y  = ' + func_name + '_ineq(x); ...
                                   dy = derivative(' + func_name + '_ineq,x);');
end

execstr('res   = ' + func_name + '_has_df_eq();');
if res then
  deff('[y,dy]=eqconstraint(x)','y  = ' + func_name + '_eq(x); ...
                                 dy = df_' + func_name + '_eq,x);');
else
  deff('[y,dy]=eqconstraint(x)','y  = ' + func_name + '_eq(x); ...
                                 dy = derivative(' + func_name + '_eq,x);');
end

execstr('upper = max_bd_' + func_name + '();');
execstr('lower = min_bd_' + func_name + '();');
execstr('x0    = ' + func_name + '_init();');

execstr('res   = ' + func_name + '_has_ineq();');
if res then
  [y_tmp,dy_tmp] = ineqconstraint(x0);
  nb_cineq       = length(y_tmp);
else
  clear ineqconstraint;
  ineqconstraint = [];
  nb_cineq       = 0;
end

execstr('res   = ' + func_name + '_has_eq();');
if res then
  [y_tmp,dy_tmp] = eqconstraint(x0);
  nb_ceq         = length(y_tmp);
else
  clear eqconstraint;
  eqconstraint = [];
  nb_ceq       = 0;
end

/////////////////////////////////////////////////////////////////////

ItMX = 2000; // 40
Log  = %T;
Plot = %T;

// Plot the graph

if Plot & length(x0)==2 then
  scf(10001);
  drawlater;
  Z_graph = [];
  x_graph = lower(1):(upper(1) - lower(1))/20:upper(1);
  y_graph = lower(2):(upper(2) - lower(2))/20:upper(2);
  for i=1:size(x_graph,2)
    for j=1:size(y_graph,2)
      Z_graph(i,j) = f([x_graph(i) y_graph(j)]');
    end
  end
  xset('fpf',' ');
  contour(x_graph,y_graph,Z_graph, 10);

  // We plot the inequality constraints
  if (ineqconstraint~=[]) then
    Aux = ineqconstraint(x0);
    for k=1:length(Aux)
      for i=1:size(x_graph,2)
        for j=1:size(y_graph,2)
          tmp = ineqconstraint([x_graph(i) y_graph(j)]');
          Z_ineq(i,j) = tmp(k);
        end
      end
      xset('fpf',' ');
      contour(x_graph,y_graph,Z_ineq,[0 0]);
    end
  end

  // We plot the equality constraints
  if (eqconstraint~=[]) then
    Aux = eqconstraint(x0);
    for k=1:length(Aux)
      for i=1:size(x_graph,2)
        for j=1:size(y_graph,2)
          tmp = eqconstraint([x_graph(i) y_graph(j)]');
          Z_eq(i,j) = tmp(k);
        end
      end
      xset('fpf',' ');
      contour(x_graph,y_graph,Z_eq,[0 0]);
    end
  end
  drawnow;
end

params = init_param();
params = add_param(params,'srand',12345);
params = add_param(params,'nb_ceq',nb_ceq);
params = add_param(params,'nb_cineq',nb_cineq);

// Naming conventions for the optimization methods:
//
// NLOPT_{G/L}{D/N}_* 
//    = global/local derivative/no-derivative optimization, respectively 
//
// * _RAND algorithms involve some randomization.
// * _NOSCAL algorithms are *not* scaled to a unit hypercube (i.e. they are sensitive to the units of x)
//
// 0  - NLOPT_GN_DIRECT = 0,
// 1  - NLOPT_GN_DIRECT_L,
// 2  - NLOPT_GN_DIRECT_L_RAND,
// 3  - NLOPT_GN_DIRECT_NOSCAL,
// 4  - NLOPT_GN_DIRECT_L_NOSCAL,
// 5  - NLOPT_GN_DIRECT_L_RAND_NOSCAL,
// 6  - NLOPT_GN_ORIG_DIRECT,
// 7  - NLOPT_GN_ORIG_DIRECT_L,
// 8  - NLOPT_GD_STOGO,
// 9  - NLOPT_GD_STOGO_RAND,
// 10 - NLOPT_LD_LBFGS_NOCEDAL,
// 11 - NLOPT_LD_LBFGS,
// 12 - NLOPT_LN_PRAXIS,
// 13 - NLOPT_LD_VAR1,
// 14 - NLOPT_LD_VAR2,
// 15 - NLOPT_LD_TNEWTON,
// 16 - NLOPT_LD_TNEWTON_RESTART,
// 17 - NLOPT_LD_TNEWTON_PRECOND,
// 18 - NLOPT_LD_TNEWTON_PRECOND_RESTART,
// 19 - NLOPT_GN_CRS2_LM,
// 20 - NLOPT_GN_MLSL,
// 21 - NLOPT_GD_MLSL,
// 22 - NLOPT_GN_MLSL_LDS,
// 23 - NLOPT_GD_MLSL_LDS,
// 24 - NLOPT_LD_MMA,
// 25 - NLOPT_LN_COBYLA,
// 26 - NLOPT_LN_NEWUOA,
// 27 - NLOPT_LN_NEWUOA_BOUND,
// 28 - NLOPT_LN_NELDERMEAD,
// 29 - NLOPT_LN_SBPLX,
// 30 - NLOPT_LN_AUGLAG,
// 31 - NLOPT_LD_AUGLAG,
// 32 - NLOPT_LN_AUGLAG_EQ,
// 33 - NLOPT_LD_AUGLAG_EQ,
// 34 - NLOPT_LN_BOBYQA,
// 35 - NLOPT_GN_ISRES,
// 
// new variants that require local_optimizer to be set, not with older constants for backwards compatibility
// 
// 36 - NLOPT_AUGLAG,
// 37 - NLOPT_AUGLAG_EQ,
// 38 - NLOPT_G_MLSL,
// 39 - NLOPT_G_MLSL_LDS,
// 40 - NLOPT_LD_SLSQP,

params = add_param(params,'method',31);
//params = add_param(params,'max',1);
params = add_param(params,'cineq_tol',[1e-4]);
params = add_param(params,'ceq_tol',[1e-4]);
//params = add_param(params,'stopval',1);
params = add_param(params,'ftol_rel',1e-12);
params = add_param(params,'ftol_abs',1e-12);
params = add_param(params,'xtol_rel',1e-6);
//params = add_param(params,'xtol_abs1',1);
//params = add_param(params,'xtol_abs',1*ones(x0));
params = add_param(params,'maxtime',1000);
//params = add_param(params,'force_stop',0);
//params = add_param(params,'population',20);
//params = add_param(params,'default_initial_step',1e-3*ones(x0));
//params = add_param(params,'initial_step',1e-3*ones(x0));
//params = add_param(params,'initial_step1',1e-3);

printf('Starting optimization, be patient ...\n');
t_start = timer();
[f_opt, x_opt, status, meth_name] = nlopt(x0, f, ineqconstraint, eqconstraint, lower, upper, ItMX, params);
t_end = timer();

printf('Time required for optimization = %f\n',t_end - t_start);
printf('name of the optimization method: %s\n\n', meth_name);

printf('Status possible values:\n');
printf('-----------------------\n\n');
printf('In case of success:\n\n');
printf('1: NLOPT_SUCCESS Generic success return value.\n');
printf('2: NLOPT_STOPVAL_REACHED Optimization stopped because stopval (above) was reached.\n');
printf('3: NLOPT_FTOL_REACHED Optimization stopped because ftol_rel or ftol_abs (above) was reached.\n');
printf('4: NLOPT_XTOL_REACHED Optimization stopped because xtol_rel or xtol_abs (above) was reached.\n');
printf('5: NLOPT_MAXEVAL_REACHED Optimization stopped because maxeval (above) was reached.\n');
printf('6: NLOPT_MAXTIME_REACHED Optimization stopped because maxtime (above) was reached.\n\n');

printf('In case of failure\n\n');
printf('-1: NLOPT_FAILURE Generic failure code.\n');
printf('-2: NLOPT_INVALID_ARGS Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etcetera).\n');
printf('-3: NLOPT_OUT_OF_MEMORY Ran out of memory.\n');
printf('-4: NLOPT_ROUNDOFF_LIMITED Halted because roundoff errors limited progress.\n');
printf('-5: NLOPT_FORCED_STOP\n');

printf('status = %d\n', status);

printf('value of the objective function: %f\n', f_opt);
if ineqconstraint~=[] then
  printf('value of the inequality constraint:');
  [y,dy] = ineqconstraint(x_opt); disp(y);
end
if eqconstraint~=[] then
  printf('value of the equality constraint:');
  [y,dy] = eqconstraint(x_opt); disp(y);
end
printf('value of the solution:'); disp(x_opt');

if (Plot) & (length(x0)==2) then
  printf('plot the points\n');
  plot(x0(1),x0(2),'go');
  plot(x_opt(1),x_opt(2),'ro');
  legends(['Starting pt';'Final pt'],[3 5],4);
end
