lines(0);

////////////////////////////////////////////////////////////////////////

current_dir = pwd();

// double sci_fobj(int size_x, const double * x, double * gradient, double * con, void * param);

f_C = ['#include <math.h>'
       '#include <stdlib.h>'
       'int f_C(int size_x, const double * x, double * gradient, double * con, void * param)'
       '{'
       '  double tmp;'
       '  tmp = 5*pow(x[0],2) - 3*pow(x[1],2);'
       '  if (gradient!=NULL)'
       '    {'
       '      gradient[0] = 10*x[0];'
       '      gradient[1] = -6*x[1];'
       '    }'
       '  return tmp;'
       '}'];
       
// void sci_constr_ineq(unsigned int size_constr, double * result, unsigned int size_x, const double * x, double * gradient, void * param);

g_C = ['#include <math.h>'
       '#include <stdlib.h>'
       'void g_C(unsigned int size_constr, double * result, unsigned int size_x, const double * x, double * gradient, void * param)'
       '{'
       '  result[0] = -x[0];'
       '  result[1] = -x[1];'
       '  if (gradient!=NULL)'
       '    {'
       '      gradient[0] = -1;'
       '      gradient[1] = 0;'
       '      gradient[2] = 0;'
       '      gradient[3] = -1;'
       '    }'
       '}'];
        
cd TMPDIR;
mputl(f_C, TMPDIR+'/nlopt_demo_f_C.c');
mputl(g_C, TMPDIR+'/nlopt_demo_g_C.c');

// compile the C code
printf('Compilation of the f_C function\n');
f_C_ilib_handle  = ilib_for_link('f_C', 'nlopt_demo_f_C.c', [],'c');
printf('Compilation of the g_C function\n');
g_C_ilib_handle  = ilib_for_link('g_C', 'nlopt_demo_g_C.c', [],'c');

// incremental linking
f_C_link_handle  = link(f_C_ilib_handle, 'f_C', 'c');
g_C_link_handle  = link(g_C_ilib_handle, 'g_C', 'c');

cd(current_dir);

upper = [4;4];
lower = [-4;-4];
x0    = [1;1];
  
////////////////////////////////////////////////////////////////////////

params = init_param();
ItMX   = 2000; // 40
Log    = %T;

params = init_param();
params = add_param(params,'srand',12345);
params = add_param(params,'nb_ceq',0);
params = add_param(params,'nb_cineq',2);

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
[f_opt, x_opt, status, meth_name] = nlopt(x0, 'f_C', 'g_C', [], lower, upper, ItMX, params);
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

// unlink the optimization problem

ulink([f_C_link_handle, g_C_link_handle]);
