lines(0);

clear f;
clear ineqconstraint;
clear eqconstraint;

f = [];
ineqconstraint = [];
eqconstraint = [];
x0 = [];

////////////////////////////////////////////////////////////////////////
//constr_pb_1
deff('[y,dy]=f(x)','y=4*x(1) - x(2)^2 - 12; ...
                    dy(1,1) = 4; ...
                    dy(2,1) = -2*x(2);');

deff('[y,dy]=ineqconstraint(x)','y(1,1) = - 10*x(1) + x(1)^2 - 10*x(2) + x(2)^2 + 34; ...
                                 dy(1,1) = -10 + 2*x(1); ...
                                 dy(2,1) = -10 + 2*x(2);');
 
deff('[y,dy]=eqconstraint(x)','y(1) = 20 - x(1)^2 - x(2)^2; ...
                               dy(1,1) = -2*x(1); ...
                               dy(2,1) = -2*x(2);');

upper = [15;15];
lower = [-15;-15];
x0    = [1.5;1.5]; // Feasible starting point

/////////////////////////////////////////////////////////////////////

ItMX     = 20; // 40
Log      = %T;

params = init_param();
params = add_param(params,'srand',12345);
//params = add_param(params,'nb_ceq',1);
//params = add_param(params,'nb_cineq',1);

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

ierr = [];

for i=0:40
  printf('Selecting method %d - ', i);
  params = add_param(params,'method',i);

  ierr = execstr('[f_opt, x_opt, status, meth_name] = nlopt(x0, f, ineqconstraint, eqconstraint, lower, upper, ItMX, params);','errcatch');
 
  if (ierr) then
    printf('Method not available\n');
  else
    printf('Method %s available.\n', meth_name);
  end
end
