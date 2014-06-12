lines(0);

// constr_pb_2
deff('y=f(x)','y=x(1)^2 + x(2)^2 - 16*x(1) - 10*x(2)');
deff('y=df(x)','y = derivative(f,x)''');
deff('[y,dy] = fobj(x)','y = f(x); dy = df(x);');

deff('[y,dy,ic]=ineqconstraint(x,ct)','Index = 1; ...
                                       y   = 0; ...
                                       dy  = 0; ...
                                       ic  = []; ...
                                       tmp = - 11 + x(1)^2 - 6*x(1) + 4*x(2) - ct; ...
                                       if tmp>0 then ...
                                         y(Index)  = tmp; ...
                                         dy(1,Index) = 2*x(1) - 6; ...
                                         dy(2,Index) = 4; ...
                                         ic(Index) = 1; ...
                                         Index = Index + 1; ...
                                       end ...
                                       tmp = - x(1)*x(2) + 3*x(2) + exp(x(1) - 3) - 1 - ct; ...
                                       if tmp>0 then ...
                                         y(Index)  = tmp; ... 
                                         dy(1,Index) = -x(2) + exp(x(1) - 3); ...
                                         dy(2,Index) = -x(1) + 3; ...
                                         ic(Index) = 2; ...
                                         Index = Index + 1; ...
                                       end ...
                                       tmp = - x(1) - ct; ...
                                       if tmp>0 then ...
                                         y(Index)  = tmp; ...
                                         dy(1,Index) = -1; ...
                                         dy(2,Index) = 0; ...
                                         ic(Index) = 3; ...
                                         Index = Index + 1; ...
                                       end ...
                                       tmp = - x(2) - ct; ...
                                       if tmp>0 then ...
                                         y(Index)  = tmp; ...
                                         dy(1,Index) = 0; ...
                                         dy(2,Index) = -1; ...
                                         ic(Index) = 4; ...
                                         Index = Index + 1; ...
                                       end;','c');                                     

x0    = [4; 3];
ncon  = 4;
upper = [15;9];
lower = [0; 0];
ItMX  = 40;


clear param;
param = init_param();

param = add_param(param,'isc',   [0; 0; 0; 0]);

// ISC(N2) Not used if NCON = 0.  Linear constraint identification vector. If constraint G(J) is known to be a linear function of the decision variables, X(I), ISC(I) should be initialized to
//         ISC(I) = 1.  If constraint G(J) is nonlinear ISC(I) is initialized to ISC(I) = 0.  Identification of linear constraints may improve
//         efficiency of the optimization process and is therefore desirable, but is not essential.  If G(J) is not specifically known to be linear, set ISC(I) = 0.
// isc = [0 0 0 0];

// NSCAL   Scaling control parameter.  The decision variables will be scaled linearly.
//         NSCAL.LT.0:  Scale variables X(I) by dividing by SCAL(I), where vector SCAL is defined by the user.
//         NSCAL.EQ.0:  Do not scale the variables.
//         NSCAL.GT.0:  Scale the variables every NSCAL iterations. Variables are normalized so that scaled X(I) = X(I)/ABS(X(I)).  When using this option, it
//                      is desirable that NSCAL = ICNDIR if ICNDIR is input as nonzero, and NSCAL = NDV + 1 in ICNDIR is input as zero.
param = add_param(param,'nscal',0);

// SCAL(N5) Not used if NSCAL = 0.  Vector of scaling parameters.  If NSCAL>0 vector SCAL need not be initialized since SCAL will be defined in CONMIN and its associated routines.
//          If NSCAL<0, vector SCAL is initialized in the main program, and the scaled variables X(I) = X(I)/SCAL(I).  Efficiency of the optimization
//          process can sometimes be improved if the variables are either normalized or are scaled in such a way that the partial derivative of the objective function, OBJ, 
//          with respect to variable X(I) is of the same order of magnitude for all X(I).  SCAL(I) must be greater than zero because a negative value of SCAL(I)
//          will result in a change of sign of X(I) and possibly yield erroneous optimization results.  The decision of if, and how, the variables should be scaled is highly 
//          problem dependent, and some experimentation is desirable for any given class of problems.
param = add_param(param,'scal',ones(size(x0,1),size(x0,2)));

// INFOG = 0:   same as when INFOG was not used.
// INFOG = 1:   only those constraints identified as active or violated in array IC(I), I = 1, NAC need be evaluated.  This is only meaningful if finite difference gradients are
//              calculated, and allows the user to avoid calculating non-essential information.  If it is convenient to evaluate all constraints each time, variable INFOG may be ignored.
param = add_param(param,'infog',0);

// INFO = 1   Calculate objective function value, OBJ, for current variables X.
// INFO = 2   Calculate objective function value, OBJ, and constraint values, G(J), J = 1, NCON for current variables, X.
// INFO = 3   Calculate analytic gradient of objective function corresponding to current variables, X.  The objective function and constraint values already correspond to the
//            current values of X and need not be recalculated. However, other information obtained in SUB1 when calculating OBJ and G(J) may not correspond to X and must
//            be calculated again here if it is used in gradient computations.  If finite difference control parameter, NFDG, is set to NFDG = 1 in the main program this value
//            of INFO will never be considered.
// INFO = 4   For current variables, X, determine which constraints are active and which are violated (G(J).GE.CT) and how many such constraints there are (NAC = Number of active
//            and violated constraints).  Calculate the analytic gradients of the objective function and all active or violated constraints.  Values of the objective function,
//            OBJ, and constraints, G(J), already correspond to the current variables, X, and need not be recalculated. As in the case of INFO = 3, all other information used
//            in gradient computations must be calculated for the current variables, X.  If finite difference control parameter NFDG, defined in the main program, is not zero,
//            this value of INFO will never be considered.
param = add_param(param,'info',0);

// NFDG = 0:  all gradient information is calculated by finite difference within CONMIN.
// NFDG = 1:  all gradient information is supplied by the user.
// NFDG = 2:  the gradient of OBJ is supplied by the user and the gradients of constraints are calculated by finite difference within CONMIN.
param = add_param(param,'nfdg',1);

// IPRINT Print control.  All printing is done on unit number 6.
//         0:  Print nothing.
//         1:  Print initial and final function information.
//         2:  1st debug level.  Print all of above plus control parameters. Print function value and X-vector at each iteration.
//         3:  2nd. debug level.  Print all of above plus all constraint values, numbers of active or violated constraints, direction vectors, move parameters and miscellaneous information.
//             The constraint parameter, BETA, printed under this option approaches zero as the optimum objective is achieved.
//         4:  Complete debug.  Print all of above plus gradients of objective function, active or violated constraint functions and miscellaneous information.Iprint = [];
param = add_param(param,'iprint',4);

// ICNDIR   Default value = NDV + 1.  Conjugate direction restart parameter. If the function is currently unconstrained, (all G(J)<CT or NCON = NSIDE = 0), 
//          Fletcher-Reeves conjugate direction method will be restarted with a steepest descent direction every ICNDIR iterations.  If ICNDIR = 1 only steepest descent will be used.
param = add_param(param,'icndir',length(x0)+1);

// FDCH    Default value = 0.01.  Not used if NFDG = 0.  Relative change in decision variable X(I) in calculating finite difference gradients.  For example, FDCH = 0.01 corresponds to a finite
//         difference step of one percent of the value of the decision variable.
param = add_param(param,'fdch',0.01);

// FDCHM   Default value = 0.01.  Not used if NFDG = 0.  Minimum absolute step in finite difference gradient calculations.  FDCHM applies to the unscaled variable values.
param = add_param(param,'fdchm',0.01);

// CT      Default value = -0.1.  Not used if NCON = NSIDE = 0. Constraint thickness parameter.  If CT<=G(J)<=ABS(CT), G(J) is defined as active.  If G(J)>ABS(CT), G(J) is said to
//         be violated.  If G(J)<CT, G(J) is not active.  CT is sequentially reduced in magnitude during the optimization process.  If ABS(CT) is very small, one or more constraints
//         may be active on one iteration and inactive on the next, only to become active again on a subsequent iteration. This is often referred to as "zigzagging" between constraints.
//         A wide initial value of the constraint thickness is desirable for highly nonlinear problems so that when a constraint becomes active it tends to remain active, thus reducing the
//         zigzagging problem.  The default value is usually adequate.
param = add_param(param,'ct',-0.1);

// CTMIN   Default value = 0.004.  Not used if NCON = NSIDE = 0.  Minimum absolute value of CT considered in the optimization process. CTMIN may be considered as "numerical zero" 
//         since it may not be meaningful to compare numbers smaller than CTMIN.  The value of CTMIN is chosen to indicate that satisfaction of a constraint
//         within this tolerance is acceptable.  The default value is usually adequate.
param = add_param(param,'ctmin',0.004);

// CTL     Default value = -0.01.  Not used if NCON = NSIDE = 0. Constraint thickness parameter for linear and side constraints. CTL is smaller in magnitude than CT because the zigzagging
//         problem is avoided with linear and side constraints.  The default value is usually adequate.
param = add_param(param,'ctl',-0.01);

// CTLMIN   Default value = 0.001.  Not used if NCON = NSIDE = 0.  Minimum absolute value of CTL considered in the optimization process. The default value is usually adequate.
param = add_param(param,'ctlmin',0.001);

// THETA   Default value = 1.0.  Not used if NCON = NSIDE = 0.  Mean value of the push-off factor in the method of feasible directions. A larger value of THETA is desirable 
//         if the constraints, G(J), are known to be highly nonlinear, and a smaller value may be used if all G(J) are known to be nearly linear.  The actual
//         value of the push-off factor used in the program is a quadratic function of each G(J), varying from 0.0 for G(J) = CT to 4.0*THETA for G(J) = ABS(CT).  
//         A value of THETA = 0.0 is used in the program for constraints which are identified by the user to be strictly linear. THETA is called a "push-off" factor because
//         it pushes the design away from the active constraints into the feasible region.  The default value is usually adequate.
param = add_param(param,'theta',1.0);

// DELFUN  Default value = 0.001.  Minimum relative change in the objective function to indicate convergence.  If in ITRM consecutive iterations, ABS(1.0-OBJ(J-1)/OBJ(J))<DELFUN 
//         and the current design is feasible (all G(J)<=ABS(CT)), the minimization process is terminated.  If the current design is infeasible
//         (some G(J)>ABS(CT)), five iterations are required to terminate and this situation indicates that a feasible design may not exist.
param = add_param(param,'delfun',0.001);

// DABFUN   Default value = 0.001 times the initial function value.  Same as DELFUN except comparison is on absolute change in the objective function, ABS(OBJ(J)-OBJ(J-1)), 
//          instead of relative change.
param = add_param(param,'dabfun',0.001);

// LINOBJ  Not used if NCON = NSIDE = 0.  Linear objective function identifier.  If the objective, OBJ, is specifically known to be a strictly linear function of the decision variables, X(I),
//         set LINOBJ = 1.  If OBJ is a general nonlinear function, set LINOBJ = 0.
param = add_param(param,'linobj',0);

// ITRM    Default value = 3.  Number of consecutive iterations to indicate convergence by relative or absolute changes, DELFUN or DABFUN.
param = add_param(param,'itrm',3);

// ALPHAX (default = 0.1) is the maximum fractional change in any component of X as an initial estimate for ALPHA in the one-dimensional search. 
//                        That is, the initial ALPHA will be such that no component of X is changed by more than this amount. 
//                        This only applies to those X(i) of magnitude greater than 0.1. If an optimization run shows numerous ALPHA = 0 results for the one-dimensional search, 
//                        it may help to try ALPHAX less than the default. ALPHAX is changed by CONMIN depending on the progress of the optimization.
param = add_param(param,'alphax',0.1);

// ABOBJ1 (default = 0.1) is the fractional change attempted as a first step in the one-dimensional search and is based on a linear approximation.
//                        ABOBJ1 is updated during the optimization, depending on progress. The initial step in the one-dimensional search is taken as the amount necessary 
//                        to change OBJ by ABOBJ1*ABS(OBJ) or to change some X(i) by ALPHAX*ABS( X(i) ), whichever is less.
param = add_param(param,'abobj1',0.1);

x_opt = conmin_optim(x0,fobj,ineqconstraint,ncon,upper,lower,ItMX,param);

printf('Initial values\n');
printf('Value of the starting solution: ');
disp(x0');
printf('Value of the objective function: %f\n',f(x0));
printf('Value of the constraints:');
disp(ineqconstraint(x0,0)');

printf('Final values\n');
printf('Value of the optimal solution: ');
disp(x_opt');
printf('Value of the objective function: %f\n',f(x_opt));
printf('Value of the constraints:');
disp(ineqconstraint(x_opt,0)');
