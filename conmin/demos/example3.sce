lines(0);

// Constr pb from example 3
deff('y=f(x)','y = 10*0.1*(2*sqrt(2)*x(1) + x(2));');
deff('y=df(x)','y(1) = 10*0.1*2*sqrt(2); ...
                y(2) = 10*0.1;');

deff('[y,dy] = fobj(x)','y = f(x); dy = df(x);');

deff('[y,dy,ic]=ineqconstraint(x,ct)','Index = 1; ...
                                       y   = 0; ...
                                       dy  = []; ...
                                       ic  = []; ...
				                        Denom = 2*x(1)*x(2) + sqrt(2)*x(1)*x(1); ...
                                       Sig11 = 20*(sqrt(2)*x(1) + x(2)) / Denom; ...
                                       Sig21 = 20*sqrt(2)*x(1) / Denom; ...
                                       Sig31 = -20*x(2) / Denom; ...
                                       y(1) = - Sig11 / 15 - 1 - ct; ...
                                       y(2) = Sig11 / 20 - 1 - ct; ...
                                       y(3) = - Sig21 / 15 - 1 - ct; ...
                                       y(4) = Sig21 / 20 - 1 - ct; ...
                                       y(5) = - Sig31 / 15 - 1 - ct; ...
                                       y(6) = Sig31 / 20 - 1 - ct;');                                     

x0   = ones(2,1);
ncon = 6;
upper = 1e10  * ones(2,1);
lower = 0.001 * ones(2,1);
ItMX  = 40;

clear param;
param = init_param();

param = add_param(param,'iprint',2);
param = add_param(param,'linobj',1);

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
