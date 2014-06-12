lines(0);

// Constr pb from example 2
deff('y=f(x)','y=x(1)^2 - 5*x(1) + x(2)^2 - 5*x(2) + 2*x(3)^2 - 21*x(3) + x(4)^2 + 7*x(4) + 50');
deff('y=df(x)','y(1) = 2*x(1) - 5; ...
                y(2) = 2*x(2) - 5; ...
                y(3) = 4*x(3) - 21; ...
                y(4) = 2*x(4) + 7;');

deff('[y,dy] = fobj(x)','y = f(x); dy = df(x);');

deff('[y,dy,ic]=ineqconstraint(x,ct)','Index = 1; ...
                                       y   = 0; ...
                                       dy  = 0; ...
                                       ic  = []; ...
                                       y(1) = x(1)^2 + x(1) + x(2)^2 - x(2) + x(3)^2 + x(3) + x(4)^2 - x(4) - 8 - ct; ...
                                       if y(1)>0 then ...
                                         dy(1,Index) = 2*x(1) + 1; ...
                                         dy(2,Index) = 2*x(2) - 1; ...
                                         dy(3,Index) = 2*x(3) + 1; ...
                                         dy(4,Index) = 2*x(4) - 1; ...
                                         ic(Index) = 1; ...
                                         Index = Index + 1; ...
                                       end ...
                                       y(2) = x(1)^2 - x(1) + 2*x(2)^2 + x(3)^2 + 2*x(4)^2 - x(4) - 10.0 - ct; ...
                                       if y(2)>0 then ...
                                         dy(1,Index) = 2*x(1) - 1; ...
                                         dy(2,Index) = 4*x(2); ...
                                         dy(3,Index) = 2*x(3); ...
                                         dy(4,Index) = 4*x(4) - 1; ...
                                         ic(Index) = 2; ...
                                         Index = Index + 1; ...
                                       end ...
                                       y(3) = 2*x(1)^2 + 2*x(1) + x(2)^2 - x(2) + x(3)^2 - x(4) - 5.0 - ct; ...
                                       if y(3)>0 then ...
                                         dy(1,Index) = 4*x(1) + 2; ...
                                         dy(2,Index) = 2*x(2) - 1; ...
                                         dy(3,Index) = 2*x(3); ...
                                         dy(4,Index) = -1; ...
                                         ic(Index) = 3; ...
                                       end;','c');                                     

x0   = ones(4,1);
ncon = 3;
upper = 99999*ones(4,1);
lower = -99999*ones(4,1);
ItMX = 40;

clear param;
param = init_param();

param = add_param(param,'iprint',4);
param = add_param(param,'nfdg',1);
param = add_param(param,'nacmx',4); // Max number of active constraints at a time

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
