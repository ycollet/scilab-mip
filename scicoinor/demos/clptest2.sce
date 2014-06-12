// Two examples to show how to solve an MILP and an LP with interior method

printf('-- Integer problem --\n');

c = [-1,-1];
a = [-2, 5; ...
      2,-2];
b = [5,1]';

lb      = [0,0];
ub      = 100*[1,1];
vartype = 'II';
//vartype = 'CC';
constrtype = 'UU';

param = init_param();
param = add_param(param,'maxnumiterations',10000);
param = add_param(param,'maxnumseconds',10000);
param = add_param(param,'primaltolerance',1e-7);
param = add_param(param,'dualtolerance',1e-7);
param = add_param(param,'verbose',1);
param = add_param(param,'solver',1); // 3 interior - other simplex
param = add_param(param,'optim_dir', 1); // optimisation direction: 1 - minimize, -1 - maximize, 0 - ignore

[xmin,lambda,status] = clp([],c,a,b,b,lb,ub,constrtype,vartype,param);

printf('solution found: \n');disp(xmin);
printf('status = %d\n',status); disp(lambda);
pause;

printf('3rd problem\n');
c = [ 0 0 0 -1 -1];
a = [-2 0 0  1  0;...
      0 1 0  0  2;...
      0 0 1  3  2];
b       = [4,12,18]';
lb      = [0,0,0,0,0];
ub      = 100*[1,1,1,1,1];
vartype = 'IIIII';
//vartype = 'CCCCC';
constrtype = 'UUU';

[xmin,lambda,status] = clp([],c,a,b,b,lb,ub,constrtype,vartype,param);

param = init_param();
param = add_param(param,'maxnumiterations',10000);
param = add_param(param,'maxnumseconds',10000);
param = add_param(param,'primaltolerance',1e-7);
param = add_param(param,'dualtolerance',1e-7);
param = add_param(param,'verbose',1);
param = add_param(param,'solver',1); // 3 interior - other simplex
param = add_param(param,'optim_dir', 1); // optimisation direction: 1 - minimize, -1 - maximize, 0 - ignore

printf('solution found: \n');disp(xmin);
printf('status = %d\n',status); disp(lambda);

	       
