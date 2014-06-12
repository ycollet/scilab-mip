// A LP example which shows all the potentials of CLP Scilab

printf('LP problem\n');

c = -[10, 6, 4];
a = [1,  1, 1;...
     10, 4, 5;...
     2,  2, 6];
b = [100,600,300]';
lb    = [0,0,0];
ub    = [100,100,100];
vartype = 'CCC';
constrtype = 'LLL';

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

