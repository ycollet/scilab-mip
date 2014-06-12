// A LP example which defines A as a sparse matrix

printf('Solve problem with sparse matrix\n');

c = [0 0 0 -1 -1];
a = sparse([-2 0 0 1 0;...
             0 1 0 0 2;...
             0 0 1 3 2]);
//a = [-2 0 0 1 0;...
//      0 1 0 0 2;...
//      0 0 1 3 2];
rhs       = [4 12 18]';
lhs       = [-4 -12 -18]';
lb      = [0,0,0,0,0];
ub      = 100*[1,1,1,1,1];
vartype = 'IIIII';
//vartype = 'CCCCC';
constrtype = 'UUU';

param = init_param();
param = add_param(param,'maxnumiterations',100000);
param = add_param(param,'maxnumseconds',100000);
param = add_param(param,'primaltolerance',1e-12);
param = add_param(param,'dualtolerance',1e-12);
param = add_param(param,'verbose',32);
param = add_param(param,'solver',5); // 6 interior - other simplex 
param = add_param(param,'optim_dir', 1); // optimisation direction: 1 - minimize, -1 - maximize, 0 - ignore
param = add_param(param,'fact_freq', 10); 
param = add_param(param,'perturb', 200);
param = add_param(param,'presolve', 0); 
param = add_param(param,'red_grad', 2); 
param = add_param(param,'writemps','test.mps');

// Probleme avec 1
[xmin,lambda,status] = clp([],c,a,lhs,rhs,lb,ub,constrtype,vartype,param);

printf('solution found: \n');disp(xmin);
printf('status = %d\n',status); disp(lambda);
