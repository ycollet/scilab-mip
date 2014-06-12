// A LP example which defines A as a sparse matrix

lines(0);

printf('Solve problem with sparse matrix\n');

c = [0 0 0 -1 -1]';
a = sparse([-2 0 0 1 0;...
             0 1 0 0 2;...
             0 0 1 3 2]);
b       = [4 12 18]';
ctype   = "GGG";
lb      = [0,0,0,0,0]';
ub      = 100*[1,1,1,1,1]';
vartype = "CCCCC";

param = init_param();
param = add_param(param,'msglev',3);

[xmin,fmin,status,extra] = glpk(c,a,b,b,lb,ub,ctype,vartype,param);
printf('solution found: fmin = %f\n', fmin);disp(xmin);
printf('status = %d\n',status); disp(extra);
