// A LP example which shows all the potentials of Scilab GLPK

lines(0);

printf('LP problem\n');

c = [10, 6, 4]';
a = [1,  1, 1;...
     10, 4, 5;...
     2,  2, 6];
b = [100,600,300]';
ctype = "LLL";
lb    = [0,0,0]';
//ub    = []';
ub    = %inf*[1,1,1]';
vartype = "CCC";
// Output all GLPK messages on workspace
param = init_param();
param = add_param(param,'msglev',3);
param = add_param(param,'lpsolver',1);
// Set save options
param = add_param(param,'save',1);
param = add_param(param,'savefilename','SimpleLP');
param = add_param(param,'savefiletype','fixedmps');

[xmin,fmin,status,extra] = glpk(c,a,b,b,lb,ub,ctype,vartype,param);

printf('solution found: fmin = %f\n', fmin);disp(xmin);
printf('status = %d\n',status); disp(extra);

