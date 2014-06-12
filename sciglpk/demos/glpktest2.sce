// Two examples to show how to solve an MILP and an LP with interior method

printf('-- Integer problem --\n');

lines(0);

c = [-1,-1]';
a = [-2,5;2,-2];
b = [5;1];

ctype   = "LL";
lb      = [0;0];
ub      = [];
vartype = "BB";

param = init_param();
param = add_param(param,'msglev',3);
param = add_param(param,'presol',1);
param = add_param(param,'sense',1);

[xmin,fmin,status,extra] = glpk(c,a,b,b,lb,ub,ctype,vartype,param);

printf('solution found: fmin = %f\n', fmin);disp(xmin);
printf('status = %d\n',status); disp(extra);
pause;

printf('3rd problem\n');
s = 1;
c = [ 0 0 0 -1 -1]';
a = [-2 0 0  1  0;...
      0 1 0  0  2;...
      0 0 1  3  2];
b       = [4 12 18]';
ctype   = "GGG";
lb      = [0,0,0,0,0]';
ub      = 100*[1,1,1,1,1]';
vartype = "CCCCC";

[xmin,fmin,status,extra] = glpk(c,a,b,b,lb,ub,ctype,vartype,param);
printf('solution found: fmin = %f\n', fmin);disp(xmin);
printf('status = %d\n',status); disp(extra);

