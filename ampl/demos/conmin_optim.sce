lines(0);

stacksize('max');

path = get_absolute_file_path('conmin_optim.sce');

exec(path + 'nl_data.sce');

Solve_macminlp = %F; // 43  files
Solve_coinor   = %F; // 24  files
Solve_asl      = %F; // 4   files
Solve_modnl    = %T; // 898 files

AddDeltaToX0 = %F;
Restart      = %F;

///////////////////////////
// Set global paremeters //
///////////////////////////

nl_index = 160;// 18 // Bug avec bonmin_optim.sce en premier + pb 18 de macminlp; 13, 14 17 (OK with ipopt NOK with bonmin): pb 
             // CoinOR + Index 13  = hs100 
             // ModNL: 52 - 62 - 462 - 489 ??
             // macminlp: 8 - 24 - 26 - [37 - 43]??

if Solve_macminlp then nl_filename = MacMINLP(nl_index); end
if Solve_coinor   then nl_filename = CoinOR(nl_index);   end
if Solve_asl      then nl_filename = ASL(nl_index);      end
if Solve_modnl    then nl_filename = ModNL(nl_index);    end

ItMX = 1000;

///////////////////////////////
// Load and test the problem //
///////////////////////////////

printf('\nOptimization of the %s problem.\n\n',basename(nl_filename));

[asl, x0, lower, upper, v, constr_lhs, constr_rhs] = ampl_init(nl_filename);

if AddDeltaToX0 then
  x0 = x0 + 0.1*rand(x0);
end

upper(find(upper==%inf))  =  1e6;
upper(find(upper==-%inf)) = -1e6;
lower(find(lower==%inf))  =  1e6;
lower(find(lower==-%inf)) = -1e6;

ncon = length(constr_rhs);

if ncon==0 then 
  printf('\nno constraints. Stops\n');
  return
end

////////////////////////////////////////////
// Dense objective and sparse constraints //
////////////////////////////////////////////

deff('[y,dy]=f(x)','[y,tmp] = ampl_evalf(asl,x); ...
                    [dy,tmp] = ampl_evalg(asl,x);');

deff('[y,dy,ic]=g(x,ct)','[tmp,y] = ampl_evalf(asl,x);...
                          [tmp,dy] = ampl_eval_sp_g(asl,x); ...
                          ic = find(y>ct); ...
                          dy = full(dy(:,ic))'';');

param = init_param();

param = add_param(param,'nacmx',ncon+2*length(x0)+1); // Number of constraints + side constraints + 1
param = add_param(param,'isc',zeros(size(x0,1),size(x0,2)));
param = add_param(param,'scal',ones(size(x0,1),size(x0,2)));
param = add_param(param,'nscal',0);
param = add_param(param,'nfdg',1);
param = add_param(param,'icndir',length(x0)+1);
param = add_param(param,'fdch',0.01);
param = add_param(param,'fdchm',0.01);
param = add_param(param,'ct',-0.1);
param = add_param(param,'ctmin',0.004);
param = add_param(param,'ctl',-0.01);
param = add_param(param,'ctlmin',0.001);
param = add_param(param,'theta',1.0);
param = add_param(param,'delfun',0.001);
param = add_param(param,'dabfun',0.001);
param = add_param(param,'linobj',0);
param = add_param(param,'itrm',3);
param = add_param(param,'alphax',0.3);
param = add_param(param,'abobj1',0.2);
param = add_param(param,'infog',1);
param = add_param(param,'info',1);
param = add_param(param,'iprint',2);

[x_opt,f_opt,df_opt,g_opt,dg_opt,ic_res] = conmin_optim(x0,f,g,ncon,upper,lower,ItMX,param);

if Restart then
  printf('Restart activated\n');
  
  printf('Initial values\n');
  g_res = g(x0,0);
  if isempty(g_res) then g_res = 0; end
  printf('Value of the objective function: %f\n',f(x0));
  printf('Maximum value of the constraints: %f\n',max(g_res));
  printf('Minimum value of the constraints: %f\n',min(g_res));

  printf('Final values\n');
  g_res = g(x_opt,0);
  if isempty(g_res) then g_res = 0; end
  printf('Value of the objective function: %f\n',f(x_opt));
  printf('Maximum value of the constraints: %f\n',max(g_res));
  printf('Minimum value of the constraints: %f\n',min(g_res));

  [x_opt,f_opt,df_opt,g_opt,dg_opt,ic_res] = conmin_optim(x_opt,f,g,ncon,upper,lower,ItMX,param);
end

printf('Initial values\n');
//printf('Value of the starting solution: ');
//disp(x0');
g_res = g(x0,0);
if isempty(g_res) then g_res = 0; end
printf('Value of the objective function: %f\n',f(x0));
printf('Maximum value of the constraints: %f\n',max(g_res));
printf('Minimum value of the constraints: %f\n',min(g_res));

printf('Final values\n');
//printf('Value of the optimal solution: ');
//disp(x_opt');
g_res = g(x_opt,0);
if isempty(g_res) then g_res = 0; end
printf('Value of the objective function: %f\n',f(x_opt));
printf('Maximum value of the constraints: %f\n',max(g_res));
printf('Minimum value of the constraints: %f\n',min(g_res));
//disp([string(lower) + ' <= ' + string(x_opt) + ' <= ' + string(upper)]);

