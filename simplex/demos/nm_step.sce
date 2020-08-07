lines(0);

path = get_absolute_file_path('nm_step.sce');

getf(path+'/cont_funcs.sci');

list_of_testpb = get_function_list();

tmp = getdate();
rand('seed',sum(tmp));

// There are 67 test problems
Index = 1;

Dimension = 2;

Max = eval('max_bd_' + list_of_testpb(Index) + '(Dimension);');
Min = eval('min_bd_' + list_of_testpb(Index) + '(Dimension);');
deff('y=fobj(x)','y = ' + list_of_testpb(Index) + '(x);');

x0 = (Max-Min).*rand(size(Min,1),size(Min,2)) + Min;

ItMX            = 1500;  // maximum number of descent steps
TOL             = 1e-6; // accuracy for convergence test - derivatives. 
                        // The tolerance is tested on the absolute variation of objective function found on the simplex.
                        // if abs(f_Good - f_Worth)<TOL then stop
Log             = %T;   // print info from line step, 1 = on / 0 = off
MaxEvalFunc     = 1500;
KelleyRestart   = %T;
KelleyAlpha     = 1e-7;
SimplexRelSize  = 0.1;
Plot            = %T;

// We set the initial simplex

x_init(:,1) = x0 + SimplexRelSize*0.5*((Max - Min).*grand(size(x0,1),size(x0,2),'unf',0,1) + Min);
x_init(:,2) = x0 + SimplexRelSize*0.5*((Max - Min).*grand(size(x0,1),size(x0,2),'unf',0,1) + Min);
x_init(:,3) = x0 + SimplexRelSize*0.5*((Max - Min).*grand(size(x0,1),size(x0,2),'unf',0,1) + Min);

printf('The initial simplex:\n');

disp(x_init)

//
// We start the Step Nelder Mead algorithm
//

printf('Test of the Step Nelder and Mead optimization method.\n');

x_hist = list();
f_hist = list();

f_init(1) = fobj(x_init(:,1));
f_init(2) = fobj(x_init(:,2));
f_init(3) = fobj(x_init(:,3));

[x_next, data_next, eval_Func, f_hist, x_hist] = step_nelder_mead(f_init, x_init, [], 'init', Log, KelleyRestart, KelleyAlpha);

f_current = fobj(x_next);

if Log then printf('step_nelder_mead - Initial iteration: f = %f\n', f_current); end

while eval_Func<MaxEvalFunc & (abs((max(f_hist($)(:)) - min(f_hist($)(:))))>TOL)
  [x_next, data_next, eval_Func, f_hist, x_hist] = step_nelder_mead(f_current, x_next, data_next, 'run', Log, KelleyRestart, KelleyAlpha);

  if typeof(x_next)=='list' then
    f_current = [];
    for i=1:length(x_next)
      f_current(i) = fobj(x_next(i));
    end
  else
    f_current = fobj(x_next);
  end
  
  if Log then printf('step_nelder_mead - Function evaluation %d, Iteration %d: f = %f\n', eval_Func, data_next('Iteration'), f_current); end
end

[x_best, f_best, eval_Func, f_hist, x_hist] = step_nelder_mead(f_current, x_next, data_next, 'exit', Log, KelleyRestart, KelleyAlpha);

printf('step_nelder_mead: best value found: %f\n', f_best);
printf('step_nelder_mead: nb of function evaluation: %d\n', eval_Func);
printf('step_nelder_mead: nb of iterations: %d\n', data_next('Iteration'));
  
if (Plot) then
  scf;
  drawlater;

  subplot(1,2,1);
  xgrid(2);
  
  X = Min(1):(Max(1)-Min(1))/10:Max(1);
  Y = Min(2):(Max(2)-Min(2))/10:Max(2);
  for i=1:length(X)
    for j=1:length(Y)
      Z(i,j) = fobj([X(i),Y(j)]);
    end
  end
  xset('fpf',' ');
  contour(X,Y,Z,10);
  
  xtitle('Evolution of the step-simplex','x1','x2');
  
  for i=1:length(x_hist)
    plot(x_hist(i)(1)(1), x_hist(i)(1)(2), 'ko');
    plot(x_hist(i)(2)(1), x_hist(i)(2)(2), 'ko');
    plot(x_hist(i)(3)(1), x_hist(i)(3)(2), 'ko');
    plot([x_hist(i)(1)(1) x_hist(i)(2)(1)], [x_hist(i)(1)(2) x_hist(i)(2)(2)], 'k-');
    plot([x_hist(i)(2)(1) x_hist(i)(3)(1)], [x_hist(i)(2)(2) x_hist(i)(3)(2)], 'k-');
    plot([x_hist(i)(3)(1) x_hist(i)(1)(1)], [x_hist(i)(3)(2) x_hist(i)(1)(2)], 'k-');
  end

  subplot(1,2,2);
  for i=1:length(f_hist)
    F(i) = f_hist(i)(1);
  end
  plot(F)
  xtitle('Evolution of the objective function','Iteration','FObj');
  drawnow;
end


//
// Plot again a surface of the search domain
//

//
// We start the Step Nelder Mead algorithm
//

printf('Test of the Step Nelder and Mead optimization method. Automatic initialization\n');

x_hist = list();
f_hist = list();

f_init = [];
x_init = x0;

[x_next, data_next, eval_Func, f_hist, x_hist] = step_nelder_mead(f_init, x_init, [], 'build_simplex', Log, KelleyRestart, KelleyAlpha, SimplexRelSize);

for i=1:size(x_next,2)
  f_init(i) = fobj(x_next(:,i));
end
x_init = x_next;

[x_next, data_next, eval_Func, f_hist, x_hist] = step_nelder_mead(f_init, x_init, [], 'init', Log, KelleyRestart, KelleyAlpha, SimplexRelSize);
f_current = fobj(x_next);

if Log then printf('step_nelder_mead - Initial iteration: f = %f\n', f_current(1)); end
while eval_Func<MaxEvalFunc & (abs((max(f_hist($)(:)) - min(f_hist($)(:))))>TOL)
  [x_next, data_next, eval_Func, f_hist, x_hist] = step_nelder_mead(f_current, x_next, data_next, 'run', Log, KelleyRestart, KelleyAlpha);
  f_current = fobj(x_next);
  if Log then printf('step_nelder_mead - Function evaluation %d, Iteration %d: f = %f\n', eval_Func, data_next('Iteration'), f_current); end
end
[x_best, f_best, eval_Func, f_hist, x_hist] = step_nelder_mead(f_current, x_next, data_next, 'exit', Log, KelleyRestart, KelleyAlpha);

printf('step_nelder_mead: best value found: %f\n', f_best);
printf('step_nelder_mead: nb of function evaluation: %d\n', eval_Func);
printf('step_nelder_mead: nb of iterations: %d\n', data_next('Iteration'));
  
if (Plot) then
  scf;
  drawlater;

  subplot(1,2,1);
  xgrid(2);
  
  X = Min(1):(Max(1)-Min(1))/10:Max(1);
  Y = Min(2):(Max(2)-Min(2))/10:Max(2);
  for i=1:length(X)
    for j=1:length(Y)
      Z(i,j) = fobj([X(i),Y(j)]);
    end
  end
  xset('fpf',' ');
  contour(X,Y,Z,10);
  
  xtitle('Evolution of the step-simplex','x1','x2');
  
  for i=1:length(x_hist)
    plot(x_hist(i)(1)(1), x_hist(i)(1)(2), 'ko');
    plot(x_hist(i)(2)(1), x_hist(i)(2)(2), 'ko');
    plot(x_hist(i)(3)(1), x_hist(i)(3)(2), 'ko');
    plot([x_hist(i)(1)(1) x_hist(i)(2)(1)], [x_hist(i)(1)(2) x_hist(i)(2)(2)], 'k-');
    plot([x_hist(i)(2)(1) x_hist(i)(3)(1)], [x_hist(i)(2)(2) x_hist(i)(3)(2)], 'k-');
    plot([x_hist(i)(3)(1) x_hist(i)(1)(1)], [x_hist(i)(3)(2) x_hist(i)(1)(2)], 'k-');
  end

  subplot(1,2,2);
  for i=1:length(f_hist)
    F(i) = f_hist(i)(1);
  end
  plot(F)
  xtitle('Evolution of the objective function','Iteration','FObj');
  drawnow;
end


