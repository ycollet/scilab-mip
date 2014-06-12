lines(0);

path = get_absolute_file_path('nm_optim.sce');

getf(path+'/cont_funcs.sci');

tmp = getdate();
rand('seed',sum(tmp));

list_of_testpb = get_function_list();

// There are 67 test problems
Index = 1;

Dimension = 2;

Max = eval('max_bd_'+list_of_testpb(Index)+'(Dimension)');
Min = eval('min_bd_'+list_of_testpb(Index)+'(Dimension)');

deff('y=fobj(x)','y = ' + list_of_testpb(Index) + '(x);');

x0 = (Max-Min).*rand(Min,'uniform') + Min;

ItMX            = 1500;  // maximum number of descent steps
TOL             = 1e-6; // accuracy for convergence test - derivatives. 
                        // The tolerance is tested on the absolute variation of objective function found on the simplex.
                        // if abs(f_Good - f_Worth)<TOL then stop
Log             = %T;   // print info from line step, 1 = on / 0 = off
MaxEvalFunc     = 1500;
KelleyRestart   = %F;
KelleyAlpha     = 1e-7;
SimplexRelSize  = 0.1;
Plot            = %T;

// We set the initial simplex

x_init(:,1) = x0 + SimplexRelSize*0.5*((Max - Min).*grand(x0,'unf',0,1) + Min);
x_init(:,2) = x0 + SimplexRelSize*0.5*((Max - Min).*grand(x0,'unf',0,1) + Min);
x_init(:,3) = x0 + SimplexRelSize*0.5*((Max - Min).*grand(x0,'unf',0,1) + Min);

disp(x_init)

//
// Nelder and Mead optimization
//

printf('Test of the Nelder and Mead optimization method.\n');
[x_opt, x_hist] = optim_nelder_mead(fobj, x_init, ItMX, TOL, MaxEvalFunc, Log, KelleyRestart, KelleyAlpha);

printf('Number of function evaluation: %d\n',length(x_hist));

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
  
  xtitle('Evolution of the classical simplex','x1','x2');
  
  for i=1:length(x_hist)
    plot(x_hist(i)(1)(1), x_hist(i)(1)(2), 'ko');
    plot(x_hist(i)(2)(1), x_hist(i)(2)(2), 'ko');
    plot(x_hist(i)(3)(1), x_hist(i)(3)(2), 'ko');
    plot([x_hist(i)(1)(1) x_hist(i)(2)(1)], [x_hist(i)(1)(2) x_hist(i)(2)(2)], 'k-');
    plot([x_hist(i)(2)(1) x_hist(i)(3)(1)], [x_hist(i)(2)(2) x_hist(i)(3)(2)], 'k-');
    plot([x_hist(i)(3)(1) x_hist(i)(1)(1)], [x_hist(i)(3)(2) x_hist(i)(1)(2)], 'k-');
  end

  subplot(1,2,2);
  F = [];
  for i=1:length(x_hist)
    F(i) = fobj(x_hist(i)(1));
  end
  plot(F)
  xtitle('Evolution of the objective function','Iteration','FObj');
  drawnow;
end

//
// Nelder and Mead optimization
//

rand('seed',sum(tmp));

printf('Test of the Nelder and Mead optimization method. Automatic initialization\n');

[x_opt, x_hist] = optim_nelder_mead(fobj, x0, ItMX, TOL, MaxEvalFunc, Log, KelleyRestart, KelleyAlpha, SimplexRelSize*(Max-Min));

printf('Number of function evaluation: %d\n',length(x_hist));

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
  
  xtitle('Evolution of the classical-simplex','x1','x2');
  
  for i=1:length(x_hist)
    plot(x_hist(i)(1)(1), x_hist(i)(1)(2), 'ko');
    plot(x_hist(i)(2)(1), x_hist(i)(2)(2), 'ko');
    plot(x_hist(i)(3)(1), x_hist(i)(3)(2), 'ko');
    plot([x_hist(i)(1)(1) x_hist(i)(2)(1)], [x_hist(i)(1)(2) x_hist(i)(2)(2)], 'k-');
    plot([x_hist(i)(2)(1) x_hist(i)(3)(1)], [x_hist(i)(2)(2) x_hist(i)(3)(2)], 'k-');
    plot([x_hist(i)(3)(1) x_hist(i)(1)(1)], [x_hist(i)(3)(2) x_hist(i)(1)(2)], 'k-');
  end

  subplot(1,2,2);
  F = [];
  for i=1:length(x_hist)
    F(i) = fobj(x_hist(i)(1));
  end
  plot(F)
  xtitle('Evolution of the objective function','Iteration','FObj');
  drawnow;
end

