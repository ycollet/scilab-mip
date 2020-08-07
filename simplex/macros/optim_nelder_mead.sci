function [x_opt, x_history] = optim_nelder_mead(onm_f, x0, ItMX, Tol, MaxEvalFunc, Log, kelley_restart, kelley_alpha, rel_size)

[nargout, nargin] = argn();

// x0:
// We must store a n+1 x n matrix which handles the initial simplex
// We can store a n x 1 vector. In this case, an initialization function will be called which computes an initial simplex based upon the rel_size vector
// rel_size:
// optional parameter which handles a vector of a box sizes which will contains the initial simplex. The initial simplex will be randomly generated in this box.

nm_alpha = 1;
nm_gamma = 2;
nm_rho   = 0.5;
nm_sigma = 0.5;

if size(x0,2)==size(x0,1)+1 then
  n = size(x0,1);
else
  n = length(x0);
end

if ~isdef('rel_size','local') then
  rel_size = 0.1*ones(size(x0,1),size(x0,2));
end

if min(size(x0))==1 then
  // We call the init function
  x0 = nm_init(x0,rel_size);
end

x_history_defined = (nargout==2);

if (x_history_defined) then
  x_history = list();
  x_history($+1) = list();
  for i=1:n+1
    x_history($)($+1) = x0(:,i);
  end
end

if (~isdef('ItMX','local')) then
  ItMX = 100;
end
if (~isdef('Tol','local')) then
  Tol = 0.0;
end
if (~isdef('MaxEvalFunc','local')) then
  MaxEvalFunc = 10*ItMX;
end
if (~isdef('kelley_restart','local')) then
  kelley_restart = %F;
end
if (~isdef('kelley_alpha','local')) then
  kelley_alpha = 1e-4;
end
if (~isdef('Log','local')) then
  Log = %F;
end

if ~isdef('onm_f','local') then
  error('optim_nelder_mead: onm_f is mandatory');
else
  if typeof(onm_f)=='list' then
    deff('y=_nm_f(x)','y=onm_f(1)(x, onm_f(2:$))');
  else
    deff('y=_nm_f(x)','y=onm_f(x)');
  end
end

// We compute and sort the objective function values of these 3 points

x_nm = list();

for i=1:n+1
  f_nm(i)   = _nm_f(x0(:,i));
  x_nm($+1) = x0(:,i);
end


[f_sort, I] = gsort(f_nm);
I_Worth     = 1;
I_2nd_Worth = 2;
I_Best      = length(f_nm);

if (Log) then
  printf('optim_nelder_mead: initialization\n');
  printf('f_Worth = %f f_Best = %f\n', f_nm(I(I_Worth)), f_nm(I(I_Best)));
end

FirstIteration = %T;
Iteration = 0;
eval_Func = 0;

// The main loop
while(((eval_Func < MaxEvalFunc) & (Iteration<ItMX) & (abs((f_nm(I(I_Best)) - f_nm(I(I_Worth))))>Tol)) | (FirstIteration))
  if (Log) then
    printf('optim_nelder_mead: test eval_Func < MaxEvalFunc: %d < %d\n',eval_Func, MaxEvalFunc);
    printf('optim_nelder_mead: Iteration < ItMX: %d < %d\n', Iteration, ItMX);
    printf('optim_nelder_mead: abs(f_Best - f_Worth)>Tol: %f > %f\n', abs(f_nm(I(I_Best)) - f_nm(I(I_Worth))),Tol);
  end
      
  f_Old          = f_nm(I(I_Best));
  FirstIteration = %F;
  Iteration = Iteration + 1;
  x_Median  = zeros(x_nm(1));
  for i=2:length(x_nm)
    x_Median = x_Median + x_nm(I(i));
  end
  x_Median = x_Median / (length(x_nm) - 1);
  
  // First, we reflect the worth point
  if (Log) then
    printf('optim_nelder_mead: iteration = %d - reflexion\n', Iteration);
  end
  
  x_Reflexion = x_Median + nm_alpha*(x_Median - x_nm(I(I_Worth)));  
  f_Reflexion = _nm_f(x_Reflexion);
  eval_Func   = eval_Func + 1;
    
  if (f_Reflexion<f_nm(I(I_2nd_Worth))) then
    // Case I: Reflexion or expansion
    if (Log) then
      printf('optim_nelder_mead: reflexion or expansion - iteration = %d\n', Iteration);
    end
    if (f_nm(I(I_Best)) <= f_Reflexion) then
      f_nm(I(I_Worth)) = f_Reflexion;
      x_nm(I(I_Worth)) = x_Reflexion;
    elseif (f_Reflexion < f_nm(I(I_Best))) then
      if (Log) then
        printf('optim_nelder_mead: expansion\n');
      end
      x_Expansion = x_Median + nm_gamma*(x_Median - x_nm(I(I_Worth)));
      f_Expansion = _nm_f(x_Expansion);
      eval_Func   = eval_Func + 1;
      if (f_Expansion<f_Reflexion) then
        f_nm(I(I_Worth)) = f_Expansion;
        x_nm(I(I_Worth)) = x_Expansion;
      else
        f_nm(I(I_Worth)) = f_Reflexion;
        x_nm(I(I_Worth)) = x_Reflexion;
      end  
    end
  else
    // Case II: Expansion or Shrink
    if (Log) then
      printf('optim_nelder_mead: contraction or shrink - iteration = %d\n', Iteration);
    end

    if (Log) then
      printf('optim_nelder_mead: contraction\n');
    end
    x_Contr1  = x_nm(I(I_Worth)) + nm_rho*(x_Median - x_nm(I(I_Worth)));
    f_Contr1  = _nm_f(x_Contr1);
    eval_Func = eval_Func + 1;
    x_Contr2  = x_Median - nm_rho*(x_Median - x_nm(I(I_Worth)));
    f_Contr2  = _nm_f(x_Contr2);
    eval_Func = eval_Func + 1;
    if (f_Contr1<=f_nm(I(I_Worth))) then
      f_nm(I(I_Worth)) = f_Contr1;
      x_nm(I(I_Worth)) = x_Contr1;
    elseif (f_Contr2<=f_nm(I(I_Worth))) then
      f_nm(I(I_Worth)) = f_Contr2;
      x_nm(I(I_Worth)) = x_Contr2;
    else
      if (Log) then
        printf('optim_nelder_mead: shrink - iteration = %d\n', Iteration);
      end
      for i=2:length(x_nm)
        x_nm(I(i)) = x_nm(I(I_Best)) + nm_rho*(x_nm(I(i)) - x_nm(I(I_Best)));
        f_nm(I(i)) = _nm_f(x_nm(I(i)));
        eval_Func = eval_Func + 1;
      end
    end
  end

  // Ranking the points  
  [f_sort, I] = gsort(f_nm);
  
  if (kelley_restart) then
    // Kelley restart
    // Computation of the simplex gradient
    df = f_nm(I(1:$-1)) - f_nm(I(I_Best));
    for i=1:length(x_nm)-1
      V     = x_nm(I(i)) - x_nm(I(I_Best));
      Df(i) = V'*df;
    end
    // Test
    if (f_nm(I(I_Best)) - f_Old > -kelley_alpha * norm(Df)^2) then
      if (Log) then
        printf('optim_nelder_mead: Kelley restart\n');
      end
      // Simplex restart
      for i=1:n-1
        Aux(i) = norm(x_nm(I(I_Best))-x_nm(I(i)))
      end
      sigma_m_V = min(Aux);
      _beta = 0.5*sigma_m_V*sign(Df);
      y_nm = list();
      y_nm(1) = x_nm(1);
      for i=1:length(x_nm(I(I_Best)))
        y_nm(i+1)    = x_nm(I(I_Best));
        y_nm(i+1)(i) = y_nm(i+1)(i) + _beta(i);
      end
      x_nm = list();
      for i=1:length(y_nm)
        x_nm(i)   = y_nm(i);
        f_nm(i)   = _nm_f(x_nm(i));
        eval_Func = eval_Func + 1;
      end
      y_nm = list();
    end
  end
  
  if (Log) then
    printf('optim_nelder_mead: f_Worth = %f f_Best = %f\n', f_nm(I(I_Worth)), f_nm(I(I_Best)));
  end

  if (x_history_defined) then
    x_history($+1) = list();
    for i=1:n+1
      x_history($)($+1) = x_nm(i);
    end
  end
end

if (Log) then
  printf('optim_nelder_mead: test evall_Func < MaxEvalFunc: %d < %d\n',eval_Func, MaxEvalFunc);
  printf('optim_nelder_mead: Iteration < ItMX: %d < %d\n', Iteration, ItMX);
  printf('optim_nelder_mead: abs(f_Best - f_Worth)>Tol: %f > %f\n', abs(f_nm(I(I_Best)) - f_nm(I(I_Worth))),Tol);
end
    
x_opt = x_nm(I(I_Best));
endfunction
