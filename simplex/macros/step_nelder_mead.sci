function [x_next, data_next, eval_Func, f_hist, x_hist] = step_nelder_mead(f_current, x_current, data_current, nm_mode, Log, kelley_restart, kelley_alpha, rel_size)

[nargout, nargin] = argn();

if typeof(x_current)=='list' then
  n = length(x_current)-1;
  init_simplex_good = %F;
else
  if size(x_current,2)==size(x_current,1)+1 then
    n = size(x_current,1);
    init_simplex_good = %T;
  else
    n = length(x_current);
    init_simplex_good = %F;
  end
end

if nm_mode=='build_simplex' then
  if ~isdef('rel_size','local') then
    rel_size = 0.1*ones(size(x_current,1),size(x_current,2));
  end

  // We call the init function
  x_current = nm_init(x_current,rel_size);
  x_next    = x_current;
  eval_Func = 0;
  f_hist = [];
  x_hist = [];
  return;
end

if (~isdef('Log','local')) then
  Log = %F;
end
if (~isdef('kelley_restart','local')) then
  kelley_restart = %F;
end
if (~isdef('kelley_alpha','local')) then
  kelley_alpha = 1e-4;
end

f_hist_is_defined = (nargout>=4);
x_hist_is_defined = (nargout==5);

Compute = %F;
// Initialisation of the mlist data_next structure

data_next = mlist(['psimplex','f_Worth','f_Best','x_Worth','x_Best','Case','f_Reflexion','x_Reflexion','f_Expansion','x_Expansion', ...
                  'f_C1','x_C1','f_C2','x_C2','x_Median','f_Hist','x_Hist','eval_Func', 'f_Old', 'I', 'Iteration', ...
                  'f_Simplex','x_Simplex']);

nm_alpha = 1;
nm_gamma = 2;
nm_rho   = 0.5;
nm_sigma = 0.5;

I_Best      = n+1;
I_Worth     = 1;
I_2nd_Worth = 2;

if length(data_current('x_Simplex')(1))==1 then
  printf('ERROR: Case = %d\n', data_current('Case'));
end

while ~Compute
  ///////////////////
  // 'init' option //
  ///////////////////
  if (nm_mode=='init') then
    n = size(x_current,1);
    if (size(x_current,2)~=n+1) then
      error('nelder_mead: you must give %d starting points\n',n+1);
    end
    // Initialisation of the state structure
    data_next('f_Simplex') = [];
    data_next('x_Simplex') = list();
    for i=1:n+1
      data_next('f_Simplex')(i) = f_current(i);
      data_next('x_Simplex')(i) = x_current(:,i);
    end

    // We order the values of the simplex
    [f_sort, data_next('I')] = gsort(f_current);
    
    data_next('f_Worth')     = f_current(data_next('I')(I_Worth));
    data_next('f_Best')      = f_current(data_next('I')(I_Best));
    data_next('x_Worth')     = x_current(:,data_next('I')(I_Worth));
    data_next('x_Best')      = x_current(:,data_next('I')(I_Best));
    data_next('f_Old')       = data_next('f_Best');
    data_next('Case')        = 1;
    data_next('f_Reflexion') = f_current(data_next('I')(I_Best)); 
    data_next('x_Reflexion') = x_current(:,data_next('I')(I_Best)); 
    data_next('f_Expansion') = f_current(data_next('I')(I_Best)); 
    data_next('x_Expansion') = x_current(:,data_next('I')(I_Best)); 
    data_next('f_C1')        = f_current(data_next('I')(I_Best)); 
    data_next('x_C1')        = x_current(:,data_next('I')(I_Best)); 
    data_next('f_C2')        = f_current(data_next('I')(I_Best)); 
    data_next('x_C2')        = x_current(:,data_next('I')(I_Best)); 
    data_next('x_Median')    = x_current(:,data_next('I')(I_Best)); 
    data_next('f_Hist')      = list(); 
    data_next('x_Hist')      = list(); 
    data_next('eval_Func')   = 0; 
    data_next('Iteration')   = 0; 
        
    // We order the values of the simplex
    f_sort = [];
    for i=1:length(data_next('f_Simplex'))
      f_sort(i) = data_next('f_Simplex')(i);
    end
    [f_sort, data_next('I')] = gsort(f_sort);
    
    if (x_hist_is_defined) then
      data_next('x_Hist')($+1)  = data_next('x_Simplex');
    end

    if (f_hist_is_defined) then
      data_next('f_Hist')($+1)  = data_next('f_Simplex');
    end

    if (Log) then
      printf('step_nelder_mead: init\n');
      printf('f_Worth = %f f_Best = %f\n', data_next('f_Worth'), data_next('f_Best'));
    end
    
    data_next('Case') = 1; // Starting position of the algorithm
  end

  //////////////////
  // 'run' option //
  //////////////////
  if (nm_mode=='run') then
    if typeof(data_current)~='psimplex' then
      error('step_nelder_mead: data_current must be a psimplex structure');
    end
  
    data_next = data_current;
        
    data_next('f_Worth') = data_current('f_Simplex')(data_next('I')(I_Worth));
    data_next('f_Best')  = data_current('f_Simplex')(data_next('I')(I_Best));
    data_next('x_Worth') = data_current('x_Simplex')(data_next('I')(I_Worth));
    data_next('x_Best')  = data_current('x_Simplex')(data_next('I')(I_Best));
  end

  ///////////////////
  // 'exit' option //
  ///////////////////  
  if (nm_mode=='exit') then
    x_next    = data_current('x_Best');
    data_next = data_current('f_Best');
    f_hist    = data_current('f_Hist');
    x_hist    = data_current('x_Hist');
    return;
  end
      
  if (Log) then
    printf('step_nelder_mead: run - Case = %d\n', data_next('Case'));
    printf('f_Worth = %f f_Best = %f\n', data_next('f_Worth'), data_next('f_Best'));
  end

  select(data_next('Case'))
  case 1 then
    // Computation of R
    if (Log) then
      printf('step_nelder_mead: Reflexion\n');
    end
    
    // We order the values of the simplex
    for i=1:length(data_next('f_Simplex'))
      f_sort(i) = data_next('f_Simplex')(i);
    end
    [f_sort, data_next('I')] = gsort(f_sort);
    
    data_next('f_Old')     = data_next('f_Best');
    data_next('Iteration') = data_next('Iteration') + 1; 
        
    data_next('x_Median')  = zeros(size(data_next('x_Best'),1),size(data_next('x_Best'),2));
    for i=2:length(data_next('x_Simplex'))
      data_next('x_Median') = data_next('x_Median') + data_next('x_Simplex')(data_next('I')(i));
    end
    data_next('x_Median') = data_next('x_Median') / (length(data_next('x_Simplex')) - 1);
    
    data_next('x_Reflexion') = data_next('x_Median') + nm_alpha*(data_next('x_Median') - data_next('x_Simplex')(data_next('I')(I_Worth))); 
    x_next                   = data_next('x_Reflexion');
    data_next('Case')        = 2; // Go to block 2
    Compute                  = %T; // Ask the user for an objective function value
  case 2 then
    // We get f(R)
    data_next('f_Reflexion') = f_current;
    data_next('Case')        = 3; // Go to block 3
    data_next('eval_Func')   = data_next('eval_Func') + 1;
    Compute                  = %F;
  case 3 then
    if (data_next('f_Reflexion')<data_next('f_Simplex')(data_next('I')(I_2nd_Worth))) then
      data_next('Case') = 4; // Go to block 4
    else
      data_next('Case') = 13; // Go to block 13
    end
    Compute = %F;
  case 4 then
    if (data_next('f_Best')<=data_next('f_Reflexion')) then
      data_next('Case') = 5; // Go to block 5
    else
      data_next('Case') = 6; // Go to block 6
    end
    Compute = %F;
  case 5 then
    data_next('f_Simplex')(data_next('I')(I_Worth)) = data_next('f_Reflexion');
    data_next('x_Simplex')(data_next('I')(I_Worth)) = data_next('x_Reflexion');
    data_next('Case')     = 22; // Next Case (we start again), but before, we rank the points
    Compute = %F;
  case 6 then
    // Computation of E
    if (Log) then
      printf('step_nelder_mead: Expansion\n');
    end
    if (data_next('f_Reflexion')<data_next('f_Best')) then
      data_next('x_Expansion') = data_next('x_Median') + nm_gamma*(data_next('x_Median') - data_next('x_Simplex')(data_next('I')(I_Worth)));
      x_next                   = data_next('x_Expansion');
      data_next('Case')        = 7; // Go to block 7
      Compute                  = %T;
    else
      data_next('Case')        = 22; // Go to block 22
      Compute                  = %F;
    end      
  case 7 then
    // Computation of f(E)
    data_next('f_Expansion') = f_current;
    data_next('Case')        = 8; // Go to block 8
    data_next('eval_Func')   = data_next('eval_Func') + 1;
    Compute                  = %F;
  case 8 then
    if (data_next('f_Expansion')<data_next('f_Reflexion')) then
      data_next('Case') = 9; // Go to block 9
    else
      data_next('Case') = 10; // Go to block 10
    end
    Compute = %F;
  case 9 then
    // We replace W by E
    data_next('x_Simplex')(data_next('I')(I_Worth)) = data_next('x_Expansion');
    data_next('f_Simplex')(data_next('I')(I_Worth)) = data_next('f_Expansion');
    data_next('Case')    = 22; // Next Case: we start again but before, we rank the points
    Compute              = %F;
  case 10 then
    // We replace W by R
    data_next('x_Simplex')(data_next('I')(I_Worth)) = data_next('x_Reflexion');
    data_next('f_Simplex')(data_next('I')(I_Worth)) = data_next('f_Reflexion');
    data_next('Case')    = 22; // Go to block 22
    Compute              = %F;
  case 13 then
    // Computation of C1
    if (Log) then
      printf('step_nelder_mead: Contraction\n');
    end
    data_next('x_C1') = data_next('x_Simplex')(data_next('I')(I_Worth)) + nm_rho*(data_next('x_Median') - data_next('x_Simplex')(data_next('I')(I_Worth)));
    x_next            = data_next('x_C1');
    data_next('Case') = 14; // Go to block 14
    Compute = %T;
  case 14 then
    // Computation of f(C1)
    data_next('f_C1')      = f_current;
    data_next('Case')      = 15; // Go to block 15
    data_next('eval_Func') = data_next('eval_Func') + 1;
    Compute                = %F;
  case 15 then
    // Computation of C2
    data_next('x_C2') = data_next('x_Median') - nm_rho*(data_next('x_Median') - data_next('x_Simplex')(data_next('I')(I_Worth)));
    x_next            = data_next('x_C2');
    data_next('Case') = 16; // Go to block 16
    Compute           = %T;
  case 16 then
    // Computation of f(C2)
    data_next('f_C2')      = f_current;
    data_next('Case')      = 17;
    data_next('eval_Func') = data_next('eval_Func') + 1;
    Compute                = %F;
  case 17 then
    if (data_next('f_C1')<=data_next('f_Worth')) then
      data_next('f_Simplex')(data_next('I')(I_Worth)) = data_next('f_C1');
      data_next('x_Simplex')(data_next('I')(I_Worth)) = data_next('x_C1');
      data_next('Case')    = 22; // Go to block 22
    elseif (data_next('f_C2')<=data_next('f_Worth')) then
      data_next('f_Simplex')(data_next('I')(I_Worth)) = data_next('f_C2');
      data_next('x_Simplex')(data_next('I')(I_Worth)) = data_next('x_C2');
      data_next('Case')    = 22; // Go to block 22
    else
      data_next('Case') = 18; // Go to block 18
    end
    Compute = %F;
  case 18 then
    // Computation of S
    if (Log) then
      printf('step_nelder_mead: Shrink\n');
    end
    x_next = list();
    for i=2:length(data_next('x_Simplex'))
      x_next(i-1) = data_next('x_Simplex')(data_next('I')(I_Best)) + nm_rho*(data_next('x_Simplex')(data_next('I')(i)) - data_next('x_Simplex')(data_next('I')(I_Best)));
    end
        
    data_next('Case')     = 19; // Go to block 19
    Compute               = %T;
  case 19 then
    // Computation of f(S)
    for i=1:length(f_current)
      data_next('f_Simplex')(i) = f_current(i);
    end
    data_next('Case')      = 22; // Go to block 20
    data_next('eval_Func') = data_next('eval_Func') + length(f_current);
    Compute                = %F;
  case 22 then
    // Ranking the points
    for i=1:length(data_next('f_Simplex'))
      f_sort(i) = data_next('f_Simplex')(i);
    end
    [f_sort, data_next('I')] = gsort(f_sort);

    data_next('f_Worth') = data_next('f_Simplex')(data_next('I')(I_Worth));
    data_next('f_Best')  = data_next('f_Simplex')(data_next('I')(I_Best));

    // Updating the simplex structure
    data_next('x_Worth') = data_next('x_Simplex')(data_next('I')(I_Worth));
    data_next('x_Best')  = data_next('x_Simplex')(data_next('I')(I_Best));
    
    if (x_hist_is_defined) then
      data_next('x_Hist')($+1)  = data_next('x_Simplex');
    end

    if (f_hist_is_defined) then
      data_next('f_Hist')($+1)  = data_next('f_Simplex');
    end

    data_next('Case') = 23; // Go to block 23
    Compute           = %F;
  case 23 then
    if (kelley_restart) then
      // Kelley restart
      // Computation of the simplex gradient
      for i=1:length(data_next('f_Simplex'))-1
        df(i) = data_next('f_Simplex')(data_next('I')(i)) - data_next('f_Simplex')(data_next('I')(I_Best));
      end
      for i=1:length(data_next('x_Simplex'))-1
        V     = data_next('x_Simplex')(data_next('I')(i)) - data_next('x_Simplex')(data_next('I')(I_Best));
        Df(i) = V'*df;
      end
      // Test
      if (data_next('f_Simplex')(data_next('I')(I_Best)) - data_next('f_Old') > -kelley_alpha * norm(Df)^2) then
        if (Log) then
          printf('step_nelder_mead: Kelley restart\n');
        end
        // Simplex restart
        Aux = [];
        for i=1:length(data_next('f_Simplex'))-1
          Aux(i) = norm(data_next('x_Simplex')(data_next('I')(I_Best)) - data_next('x_Simplex')(data_next('I')(i)))
        end
        sigma_m_V = min(Aux);
        _beta     = 0.5*sigma_m_V*sign(Df);
        x_next    = list();
        x_next(1) = data_next('x_Simplex')(data_next('I')(1));
        for i=1:length(data_next('x_Simplex')(data_next('I')(I_Best)))
          x_next(i+1)    = data_next('x_Simplex')(data_next('I')(I_Best));
          x_next(i+1)(i) = x_next(i+1)(i) + _beta(i);
        end
      end
    end
    data_next('Case') = 24; // Go to block 24
    Compute           = %T;
  case 24 then
    for i=1:length(f_current)
      data_next('f_Simplex')(i) = f_current(i)
    end
    if typeof(x_current)=='list' then
      for i=1:length(f_current)
        data_next('x_Simplex')(i) = x_current(i)
      end
    else
      for i=1:length(f_current)
        data_next('x_Simplex')(i) = x_current(:,i)
      end
    end
    data_next('eval_Func') = data_next('eval_Func') + length(f_current);
    data_next('Case')      = 1; // Go back to block 1
    Compute                = %F;
  end
  
  data_current = data_next; // when in the while loop, we update the data_current
                            // structure to update the state of the Nelder and Mead
                            // method
end

if (f_hist_is_defined) then
  f_hist = data_next('f_Hist');
end
if (x_hist_is_defined) then
  x_hist = data_next('x_Hist');
end
eval_Func = data_next('eval_Func');
endfunction
