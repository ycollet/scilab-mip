function v = pde_sol(a,b,c)
  a = a('value');
  b = b('value');
  c = c('value');
  [ma,na] = size(a);
  [mb,nb] = size(b);
  [mc,nc] = size(c);

  if na<>1|nb<>1|nc<>1|mb<>ma|mc<>ma then
    error('arguments may be column vectors with equal sizes');
  end
  
  v = strcat('pde('+a+') '+b+'='+c+'; ')
endfunction
