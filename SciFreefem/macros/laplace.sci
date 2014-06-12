function v=laplace(u)
  u     = u('value');
  [m,n] = size(u);

  if n>1 then error('vectoriel expression expected'); end;

  for l=1:m
    v(l) = 'laplace('+u(l)+')';
  end

  v = tlist(['ffeq','value'],v);
endfunction
