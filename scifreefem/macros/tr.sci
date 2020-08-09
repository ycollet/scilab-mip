function v=tr(u)
  u     = u('value');
  [m,n] = size(u);

  if m<>n then error('square matrix expected'); end;

  v = '0';

  for i=1:m
    v=addf(v,u(i,i));
  end

  v = tlist(['ffeq','value'],v);
endfunction
