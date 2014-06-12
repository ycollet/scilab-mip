function v=grad(u)
  u     = u('value');
  [m,n] = size(u);
  v     = [];
  t     = ['x' 'y'];

  if (m*n <> 1) then
    if (n > 1) then error('scalar or vector expected'); end;

    for l=1:m*n
      for k=1:2
        v(l,k) = derive(t(k),u(l));
      end
    end
  else
    for k=1:2
      v(k) = derive(t(k),u(1));
    end
  end

  v = tlist(['ffeq','value'],v);
endfunction
