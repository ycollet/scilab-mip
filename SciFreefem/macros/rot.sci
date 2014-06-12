function v=rot(u)
  u     = u('value');
  [m,n] = size(u);

  if (m*n == 1) then error('scalar not accepted'); end;
  if (n > 1) then error('sorry, rot on matrix not implemented'); end;

  v   = [];
  t   = ['x' 'y'];
  ind = [2 1];
  moins = [' ' '-'];

  for k=1:m
    v(k) = stripblanks(moins(k)+derive(t(k),u(ind(k))));
  end

  v = tlist(['ffeq','value'],v);
endfunction
