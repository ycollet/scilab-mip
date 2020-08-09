function v = div(u)
  u     = u('value');
  [m,n] = size(u);
  v     = [];
  t     = ['x' 'y'];
  v     = '0';

  if m*n==1 then error('scalar not accepted'); end;

  if (n <> 1) then
    for i=1:m
      x = '0';

      for j=1:n
        x = addf(x,derive(t(j),u(i,j)));
      end

      v(i) = x;
    end
  else
    v = '0';

    for i=1:m
      v = addf(v,derive(t(i),u(i)));
    end
  end

  v = tlist(['ffeq','value'],v)
endfunction
