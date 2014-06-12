function w=%ffeq_s(a)
  a     = a('value');
  [m,n] = size(a);

  for i=1:m
    for j=1:n
      b = strindex(a(i,j),'-');

      if (b==[]) then
        a(i,j) = '-'+a(i,j);
      else
        a(i,j) = stripblanks(strsubst(a(i,j),'-',' '));
      end
    end
  end

  w = tlist(['ffeq','value'],a);
endfunction
