function v = %ffeq_a_ffeq(a,b)
  a       = a('value')
  b       = b('value')
  [ma,na] = size(a)
  [mb,nb] = size(b)

  if na<>nb then error('incompatible factors'); end;
  if ma<>mb then error('incompatible factors'); end;

  for i=1:ma
    for j=1:na
      v(i,j) = addf(a(i,j),b(i,j));
    end
  end

  v = tlist(['ffeq','value'],v)
endfunction

