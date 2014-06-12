function v=%ffeq_r_s(a,b)
  a       = a('value');
  [ma,na] = size(a);
  
  for i=1:ma
    for j=1:na
      v(i,j) = rdivf(a(i,j),string(b));
    end
  end
  v = tlist(['ffeq','value'],v);
endfunction
