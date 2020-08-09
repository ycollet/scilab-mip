function v = %ffeq_m_s(a,b)
  a       = a('value');
  [ma,na] = size(a);
  
  for i=1:ma
    for j=1:na
      // val='('+a(i,j)+')'
      v(i,j) = mulf(a(i,j),string(b));
    end
  end

  v = tlist(['ffeq','value'],v);
endfunction
