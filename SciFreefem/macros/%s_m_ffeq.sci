function v=%s_m_ffeq(a,b)
  b       = b('value');
  [mb,nb] = size(b);
  
  for i=1:mb
    for j=1:nb
      // val = '('+b(i,j)+')';
      v(i,j) = mulf(string(a),b(i,j));
    end
  end

  v = tlist(['ffeq','value'],v);
endfunction

