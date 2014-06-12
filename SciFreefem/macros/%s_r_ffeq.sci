function v=%s_r_ffeq(a,b)
  b       = b('value')
  [mb,nb] = size(b)
  
  for i=1:mb
    for j=1:nb
      v(i,j) = rdivff(string(a),b(i,j));
    end
  end

  v = tlist(['ffeq','value'],v);
endfunction
