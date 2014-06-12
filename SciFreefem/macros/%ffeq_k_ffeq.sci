function v=%ffeq_k_ffeq(a,b)
  a       = a('value')
  b       = b('value')
  [ma,na] = size(a)
  [mb,nb] = size(b)
  
  if (na<>nb) then error('incompatible factors'); end;
  if (ma<>mb) then error('incompatible factors'); end;

  v = '0';

  if (na>1) then
    for i=1:ma
      for j=1:na
        v = addf(v,mulf(a(i,j),b(i,j)));
      end
    end
  else
    if (ma>2) then error('sorry: not implemented'); end;

    v = subf(mulf(a(1),b(2)),mulf(a(2),b(1)))
  end

  v = tlist(['ffeq','value'],v)
endfunction

