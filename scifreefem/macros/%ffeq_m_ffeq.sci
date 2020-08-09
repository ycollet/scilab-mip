function v=%ffeq_m_ffeq(a,b)
  a       = a('value');
  b       = b('value');
  [ma,na] = size(a);
  [mb,nb] = size(b);
  
  if na==mb then
    for i=1:ma
      for j=1:nb
        x = '0';

        for k=1:na
          x = addf(x,mulf(a(i,k),b(k,j)));
        end

        v(i,j) = x;
      end
    end
  else
    if (na==nb) then
      if (ma==mb) then
        if (na==1) then
          m = ma;
        else
          if (ma==1) then
            m = na;
          else
            error('incompatible factors');
          end
        end
  
        v = '0';

        for i=1:m
          v = addf(v,mulf(a(i),b(i)));
        end
      else
        error('incompatible factors');
      end
    else
      error('incompatible factors');
    end
  end

  v = tlist(['ffeq','value'],v);
endfunction
