function v = pde_varsol(name,integral,a)
  a  = a('value');
  st = '';
  if (integral <> 'domain') then
    st = integral;
  end
  
  v = name + '=int('+st+')('+a+')';
endfunction
