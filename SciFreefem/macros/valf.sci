function [tCoor] = valf(f,t)
  st    = size(t);
  [x,y] = f(t);

  sz    = size(x);
  if (sz(2) == 1) then
    x = x * ones(1,st(2))
  end

  sz = size(y);
  if (sz(2) == 1) then
    y = y * ones(1,st(2))
  end

  tCoor = [x' y'];
endfunction
