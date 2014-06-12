//////////////////////////////////////////
// Definition of the objective function //
//////////////////////////////////////////

function [f]=calfun(n, x)
y = ones(10000,1);

for j=1:n
  y(j*10-9) = 1;
  y(j*10-8) = x(j)*2-1;
end

for i=2:n
  for j=1:n
    y(i + 1 + j * 10 - 11+1) = y(j * 10 - 9+1) * 2. * y(i+ j * 10 - 11+1) - y(i - 1 + j * 10 - 11+1);
  end
end

f = 0;
np = n+1;
iw = 1;

for i=1:np
  s = 0.;
  for j=1:n
    s = s+y(i + j * 10 - 11+1);
  end
  s = s/n;
  if iw > 0 then
    s = s + 1/ (i* i - (i*2));
  end
  iw = -iw;
  f = f + s*s;
end
endfunction

//////////////////
// Main program //
//////////////////

iprint = 2;
maxfun = 5000;
rhoend = 1e-6;
x = ones(10,1);

for n=2:2:8
  npt = 2*n + 1;
  x = ones(n,1);
  
  for i=1:n
    x(i) = i/(n+1);  
  end

  rhobeg = x(1) * 0.2;
  disp("Results with N = "+string(n)+ "and NPT = "+string(npt));
  
  xopt = newuoa(npt, x, rhobeg, rhoend, iprint, maxfun, calfun);
end
	
