function result = doe_test_var(x, y, level, operation)
// This statistical test is based on the Fisher test
// this test allows to test the variance of two samples
if (size(x,1)==1) then x = x'; end
if (size(y,1)==1) then y = y'; end
if (~isdef('operation','local')) then
  operation = '==';
end
if (~isdef('level','local')) then
  level = 0.05;
end
if ((level<0)|(level>1)) then
  error('doe_test_var: level must be comprised between 0 and 1');
end
nx = size(x,1);
ny = size(y,1);
F = (nx/(nx-1)*stdev(x)^2)/(ny/(ny-1)*stdev(y)^2);
Threshold = cdff('F', nx - 1, ny - 1, 1 - level/2, level/2);
select operation
case '==' then
  result = (F<Threshold)&(F>-Threshold);
case '<=' then
  result = (F<Threshold);
case '>=' then
  result = (F>-Threshold);
else
  error('doe_test_var: wrong operator');
end
endfunction
