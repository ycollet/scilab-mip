function result = doe_test_mean(x, y, level, operation)
// This statistical test is based on the normal test
// It allows to test for equality of mean of two samples

if (size(x,1)==1) then x = x'; end
if (size(y,1)==1) then y = y'; end
if (~isdef('operation','local')) then
  operation = '==';
end
if (~isdef('level','local')) then
  level = 0.05;
end
if ((level<0)|(level>1)) then
  error('doe_test_mean: level must be comprised between 0 and 1');
end
nx = size(x,1);
ny = size(y,1);
if (nx>30) then
  // A statistical test for big population (based on the normal distribution)
  sigmaxy = sqrt(stdev(x)^2/nx + stdev(y)^2/ny);
  Z = (mean(x) - mean(y)) / sigmaxy;
else
  // A statistical test for small population (based on the student distribution)
  sigma = sqrt((nx*stdev(x)^2+ny*stdev(y)^2)/(nx + ny - 2));
  Z = (mean(x) - mean(y))/(sigma*sqrt(1/nx + 1/ny));
end

Threshold = cdfnor('X',0,1, 1-level/2, level/2);

select operation
case '==' then
  result = (Z<Threshold)&(Z>-Threshold);
case '<=' then
  result = (Z<Threshold);
case '>=' then
  result = (Z>-Threshold);
else
  error('doe_test_mean: wrong operator');
end
endfunction
