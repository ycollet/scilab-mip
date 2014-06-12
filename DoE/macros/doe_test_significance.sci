function result = doe_test_significance(param_mean, param_var, size_stat, val_to_comp, level, operation)
// Statistical test based on the student test.
// This test allows to test the level of significance of a parameter
if (~isdef('operation','local')) then
  operation = '==';
end
if (~isdef('val_to_comp','local')) then
  val_to_comp = 0.0;
end
if (~isdef('level','local')) then
  level = 0.05;
end
T = (param_mean - val_to_comp)/param_var*sqrt(size_stat-1);
Threshold = cdft('T',size_stat-1, 1 - level/2, level/2);

select operation
case '==' then
  result = (T<Threshold)&(T>-Threshold);
case '<=' then
  result = (T<Threshold);
case '>=' then
  result = (T>-Threshold);
else
  error('doe_test_significance: wrong operator');
end
endfunction
