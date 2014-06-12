function r2 = gammatest_r2(data_in,a)
  r2 = 1 - abs(a) / stdev(data_in(:,$));
endfunction
