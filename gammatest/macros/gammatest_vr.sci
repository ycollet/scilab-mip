function vr = gammatest_vr(data_in,a)
  vr = abs(a)/stdev(data_in(:,$));
endfunction
