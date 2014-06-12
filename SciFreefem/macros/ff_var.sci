function [] = ff_var(name,value)
  str = name + ':=' + string(value) + ';';
  ff_exec(str);
endfunction
