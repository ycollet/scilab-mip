function path_out = ampl_path()
  path_out = get_function_path('ampl_path');
  path_out = strsubst(path_out,'macros/ampl_path.sci','');
endfunction
