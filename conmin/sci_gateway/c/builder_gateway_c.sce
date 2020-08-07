// ====================================================================
// Yann COLLETTE
// Copyright 2009-2012
// This file is released into the public domain
// ====================================================================

conmin_path = get_absolute_file_path('builder_gateway_c.sce');

list_add_inter = ['conmin_internal','sci_conmin'];

if getos()=='Windows' then
  files_to_compile = ['conmin_int.c','common_conmin.c'];

  // rebuild parameters.lib
  exec(path_builder + 'rebuild_lib_windows.sci');
  // We need to use Visual studio 10.0
  if win64() then
    machine = 'X64';
  else
    machine = 'X86';
  end
  status = rebuild_lib_windows(path_builder,'parameters',machine,'10.0');
  if ~status then
    printf('Error: problem while rebuilding parameters.lib\n');
    abort();
  end

  cflags  = '-D__USE_DEPRECATED_STACK_FUNCTIONS__ -I' + conmin_path;
  ldflags =  'parameters.lib ';
  libs    = ['../../src/fortran/libget_ct'];
else
  files_to_compile = ['conmin_int.c'];
  cflags  = '-D__USE_DEPRECATED_STACK_FUNCTIONS__ -I' + conmin_path;
  ldflags = '';
  libs    = ['../../src/fortran/libget_ct'];
end

tbx_build_gateway('conmin_c', list_add_inter, files_to_compile, get_absolute_file_path('builder_gateway_c.sce'), libs, ldflags, cflags);

clear tbx_build_gateway;
