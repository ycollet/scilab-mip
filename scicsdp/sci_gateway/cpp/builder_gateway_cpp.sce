// ====================================================================
// Copyright Yann COLLETTE 2010-2012
// This file is released into the public domain
// ====================================================================

path_builder = get_absolute_file_path('builder_gateway_cpp.sce');

list_add_inter   = ['csdp',            'scicsdp'; ...
                    'sdpa_read_prob',  'sci_sdpa_read_prob'; ...
                    'sdpa_write_prob', 'sci_sdpa_write_prob'];

files_to_compile = ['scicsdp.cpp', 'sdpa_read_prob.cpp', 'sdpa_write_prob.cpp', 'overload_printf.c'];

if getos()~='Windows' then
  ldflags = '-Wl,--wrap,printf';
else
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
  ldflags = 'parameters.lib';
end

cflags = '-g -Wall -I. -D__USE_DEPRECATED_STACK_FUNCTIONS__ -DNOSHORTS -DBIT64';
cflags = cflags + ' -I' + path_builder;
cflags = cflags + ' -I' + path_builder + '/../../src/c';

libs = ['../../src/c/libsdp'];

tbx_build_gateway('scicsdp_cpp', list_add_inter, files_to_compile, path_builder, libs, ldflags, cflags);

clear tbx_build_gateway;
