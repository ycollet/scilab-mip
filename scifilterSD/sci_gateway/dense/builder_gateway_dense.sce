// ====================================================================
// Yann COLLETTE
// Copyright 2012
// This file is released into the public domain
// ====================================================================

filtersd_path = get_absolute_file_path('builder_gateway_dense.sce');

// For the dense version:
dense_list_add_inter = ['filtersd_dense', 'sci_filtersd_dense'; ...
		        'glpcd_dense',    'sci_glcpd_dense'];

dense_files_to_compile = ['sci_filtersd_dense.c', 'sci_glcpd_dense.c', ...
			  'func_utils.c', 'ds_denseA.f', 'ds_denseL.f', ...
                          'ds_filterSD.f', 'sci_ds_glcpd.f', 'ds_glcpd.f', 'ds_l1sold.f', 'ds_util.f'];

if getos()=='Windows' then
  if (findmsifortcompiler()=='unknown') then
    // F2C used
    setenv('F2C_IMPORT_COMMON','YES');
  end

  // rebuild parameters.lib
  exec(path_builder + 'rebuild_lib_windows.sci');
  // We need to use Visual studio 10.0
  if win64() then
    machine = 'X64';
  else
    machine = 'X86';
  end
  status = rebuild_lib_windows(filtersd_path,'parameters',machine,'10.0');
  if ~status then
    printf('Error: problem while rebuilding parameters.lib\n');
    abort();
  end

  cflags  = '-D__USE_DEPRECATED_STACK_FUNCTIONS__ -I' + filtersd_path;
  fflags  = '-D__USE_DEPRECATED_STACK_FUNCTIONS__';
  ldflags = 'parameters.lib';
else
  cflags  = '-g -Wall -D__USE_DEPRECATED_STACK_FUNCTIONS__ -I' + filtersd_path;
  fflags  = '-g -Wall -D__USE_DEPRECATED_STACK_FUNCTIONS__';
  ldflags = '-lgfortran';
end

libs = '';

tbx_build_gateway('filtersd_dense_c',  dense_list_add_inter, dense_files_to_compile, filtersd_path, libs, ldflags, cflags, fflags);

clear tbx_build_gateway;
