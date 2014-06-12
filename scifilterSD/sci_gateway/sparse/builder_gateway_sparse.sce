// ====================================================================
// Yann COLLETTE
// Copyright 2012
// This file is released into the public domain
// ====================================================================

filtersd_path = get_absolute_file_path('builder_gateway_sparse.sce');

// For the sparse version:
sparse_list_add_inter = ['filtersd_sparse', 'sci_filtersd_sparse'; ...
		         'glcpd_sparse',    'sci_glcpd_sparse'];

sparse_files_to_compile = ['sci_filtersd_sparse.c', 'sci_glcpd_sparse.c', ...
			   'func_utils.c', 'sp_sparseA.f', 'sp_filterSD.f', ...
			   'sp_glcpd.f', 'sci_sp_glcpd.f', 'sp_l1sold.f', 'sp_schurQR.f', 'sp_util.f'];

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

tbx_build_gateway('filtersd_sparse_c', sparse_list_add_inter, sparse_files_to_compile, filtersd_path, libs, ldflags, cflags, fflags);

clear tbx_build_gateway;
