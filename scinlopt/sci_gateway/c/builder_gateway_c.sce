// ====================================================================
// Yann COLLETTE
// Copyright 2010-2012
// This file is released into the public domain
// ====================================================================

path = get_absolute_file_path('builder_gateway_c.sce');

nlopt_install_path = path + '../../thirdparty/';

if getos()~='Windows' then
  nlopt_install_path = nlopt_install_path + 'nlopt/';
else
  nlopt_install_path = nlopt_install_path + 'win/';
end

TABLE = ['nlopt',        'sci_nlopt'; ...
         'nlopt_version','sci_nlopt_version'];
FILES = ['sci_nlopt.c','sci_nlopt_version.c'];

if getos()~='Windows' then
  CFLAGS = '-ggdb -D__USE_DEPRECATED_STACK_FUNCTIONS__ -I' + nlopt_install_path + 'include -I' + path;
  LDFLAGS = '-L' + nlopt_install_path + 'lib -lnlopt_cxx'
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

  CFLAGS = '-D__USE_DEPRECATED_STACK_FUNCTIONS__ -I' + nlopt_install_path + ' -I' + path;
  LDFLAGS = nlopt_install_path + 'libnlopt-0.lib parameters.lib';
end

tbx_build_gateway('nlopt_c', TABLE, FILES, path, [], LDFLAGS, CFLAGS);

clear tbx_build_gateway;
