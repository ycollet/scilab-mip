// ====================================================================
// Copyright Yann COLLETTE 2010-2012
// This file is released into the public domain
// ====================================================================

list_add_inter      = [];
files_to_compile    = [];
include_lpsolve     = [];
lpsolve_lib         = [];
cflags              = [];
ldflags             = [];

path_builder     = get_absolute_file_path('builder_gateway_cpp.sce');

list_add_inter   = ['scilpsolve', 'scilpsolve'];

files_to_compile = ['scilpsolve.cpp'];

if getos()=='Windows' then
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

  ////////////////////////////////////////////////////
  // Definition of base paths                       //
  // This part need to be modified if you work with //
  // other version of tools                         //
  ////////////////////////////////////////////////////

  base_dir = path_builder + '..\..\thirdparty\win\';

  /////////////////////////////
  // Definition of libraries //
  /////////////////////////////

  lpsolve_lib = base_dir + 'liblpsolve55.lib';
  lpsolve_lib = lpsolve_lib + ' ' + base_dir + 'lpsolve55.lib';
  lpsolve_lib = lpsolve_lib + ' ' + 'parameters.lib ';

  /////////////////////////////
  // Definitions of includes //
  /////////////////////////////

  lpsolve_inc = ' -I ' + base_dir;
  
  ////////////////////////////////////////
  // Definition of compilator variables //
  ////////////////////////////////////////

  cflags = '-D__USE_DEPRECATED_STACK_FUNCTIONS__';  
  cflags = cflags + ' -I ' + path_builder + ' -I ' + SCI + '/modules/parameters/includes' + ' ';
  cflags = cflags + lpsolve_inc  + ' ';

  ldflags = '';
  ldflags = ldflags + lpsolve_lib  + ' ';

  libs = [];
else
  // Add CXXFLAGS=-DCOIN_DO_PDCO CFLAGS=-DCOIN_DO_PDCO in the configure parameters when compiling Bonmin
  // 
  // Before exec builder.sce, don't forget to set the LD_LIBRARY_PATH variable if necessary:
  // export LD_LIBRARY_PATH=/local/stow/Cbc-2.2.2/lib/:/local/stow/bonmin-1.0.1/lib:/local/stow/glpk-4.33/lib:$LD_LIBRARY_PATH

  /////////////////////////////////////////////////
  // Adapt the paths wrt your installation paths //
  /////////////////////////////////////////////////

  base_dir = path_builder + '/../../thirdparty/lp_solve_5.5/';

  include_lpsolve = ' -I ' + base_dir;
  include_lpsolve = include_lpsolve + ' -I ' + base_dir + '/bfp';
  include_lpsolve = include_lpsolve + ' -I ' + base_dir + '/colamd';
  include_lpsolve = include_lpsolve + ' -I ' + base_dir + '/lp_solve';
  include_lpsolve = include_lpsolve + ' -I ' + base_dir + '/lpsolve55';
  include_lpsolve = include_lpsolve + ' -I ' + base_dir + '/shared';
  
  lpsolve_lib = base_dir + '/lpsolve55/bin/ux*/liblpsolve55.a';

  cflags = '-ggdb -Wall -D__USE_DEPRECATED_STACK_FUNCTIONS__ -I.';
  cflags = cflags + ' -I' + path_builder;
  cflags = cflags + include_lpsolve;
  cflags = cflags + ' -I' + SCI + '/modules/parameters/includes';

  ldflags = lpsolve_lib;

  libs = [];
end

tbx_build_gateway('sci_lpsolve', list_add_inter, files_to_compile, path_builder, libs, ldflags, cflags);

clear tbx_build_gateway;
