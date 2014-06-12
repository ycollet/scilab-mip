// ====================================================================
// Copyright Yann COLLETTE 2011-2012
// This file is released into the public domain
// ====================================================================

list_add_inter      = [];
files_to_compile    = [];
include_clp         = [];
include_cbc         = [];
uuid_lib            = [];
coins_clp_lib       = [];
coins_cbc_lib       = [];
windows_sdk_include = [];
cflags              = [];
ldflags             = [];

path_builder = get_absolute_file_path('builder_gateway_cpp.sce');

list_add_inter   = ['scisymphony', 'scisymphony'];
files_to_compile = ['scisymphony.cpp', 'overload_printf.cpp'];


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

  base_dir        = '../../thirdparty/win/';
  //base_dir        = 'e:/Scilab/toolboxes/';
  symphony_base_dir = 'SYMPHONY-5.4.1/';
  symphony_arch_dir = base_dir + symphony_base_dir;
  cbc_dir           = symphony_arch_dir + 'Cbc/';
  clp_dir           = symphony_arch_dir + 'Clp/';
  symphony_dir      = symphony_arch_dir + 'SYMPHONY/';

  /////////////////////////////
  // Definition of libraries //
  /////////////////////////////

  symphony_lib = '';
  symphony_lib = symphony_lib + ' ' + symphony_dir + 'MSVisualStudio/v10/Release/libOsiCommonTest.lib';
  symphony_lib = symphony_lib + ' ' + symphony_dir + 'MSVisualStudio/v10/Win32/Release/libCgl.lib';
  symphony_lib = symphony_lib + ' ' + symphony_dir + 'MSVisualStudio/v10/Win32/Release/libClp.lib';
  symphony_lib = symphony_lib + ' ' + symphony_dir + 'MSVisualStudio/v10/Win32/Release/libCoinUtils.lib';
  symphony_lib = symphony_lib + ' ' + symphony_dir + 'MSVisualStudio/v10/Win32/Release/libOsi.lib';
  symphony_lib = symphony_lib + ' ' + symphony_dir + 'MSVisualStudio/v10/Win32/Release/libOsiClp.lib';
  symphony_lib = symphony_lib + ' ' + symphony_dir + 'MSVisualStudio/v10/Win32/Release/libOsiSym.lib';
  symphony_lib = symphony_lib + ' ' + symphony_dir + 'MSVisualStudio/v10/Win32/Release/libSymphony.lib';
  symphony_lib = symphony_lib + ' ' + 'parameters.lib ';

  /////////////////////////////
  // Definitions of includes //
  /////////////////////////////

  symphony_inc = ' -I ' + symphony_dir + 'Applications/MPP/include/';
  symphony_inc = symphony_inc + ' -I ' + symphony_dir + 'Applications/SPP/include/';
  symphony_inc = symphony_inc + ' -I ' + symphony_dir + 'Applications/USER/include/';
  symphony_inc = symphony_inc + ' -I ' + symphony_dir + 'Applications/SPP+CUTS/include/';
  symphony_inc = symphony_inc + ' -I ' + symphony_dir + 'Applications/VRP/include/';
  symphony_inc = symphony_inc + ' -I ' + symphony_dir + 'Applications/VRP/include/heurs/';
  symphony_inc = symphony_inc + ' -I ' + symphony_dir + 'Applications/VRP/include/';
  symphony_inc = symphony_inc + ' -I ' + symphony_dir + 'Applications/VRP/include/min_cut/';
  symphony_inc = symphony_inc + ' -I ' + symphony_dir + 'Applications/VRP/include/decomp/';
  symphony_inc = symphony_inc + ' -I ' + symphony_dir + 'Applications/MATCH/include/';
  symphony_inc = symphony_inc + ' -I ' + symphony_dir + 'Applications/CNRP/include/';
  symphony_inc = symphony_inc + ' -I ' + symphony_dir + 'include/';
  symphony_inc = symphony_inc + ' -I ' + symphony_arch_dir + 'CoinUtils/src/';
  symphony_inc = symphony_inc + ' -I ' + path_builder;

  cbc_inc = ' -I ' + cbc_dir + 'inc/';
  cbc_inc = cbc_inc + ' -I ' + cbc_dir + 'src/';
  
  clp_inc = ' -I ' + clp_dir + 'inc/';
  clp_inc = clp_inc + ' -I ' + clp_dir + 'src/';
  
  ////////////////////////////////////////
  // Definition of compilator variables //
  ////////////////////////////////////////

  cflags = '-D__USE_DEPRECATED_STACK_FUNCTIONS__';  
  cflags = cflags + ' -I ' + path_builder + ' -I ' + SCI + '/modules/parameters/includes' + ' ';
  cflags = cflags + symphony_inc + ' ';
  cflags = cflags + cbc_inc      + ' ';
  cflags = cflags + clp_inc      + ' ';

  ldflags = '';
  ldflags = ldflags + symphony_lib + ' ';

  libs = [];
else
  // Before exec builder.sce, don't forget to set the LD_LIBRARY_PATH variable if necessary:
  // export LD_LIBRARY_PATH=/local/stow/Cbc-2.2.2/lib/:/local/stow/bonmin-1.0.1/lib:/local/stow/glpk-4.33/lib:$LD_LIBRARY_PATH

  /////////////////////////////////////////////////
  // Adapt the paths wrt your installation paths //
  /////////////////////////////////////////////////

  symphony_dir     = path_builder + '../../thirdparty/symphony/';

  include_symphony = symphony_dir + 'include/coin';
  // -lgomp -lpthread
  symphony_lib     = '-L' + symphony_dir + 'lib -lSym -lOsiClp -lClp -lOsiSym -lOsi -lCgl -lCoinUtils -lbz2';

  cflags  = '-ggdb -D__USE_DEPRECATED_STACK_FUNCTIONS__ -Wall -I.';
  cflags  = cflags + ' -I' + path_builder;
  cflags  = cflags + ' -I' + include_symphony;
  cflags  = cflags + ' -I' + SCI + '/modules/parameters/includes';

  ldflags = symphony_lib + ' -Wl,--wrap,printf';

  libs = [];
end

tbx_build_gateway('sci_symphony', list_add_inter, files_to_compile, path_builder, libs, ldflags, cflags);

clear tbx_build_gateway;
