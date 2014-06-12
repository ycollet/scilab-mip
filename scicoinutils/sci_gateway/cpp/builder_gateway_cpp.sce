// ====================================================================
// Copyright Yann COLLETTE 2010-2012
// This file is released into the public domain
// ====================================================================

list_add_inter      = [];
files_to_compile    = [];
include_glpk        = [];
glpk_lib            = [];
cflags              = [];
ldflags             = [];

path_builder = get_absolute_file_path('builder_gateway_cpp.sce');

list_add_inter      = ['write_lp_file',    'write_lp_file'; ...
                       'write_mps_file',   'write_mps_file'; ...
                       'read_lp_file',     'read_lp_file'; ...
                       'read_lp_file_mp',  'read_lp_file_mp'; ...
                       'read_mps_file',    'read_mps_file'; ...
                       'read_mps_file_mp', 'read_mps_file_mp'];

files_to_compile    = ['read_mps_file.o', 'read_lp_file.o', ...
                       'write_lp_file.o', 'write_mps_file.o','helper.o'];

if getos()=='Windows' then
  ////////////////////////////////////////////////////
  // Definition of base paths                       //
  // This part need to be modified if you work with //
  // other version of tools                         //
  ////////////////////////////////////////////////////

  base_dir = path_builder + '\..\..\thirdparty\win\';

  /////////////////////////////
  // Definition of libraries //
  /////////////////////////////

  if Build_64Bits then
    coinutils_lib = base_dir + '\CoinUtils-' + VERSION + '\CoinUtils\MSVisualStudio\v10\Win64\Release\libCoinUtils.lib';
  else
    coinutils_lib = base_dir + '\CoinUtils-' + VERSION + '\CoinUtils\MSVisualStudio\v10\Win32\Release\libCoinUtils.lib';
  end

  /////////////////////////////
  // Definitions of includes //
  /////////////////////////////

  coinutils_inc = ' -I ' + base_dir + '/CoinUtils-' + VERSION + '/CoinUtils/src';
  
  ////////////////////////////////////////
  // Definition of compilator variables //
  ////////////////////////////////////////

  cflags = '-D__USE_DEPRECATED_STACK_FUNCTIONS__';  
  cflags = cflags + ' -I ' + path_builder + ' -I ' + SCI + '/modules/parameters/includes' + ' ';
  cflags = cflags + coinutils_inc  + ' ';

  ldflags = '';
  ldflags = ldflags + coinutils_lib  + ' ';

  libs = [];
else
  /////////////////////////////////////////////////
  // Adapt the paths wrt your installation paths //
  /////////////////////////////////////////////////

  base_dir = path_builder + '/../../thirdparty/coinutils/';

  include_coinutils = ' -I ' + base_dir + '/include/coin';
  
  coinutils_lib = base_dir + '/lib/libCoinUtils.a';

  cflags = '-ggdb -Wall -D__USE_DEPRECATED_STACK_FUNCTIONS__ -I.';
  cflags = cflags + ' -I' + path_builder;
  cflags = cflags + include_coinutils;
  cflags = cflags + ' -I' + SCI + '/modules/parameters/includes';

  ldflags = coinutils_lib;

  libs = [];
end

tbx_build_gateway('scicoinutils', list_add_inter, files_to_compile, get_absolute_file_path('builder_gateway_cpp.sce'), libs, ldflags, cflags);

clear tbx_build_gateway;
