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

list_add_inter = ['sciglpk', 'sciglpk'];

//if (getos()=='Windows') then
//  unix_w('copy /s /q thirdparty\win\glpk-' + VERSION + '\examples\oldapi\lpx.* sci_gateway\cpp\')
//  files_to_compile  = ['sciglpk.cpp', 'lpx.c'];
//else
//  files_to_compile  = ['sciglpk.cpp'];
//end
files_to_compile  = ['sciglpk.cpp', 'lpx.c'];

if getos()=='Windows' then
  // rebuild parameters.lib
  // We need to use Visual studio 10.0
  if win64() then
    machine = 'X64';
  else
    machine = 'X86';
  end
  //status = rebuild_lib_windows(path_builder,'parameters',machine,'12.0');
  //if ~status then
  //  printf('Error: problem while rebuilding parameters.lib\n');
  //  abort();
  //end
  ////////////////////////////////////////////////////
  // Definition of base paths                       //
  // This part need to be modified if you work with //
  // other version of tools                         //
  ////////////////////////////////////////////////////

  /////////////////////////////
  // Definition of libraries //
  /////////////////////////////
  
  base_dir = path_builder + '\..\..\thirdparty\win\';
  
  tmp_version = strsubst(VERSION,'.','_');
  if win64() then
    glpk_lib = base_dir + 'glpk-' + VERSION +'\w64\glpk_' + tmp_version + '.lib';
  else
    glpk_lib = base_dir + 'glpk-' + VERSION +'\w32\glpk_' + tmp_version + '.lib';
  end

  glpk_lib = glpk_lib + ' ' + SCI + '\bin\parameters.lib ';

  /////////////////////////////
  // Definitions of includes //
  /////////////////////////////

  glpk_inc = ' -I ' + base_dir + 'glpk-' + VERSION + '/src';
  
  ////////////////////////////////////////
  // Definition of compilator variables //
  ////////////////////////////////////////

  cflags = ' -D__USE_DEPRECATED_STACK_FUNCTIONS__';  
  cflags = cflags + ' -I ' + path_builder + ' -I ' + SCI + '\modules\parameters\includes' + ' ';
  cflags = cflags + glpk_inc  + ' ';

  ldflags = '';
  ldflags = ldflags + glpk_lib  + ' ';

  libs = [];
  
  clear tmp_version
else
  /////////////////////////////////////////////////
  // Adapt the paths wrt your installation paths //
  /////////////////////////////////////////////////

  base_dir = path_builder + '/../../thirdparty/glpk/';

  include_glpk = ' -I ' + base_dir + '/include';
  
  glpk_lib = base_dir + '/lib/libglpk.a';

  cflags = '-ggdb -Wall -I. -D__USE_DEPRECATED_STACK_FUNCTIONS__';
  cflags = cflags + ' -I' + path_builder;
  cflags = cflags + include_glpk;
  cflags = cflags + ' -I' + SCI + '/modules/parameters/includes';

  ldflags = glpk_lib;

  libs = [];
end

disp(files_to_compile)

tbx_build_gateway('sci_glpk', list_add_inter, files_to_compile, path_builder, libs, ldflags, cflags);

clear tbx_build_gateway;
