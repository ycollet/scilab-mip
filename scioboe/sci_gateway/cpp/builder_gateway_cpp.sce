// ====================================================================
// Copyright Yann COLLETTE 2010-2012
// This file is released into the public domain
// ====================================================================

list_add_inter      = [];
files_to_compile    = [];
cflags              = [];
ldflags             = [];
oldlib              = [];

path_builder     = get_absolute_file_path('builder_gateway_cpp.sce');

list_add_inter   = ['oboe', 'scioboe'];

files_to_compile = ['scioboe.o'];

if getos()=='Windows' then
  ////////////////////////////////////////////////////
  // Definition of base paths                       //
  // This part need to be modified if you work with //
  // other version of tools                         //
  ////////////////////////////////////////////////////

  base_dir     = 'z:/toolboxes_google/oboe_thirdparty/';
  oboe_dir     = base_dir + 'OBOE-1.0.2-patched/';
  lapackpp_dir = base_dir + 'lapackpp-2.5.3/';

  /////////////////////////////
  // Definition of libraries //
  /////////////////////////////

  oboe_lib     = oboe_dir + 'WIN32/Release/OBOE.vs5.lib';
  lapackpp_lib = lapackpp_dir + 'Release/lapackpp.lib';

  /////////////////////////////
  // Definitions of includes //
  /////////////////////////////

  oboe_inc = ' -I ' + oboe_dir + 'WIN32/';
  oboe_inc = oboe_inc + ' -I ' + oboe_dir;
  oboe_inc = oboe_inc + ' -I ' + oboe_dir + 'include/';
  oboe_inc = oboe_inc + ' -I ' + oboe_dir + 'src/Utilities/';
  oboe_inc = oboe_inc + ' -I ' + oboe_dir + 'src/AccpmLA/';
  oboe_inc = oboe_inc + ' -I ' + oboe_dir + 'src/UI/';
  oboe_inc = oboe_inc + ' -I ' + oboe_dir + 'src/Benchmark/QP/';
  oboe_inc = oboe_inc + ' -I ' + oboe_dir + 'src/Benchmark/NNPOL/';
  oboe_inc = oboe_inc + ' -I ' + oboe_dir + 'src/Benchmark/MCF/';
  oboe_inc = oboe_inc + ' -I ' + oboe_dir + 'src/Oracle/';
  oboe_inc = oboe_inc + ' -I ' + oboe_dir + 'src/AccpmCore/';
  oboe_inc = oboe_inc + ' -I ' + oboe_dir + 'src/ProblemInput/';
  
  lapackpp_inc = ' -I ' + lapackpp_dir;
  lapackpp_inc = lapackpp_inc + ' -I ' + lapackpp_dir + 'include/';
  lapackpp_inc = lapackpp_inc + ' -I ' + lapackpp_dir + 'matrix/src/';

  // MUMPS need ifconsol.lib
  ifortran_lib = '""C:\Program Files\Intel\Compiler\11.0\061\fortran\lib\ia32\ifconsol.lib""';
  // OBOE need also these ones
  ifortran_lib = ifortran_lib + ' ' + '""C:\Program Files\Intel\Compiler\11.0\061\fortran\lib\ia32\libifcoremt.lib""';
  ifortran_lib = ifortran_lib + ' ' + '""C:\Program Files\Intel\Compiler\11.0\061\fortran\lib\ia32\libifport.lib""';
  ifortran_lib = ifortran_lib + ' ' + '""C:\Program Files\Intel\Compiler\11.0\061\fortran\lib\ia32\libmmt.lib""';
  ifortran_lib = ifortran_lib + ' ' + '""C:\Program Files\Intel\Compiler\11.0\061\fortran\lib\ia32\libirc.lib""';
  ifortran_lib = ifortran_lib + ' ' + '""C:\Program Files\Intel\Compiler\11.0\061\fortran\lib\ia32\svml_disp.lib""';

  ////////////////////////////////////////
  // Definition of compilator variables //
  ////////////////////////////////////////

  cflags  = '-D__USE_DEPRECATED_STACK_FUNCTIONS__';  
  cflags  = cflags + ' -I ' + path_builder + ' -I ' + SCI + '/modules/parameters/includes' + ' ';
  cflags  = cflags + oboe_inc     + ' ';
  cflags  = cflags + lapackpp_inc + ' ';

  ldflags = '';
  ldflags = ldflags + lapackpp_lib + ' ';
  ldflags = ldflags + oboe_lib     + ' ';
  ldflags = ldflags + ifortran_lib + ' ';

  libs = [];
else
  /////////////////////////////////////////////////
  // Adapt the paths wrt your installation paths //
  /////////////////////////////////////////////////

  tools_path   = '/opt/stow/';
  oboe_dir     = 'oboe-1.0.2/';
  lapackpp_dir = 'lapack++-2.5.3/';

  include_oboe     = tools_path + oboe_dir + 'include';
  include_lapackpp = tools_path + lapackpp_dir + 'include/lapackpp';
  
  oboe_lib     = '-L' + tools_path + oboe_dir + 'lib -laccpm -laccpmcore -laccpmla -laccpmoracle -laccpmparam';
  lapackpp_lib = '-L' + tools_path + lapackpp_dir + 'lib -llapackpp';
 
  cflags = '-ggdb -D__USE_DEPRECATED_STACK_FUNCTIONS__ -Wall -I.';
  cflags = cflags + ' -I' + path_builder;
  cflags = cflags + ' -I' + SCI + '/modules/parameters/includes';
  cflags = cflags + ' -I' + include_oboe;
  cflags = cflags + ' -I' + include_lapackpp;

  ldflags = oboe_lib + ' ' + lapackpp_lib;

  libs = [];
  
  old_ldlib = getenv('LD_LIBRARY_PATH');

  // Set some environment variables
  setenv('LD_LIBRARY_PATH', tools_path + lapackpp_dir + 'lib:' + getenv('LD_LIBRARY_PATH'));
  setenv('LD_LIBRARY_PATH', tools_path + oboe_dir     + 'lib:' + getenv('LD_LIBRARY_PATH'));

  setenv('LDFLAGS',ldflags);
end

tbx_build_gateway('sci_coinor', list_add_inter, files_to_compile, path_builder, libs, ldflags, cflags);

clear tbx_build_gateway;

if getos()~='Windows' then
  setenv('LD_LIBRARY_PATH',old_ldlib);
  setenv('LDFLAGS','');
end
