// ====================================================================
// Copyright Yann COLLETTE 2011
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

old_ldlib  = [];

use_thread = %F;
use_mingw  = %F;

list_add_inter      = ['sciclp',      'sciclp'; ...
                       'sciosi',      'sciosi'; ...
                       'scicbc',      'scicbc'; ...
                       'scibonmin',   'scibonmin'; ...
                       'sciipopt',    'sciipopt';];

files_to_compile    = ['sciclp.cpp','sciosi.cpp','scicbc.cpp', ...
                       'scibonmin.cpp', 'sciipopt.cpp', ...
                       'scilabjournal.cpp', 'scilabexception.cpp', 'helper.cpp', ...
                       'call_function.cpp', 'manage_ipopt_params.cpp', 'manage_bonmin_params.cpp'];

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

  if ~use_mingw then
    //base_dir        = '../../../scicoinor_thirdparty/';
    base_dir        = 'e:/Scilab/toolboxes/';
    
    bonmin_base_dir = 'Bonmin-1.5.1/';
    bonmin_arch_dir = base_dir + bonmin_base_dir;
    cbc_dir         = bonmin_arch_dir + 'Cbc/';
    bonmin_dir      = bonmin_arch_dir + 'Bonmin/';
    clp_dir         = bonmin_arch_dir + 'Clp/';
    ipopt_dir       = bonmin_arch_dir + 'Ipopt/';
  
    /////////////////////////////
    // Definition of libraries //
    /////////////////////////////
  
    cbc_lib = bonmin_arch_dir + 'Cbc/MSVisualStudio/v10/win32/Debug/libCoinUtils.lib';
    
    cbc_lib = cbc_lib + ' ' + bonmin_arch_dir + 'Ipopt/MSVisualStudio/v8-ifort/CoinBlas/Debug/CoinBlas.lib';
    cbc_lib = cbc_lib + ' ' + bonmin_arch_dir + 'Ipopt/MSVisualStudio/v8-ifort/CoinLapack/Debug/CoinLapack.lib';
    cbc_lib = cbc_lib + ' ' + bonmin_arch_dir + 'Clp/MSVisualStudio/v10/Win32/Debug/libClp.lib';
    cbc_lib = cbc_lib + ' ' + bonmin_arch_dir + 'Cgl/MSVisualStudio/v10/Win32/Debug/libCgl.lib';
    cbc_lib = cbc_lib + ' ' + bonmin_arch_dir + 'Cbc/MSVisualStudio/v10/Win32/Debug/libCbc.lib';
    cbc_lib = cbc_lib + ' ' + bonmin_arch_dir + 'Cbc/MSVisualStudio/v10/Win32/Debug/libCbcSolver.lib';
    cbc_lib = cbc_lib + ' ' + bonmin_arch_dir + 'Cbc/MSVisualStudio/v10/Win32/Debug/libOsi.lib';
    cbc_lib = cbc_lib + ' ' + bonmin_arch_dir + 'Cbc/MSVisualStudio/v10/Win32/Debug/libOsiClp.lib';
    cbc_lib = cbc_lib + ' ' + bonmin_arch_dir + 'Cbc/MSVisualStudio/v10/Win32/Debug/libOsiCbc.lib';
    cbc_lib = cbc_lib + ' ' + bonmin_arch_dir + 'Ipopt/MSVisualStudio/v8-ifort/CoinMumps/Debug/CoinMumpsF90.lib';
    cbc_lib = cbc_lib + ' ' + bonmin_arch_dir + 'Ipopt/MSVisualStudio/v8-ifort/Debug/CoinMumpsC.lib';
    cbc_lib = cbc_lib + ' ' + bonmin_arch_dir + 'Ipopt/MSVisualStudio/v8-ifort/Debug/CoinMetis.lib';
    cbc_lib = cbc_lib + ' ' + bonmin_arch_dir + 'Ipopt/MSVisualStudio/v8/libIpopt/Debug/libIpopt.lib';
    cbc_lib = cbc_lib + ' ' + bonmin_arch_dir + 'Bonmin/MSVisualStudio/v8/libBonmin/Debug/libBonmin.lib';

    scilab_lib = ' parameters.lib ';

    // Fortran compiler path  
    fortran_comp_path = '""C:\Program Files (x86)\Intel\ComposerXE-2011\compiler\lib\ia32""';

    // MUMPS need ifconsol.lib
    ifortran_lib    = fortran_comp_path + '\ifconsol.lib';

    // OBOE need also these ones
    ifortran_lib    = ifortran_lib + ' ' + fortran_comp_path + '\libifcoremt.lib';
    ifortran_lib    = ifortran_lib + ' ' + fortran_comp_path + '\libifport.lib';
    ifortran_lib    = ifortran_lib + ' ' + fortran_comp_path + '\libmmt.lib';
    ifortran_lib    = ifortran_lib + ' ' + fortran_comp_path + '\libirc.lib';
    ifortran_lib    = ifortran_lib + ' ' + fortran_comp_path + '\svml_disp.lib';
  
    /////////////////////////////
    // Definitions of includes //
    /////////////////////////////

    cbc_inc = ' -I ' + cbc_dir + 'inc/';
    cbc_inc = cbc_inc + ' -I ' + cbc_dir + 'src/';
    cbc_inc = cbc_inc + ' -I ' + cbc_dir + 'src/OsiCbc';
    
    ipopt_inc = ' -I ' + ipopt_dir + 'inc/';
    ipopt_inc = ipopt_inc + ' -I ' + ipopt_dir + 'src/Apps/CompositeInterface/';
    ipopt_inc = ipopt_inc + ' -I ' + ipopt_dir + 'src/Apps/AmplSolver/';
    ipopt_inc = ipopt_inc + ' -I ' + ipopt_dir + 'src/Interfaces/';
    ipopt_inc = ipopt_inc + ' -I ' + ipopt_dir + 'src/Common/';
    ipopt_inc = ipopt_inc + ' -I ' + ipopt_dir + 'src/contrib/CGPenalty/';
    ipopt_inc = ipopt_inc + ' -I ' + ipopt_dir + 'src/contrib/LinearSolverLoader/';
    ipopt_inc = ipopt_inc + ' -I ' + ipopt_dir + 'src/LinAlg/TMatrices/';
    ipopt_inc = ipopt_inc + ' -I ' + ipopt_dir + 'src/LinAlg/';
    ipopt_inc = ipopt_inc + ' -I ' + ipopt_dir + 'src/Algorithm/';
    ipopt_inc = ipopt_inc + ' -I ' + ipopt_dir + 'src/Algorithm/Inexact/';
    ipopt_inc = ipopt_inc + ' -I ' + ipopt_dir + 'src/Algorithm/LinearSolvers/';
    ipopt_inc = ipopt_inc + ' -I ' + ipopt_dir + 'MSVisualStudio/v8-ifort/IpOpt/';
    ipopt_inc = ipopt_inc + ' -I ' + ipopt_dir + 'MSVisualStudio/v8-ifort/include/coin/';
    ipopt_inc = ipopt_inc + ' -I ' + ipopt_dir + 'MSVisualStudio/v8-ifort/CoinMetis/';
    ipopt_inc = ipopt_inc + ' -I ' + bonmin_arch_dir + 'ThirdParty/Metis/metis-4.0/Lib/';
    ipopt_inc = ipopt_inc + ' -I ' + bonmin_arch_dir + 'ThirdParty/Mumps/MUMPS/PORD/include/';
    ipopt_inc = ipopt_inc + ' -I ' + bonmin_arch_dir + 'ThirdParty/Mumps/MUMPS/include/';
    ipopt_inc = ipopt_inc + ' -I ' + bonmin_arch_dir + 'ThirdParty/Mumps/MUMPS/src/';
    ipopt_inc = ipopt_inc + ' -I ' + bonmin_arch_dir + 'ThirdParty/Mumps/MUMPS/libseq/';
    
    clp_inc = ' -I ' + clp_dir + 'inc/';
    clp_inc = clp_inc + ' -I ' + clp_dir + 'src/';
    clp_inc = clp_inc + ' -I ' + clp_dir + 'src/OsiClp';
    
    bonmin_inc = ' -I ' + bonmin_dir + 'Osi/inc/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Osi/src/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Osi/src/Osi';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Osi/src/OsiGrb/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Osi/src/OsiCpx/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Osi/src/OsiXpr/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Osi/src/OsiOsl/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Osi/src/OsiFmp/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Osi/src/OsiDylp/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Osi/src/OsiSym/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Osi/src/OsiVol/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Osi/src/OsiSpx/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Osi/src/OsiCbc/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Osi/src/OsiMsk/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Osi/src/OsiClp/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Osi/src/OsiGlpk/';
    
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Cgl/inc/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Cgl/src/CglMixedIntegerRounding2/'
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Cgl/src/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Cgl/src/CglFlowCover/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Cgl/src/CglAllDifferent/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Cgl/src/CglGomory/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Cgl/src/CglOddHole/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Cgl/src/CglSimpleRounding/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Cgl/src/CglPreProcess/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Cgl/src/CglClique/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Cgl/src/CglTwomir/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Cgl/src/CglLandP/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Cgl/src/CglProbing/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Cgl/src/CglMixedIntegerRounding/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Cgl/src/CglLiftAndProject/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Cgl/src/CglKnapsackCover/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Cgl/src/CglDuplicateRow/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Cgl/src/CglRedSplit/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'Cgl/src/CglResidualCapacity/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'BuildTools/headers/';
    
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_dir + 'inc/';	
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_dir + 'src/Apps/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_dir + 'src/CbcBonmin/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_dir + 'src/CbcBonmin/Heuristics/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_dir + 'src/Interfaces/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_dir + 'src/Interfaces/Filter/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_dir + 'src/Interfaces/Ampl/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_dir + 'src/Interfaces/Ipopt/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_dir + 'src/Interfaces/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_dir + 'src/Algorithms/QuadCuts/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_dir + 'src/Algorithms/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_dir + 'src/Algorithms/OaGenerators/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_dir + 'src/Algorithms/Ampl/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_dir + 'src/Algorithms/Branching/';
    
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'CoinUtils/inc/';
    bonmin_inc = bonmin_inc + ' -I ' + bonmin_arch_dir + 'CoinUtils/src/';
    
    ////////////////////////////////////////
    // Definition of compilator variables //
    ////////////////////////////////////////
  
    cflags = '-D__USE_DEPRECATED_STACK_FUNCTIONS__';  
    cflags = cflags + ' -I ' + path_builder + ' -I ' + SCI + '/modules/parameters/includes' + ' ';
    cflags = cflags + ipopt_inc    + ' ';
    cflags = cflags + cbc_inc      + ' ';
    cflags = cflags + clp_inc      + ' ';
    cflags = cflags + bonmin_inc   + ' ';
    cflags = cflags + "-D_ITERATOR_DEBUG_LEVEL=2"
    
    ldflags = '';
    ldflags = ldflags + cbc_lib      + ' ';
    ldflags = ldflags + ifortran_lib + ' ';
    ldflags = ldflags + scilab_lib + ' ';
  
    libs = [];
  else
    base_dir = 'e:/Scilab/toolboxes/Bonmin-compile/';
  
    /////////////////////////////
    // Definition of libraries //
    /////////////////////////////
  
    bonmin_lib = '-L' + base_dir + 'lib ';
    bonmin_lib = bonmin_lib + '-lbonmin ';
    bonmin_lib = bonmin_lib + '-lipopt ';
    bonmin_lib = bonmin_lib + '-lCbcSolver ';
    bonmin_lib = bonmin_lib + '-lCbc ';
    bonmin_lib = bonmin_lib + '-lCgl ';
    bonmin_lib = bonmin_lib + '-lClp ';
    bonmin_lib = bonmin_lib + '-lOsi ';
    bonmin_lib = bonmin_lib + '-lOsiCbc ';
    bonmin_lib = bonmin_lib + '-lOsiClp ';
    bonmin_lib = bonmin_lib + '-lCoinUtils ';
    bonmin_lib = bonmin_lib + '-lOsiCommonTests '; // ??
    bonmin_lib = bonmin_lib + '-lcoinmumps ';
    bonmin_lib = bonmin_lib + '-lcoinlapack ';
    bonmin_lib = bonmin_lib + '-lcoinblas ';
    bonmin_lib = bonmin_lib + '-lcoinmetis ';
    bonmin_lib = bonmin_lib + '-static-libgcc ';
    bonmin_lib = bonmin_lib + '-lstdc++ ';
    bonmin_lib = bonmin_lib + '-lgfortran ';
    //bonmin_lib = bonmin_lib + '-lpthread ';
    bonmin_lib = bonmin_lib + ' parameters.lib ';
  
//    bonmin_dir = '';
//    bonmin_lib = bonmin_lib + base_dir + 'lib/libgcc.lib ';
//    bonmin_lib = bonmin_lib + base_dir + 'lib/libgfortran.lib ';
//    bonmin_lib = bonmin_lib + base_dir + 'lib/libstdc++.lib ';
//    bonmin_lib = bonmin_dir + base_dir + 'lib/libbonmin_yc.lib ';
//    bonmin_lib = bonmin_lib + ' parameters.lib ';

    /////////////////////////////
    // Definitions of includes //
    /////////////////////////////
    
    bonmin_inc = base_dir + 'include/coin/';
    
    ////////////////////////////////////////
    // Definition of compilator variables //
    ////////////////////////////////////////
  
    cflags  = '-D__USE_DEPRECATED_STACK_FUNCTIONS__';
    cflags  = cflags + ' -I' + path_builder + ' -I' + SCI + '/modules/parameters/includes' + ' -I' + bonmin_inc;
  
    ldflags = '';
    ldflags = ldflags + bonmin_lib + ' ';
  
    libs = [];
  end
else
  // Add CXXFLAGS=-DCOIN_DO_PDCO CFLAGS=-DCOIN_DO_PDCO in the configure parameters when compiling Bonmin
  // 
  // Before exec builder.sce, don't forget to set the LD_LIBRARY_PATH variable if necessary:
  // export LD_LIBRARY_PATH=/local/stow/Cbc-2.2.2/lib/:/local/stow/bonmin-1.0.1/lib:/local/stow/glpk-4.33/lib:$LD_LIBRARY_PATH

  /////////////////////////////////////////////////
  // Adapt the paths wrt your installation paths //
  /////////////////////////////////////////////////

  tools_path          = '/opt/';
  bonmin_dir          = 'bonmin-1.5.0/';

  include_bonmin      = tools_path + bonmin_dir + 'include/coin';
  
  bonmin_lib = '-L' + tools_path + bonmin_dir + 'lib -lbonmin -lcoinmumps -lCbc -lCoinUtils -lCbcSolver -lipopt -lCgl';
  bonmin_lib = bonmin_lib + ' -lOsiCbc -lClp -lOsiClp -lcoinblas -lOsiCommonTests -lcoinlapack';
  bonmin_lib = bonmin_lib + ' -lOsi -lcoinmetis -L/usr/lib64 -ldl';
  thread_lib          = '-lpthread';

  cflags              = '-ggdb -Wall -D__USE_DEPRECATED_STACK_FUNCTIONS__ -fpermissive -I.';
  cflags              = cflags + ' -I' + path_builder;
  cflags              = cflags + ' -I' + include_bonmin;
  cflags              = cflags + ' -I' + SCI + '/modules/parameters/includes';
  ldflags             = bonmin_lib + ' ' + thread_lib;

  if use_thread then
    cflags = cflags + ' -DCBC_THREAD=1 ';
  end
  if use_thread then
    ldflags = ldflags + ' ' + thread_lib;
  end

  libs = [];
  
  old_ldlib = getenv('LD_LIBRARY_PATH');

  // Set some environment variables
  setenv('LD_LIBRARY_PATH', tools_path + bonmin_dir + 'lib:' + getenv('LD_LIBRARY_PATH'));

  setenv('LDFLAGS',ldflags);
end

tbx_build_gateway('sci_coinor', list_add_inter, files_to_compile, path_builder, libs, ldflags, cflags);

clear tbx_build_gateway;

if getos()~='Windows' then
  setenv('LD_LIBRARY_PATH',old_ldlib);
  setenv('LDFLAGS','');
end
