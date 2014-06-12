// ====================================================================
// Yann COLLETTE
// Christopher MAES
// Copyright 2009-2012
// This file is released into the public domain
// ====================================================================

// Trilinos 10.2.0 has been compiled using these flags:
// In the root of the trilinos package:
// export TRILINOS_HOME=`pwd`
// Then:
// cmake -D CMAKE_BUILD_TYPE:STRING=DEBUG \
//       -D Trilinos_ENABLE_amesos:BOOL=ON \
//       -D Trilinos_ENABLE_anasazi:BOOL=ON \
//       -D Trilinos_ENABLE_aztecoo:BOOL=ON \
//       -D Trilinos_ENABLE_belos:BOOL=ON \
//       -D Trilinos_ENABLE_common:BOOL=ON \
//       -D Trilinos_ENABLE_didasko:BOOL=ON \
//       -D Trilinos_ENABLE_epetra:BOOL=ON \
//       -D Trilinos_ENABLE_epetraext:BOOL=ON \
//       -D Trilinos_ENABLE_fei:BOOL=ON \
//       -D Trilinos_ENABLE_galeri:BOOL=ON \
//       -D Trilinos_ENABLE_ifpack:BOOL=ON \
//       -D Trilinos_ENABLE_intrepid:BOOL=ON \
//       -D Trilinos_ENABLE_isorropia:BOOL=ON  \
//       -D Trilinos_ENABLE_kokkos:BOOL=ON \
//       -D Trilinos_ENABLE_komplex:BOOL=ON \
//       -D Trilinos_ENABLE_meros:BOOL=ON \
//       -D Trilinos_ENABLE_mesquite:BOOL=ON \
//       -D Trilinos_ENABLE_ml:BOOL=ON \
//       -D Trilinos_ENABLE_moertel:BOOL=ON \
//       -D Trilinos_ENABLE_moocho:BOOL=ON \
//       -D Trilinos_ENABLE_nox:BOOL=ON \
//       -D Trilinos_ENABLE_optika:BOOL=ON \
//       -D Trilinos_ENABLE_pamgen:BOOL=ON \
//       -D Trilinos_ENABLE_phalanx:BOOL=ON \
//       -D Trilinos_ENABLE_phdmesh:BOOL=ON \
//       -D Trilinos_ENABLE_piro:BOOL=ON \
//       -D Trilinos_ENABLE_pliris:BOOL=ON \
//       -D Trilinos_ENABLE_PyTrilinos:BOOL=ON \
//       -D Trilinos_ENABLE_rtop:BOOL=ON \
//       -D Trilinos_ENABLE_rythmos:BOOL=ON \
//       -D Trilinos_ENABLE_sacado:BOOL=ON \
//       -D Trilinos_ENABLE_shards:BOOL=ON \
//       -D Trilinos_ENABLE_stk:BOOL=ON \
//       -D Trilinos_ENABLE_stokhos:BOOL=ON \
//       -D Trilinos_ENABLE_stratimikos:BOOL=ON \
//       -D Trilinos_ENABLE_Sundance:BOOL=ON \
//       -D Trilinos_ENABLE_teuchos:BOOL=ON \
//       -D Trilinos_ENABLE_ThreadPool:BOOL=ON \
//       -D Trilinos_ENABLE_thyra:BOOL=ON \
//       -D Trilinos_ENABLE_tifpack:BOOL=ON \
//       -D Trilinos_ENABLE_tpetra:BOOL=ON \
//       -D Trilinos_ENABLE_TriKota:BOOL=ON \
//       -D Trilinos_ENABLE_trilinoscouplings:BOOL=ON \
//       -D Trilinos_ENABLE_triutils:BOOL=ON \
//       -D Trilinos_ENABLE_zoltan:BOOL=ON \
//       -D BUILD_SHARED_LIBS:BOOL=ON \
//       -D CMAKE_INSTALL_PREFIX:PATH=/opt/stow/trilinos-10.2.0 \
//       $TRILINOS_HOME

fsolver_path = get_absolute_file_path('builder_gateway_cpp.sce');

// Trilinos-10.0.5
trilinos_base         = '/opt/stow/trilinos-10.2.0/';
trilinos_lib_inst_dir = trilinos_base + '/lib/';
trilinos_inc_inst_dir = trilinos_base + '/include/';

library_name = 'scinox';

ldflags = [];
fflags = ' -ggdb';
cflags  = '-ggdb -D__USE_DEPRECATED_STACK_FUNCTIONS__ -I' + fsolver_path;

table = ['fsolver_nox', 'sci_fsolver_nox'];

files = ['sci_fsolver_nox.cpp','parameters.cpp'];

libs  = [trilinos_lib_inst_dir + 'libteuchos', ...
         trilinos_lib_inst_dir + 'libepetra', ...
         trilinos_lib_inst_dir + 'libtriutils', ...
         trilinos_lib_inst_dir + 'libaztecoo', ...
         trilinos_lib_inst_dir + 'libsimpi', ...
         trilinos_lib_inst_dir + 'libzoltan', ...
         trilinos_lib_inst_dir + 'libepetraext', ...
         trilinos_lib_inst_dir + 'libamesos', ...
         trilinos_lib_inst_dir + 'libifpack', ...
         trilinos_lib_inst_dir + 'libisorropia', ...
         trilinos_lib_inst_dir + 'libgaleri', ...
         trilinos_lib_inst_dir + 'libml', ...
         trilinos_lib_inst_dir + 'libFEApp', ...
         trilinos_lib_inst_dir + 'libsacado', ...
         trilinos_lib_inst_dir + 'librtop', ...
         trilinos_lib_inst_dir + 'libthyra', ...
         trilinos_lib_inst_dir + 'libthyraepetra', ...
         trilinos_lib_inst_dir + 'libthyraepetraext', ...
         trilinos_lib_inst_dir + 'libtpi', ...
         trilinos_lib_inst_dir + 'libkokkos', ...
         trilinos_lib_inst_dir + 'libkokkosnodeapi', ...
         trilinos_lib_inst_dir + 'libkokkoslinalg', ...
         trilinos_lib_inst_dir + 'libtpetra', ...
         trilinos_lib_inst_dir + 'libtpetrainout', ...
         trilinos_lib_inst_dir + 'libanasazi', ...
         trilinos_lib_inst_dir + 'libanasaziepetra', ...
         trilinos_lib_inst_dir + 'libModeLaplace', ...
         trilinos_lib_inst_dir + 'libbelos', ...
         trilinos_lib_inst_dir + 'libbelosepetra', ...
         trilinos_lib_inst_dir + 'libstratimikosamesos', ...
         trilinos_lib_inst_dir + 'libstratimikosaztecoo', ...
         trilinos_lib_inst_dir + 'libstratimikosbelos', ...
         trilinos_lib_inst_dir + 'libstratimikosifpack', ...
         trilinos_lib_inst_dir + 'libstratimikosml', ...
         trilinos_lib_inst_dir + 'libstratimikos', ...
         trilinos_lib_inst_dir + 'libnox', ...
         trilinos_lib_inst_dir + 'libnoxepetra'];

cflags = '-g -Wall -I' + trilinos_inc_inst_dir + ' -I' + fsolver_path;

tbx_build_gateway(library_name, table, files, fsolver_path,libs,ldflags,cflags,fflags);

clear tbx_build_gateway;
