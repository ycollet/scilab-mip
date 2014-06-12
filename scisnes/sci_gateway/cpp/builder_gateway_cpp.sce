// ====================================================================
// Yann COLLETTE
// Christopher MAES
// Copyright 2009-2012
// This file is released into the public domain
// ====================================================================

// Petsc must be compiled as a shared library and without mpi:
// ./configure --prefix=myinstalldir --enable-shared --without-mpi

fsolver_path = get_absolute_file_path('builder_gateway_cpp.sce');

// Petsc-3
petsc_install_dir = '/opt/stow/petsc-3.1-p1/';
petsc_lib         = petsc_install_dir + '/lib/';
petsc_inc         = petsc_install_dir + '/include/';

//openmpi_lib_install_dir = '/usr/lib/'
//openmpi_inc_install_dir = '/usr/include/';
// If compiled without mpi support, petsc add a uniprocessor mpitool. We must include this ...
openmpi_inc_install_dir = petsc_inc + 'mpiuni/';

library_name = 'scisnes';

table = ['fsolver_snes', 'sci_fsolver_snes']; 

files = ['sci_fsolver_snes.cpp','parameters.cpp'];

fflags  = '';
ldflags = [];

// for petsc-3.1
libs  = [petsc_lib + 'libpetsc'];

//// For petsc-3.0
//libs  = [petsc_lib + 'libpetsc', ...
//	 petsc_lib + 'libpetscvec', ...
//	 petsc_lib + 'libpetscmat', ...
//	 petsc_lib + 'libpetscdm', ...
//	 petsc_lib + 'libpetscksp', ...
//	 petsc_lib + 'libpetscsnes', ...
//	 petsc_lib + 'libpetscts', ...
//	 petsc_lib + 'libpetsccontrib'];

cflags = '-g -D__USE_DEPRECATED_STACK_FUNCTIONS__ -I' + petsc_inc + ' -I' + openmpi_inc_install_dir + ' -I' + fsolver_path;

tbx_build_gateway(library_name, table, files, fsolver_path, libs, ldflags, cflags, fflags);

clear tbx_build_gateway;
