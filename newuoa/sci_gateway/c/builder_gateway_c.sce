// ====================================================================
// Copyright 2009-2012
// Yann COLLETTE
// This file is released into the public domain
// ====================================================================

newuoa_path = get_absolute_file_path('builder_gateway_c.sce');

ilib_name  = 'libcpowell' ; // interface library name 

// objects files 

files = ['newuoa','bigden','biglag','newuob','trsapp','update','cintpowell'];

libs  = []; // other libs needed for linking

if MSDOS then
  cflags = '/D__USE_DEPRECATED_STACK_FUNCTIONS__ /I' + newuoa_path;
else
  cflags = '-D__USE_DEPRECATED_STACK_FUNCTIONS__ -I' + newuoa_path;
end

ldflags = '';

// table of (scilab_name,interface-name) 
// for fortran coded interface use 'C2F(name)'

table =['newuoa', 'cintpowell'];
	    
tbx_build_gateway(ilib_name, table, files, get_absolute_file_path('builder_gateway_c.sce'),libs,ldflags,cflags);

clear tbx_build_gateway;
