// ====================================================================
// Copyright Y. Collette - 2012
// This file is released into the public domain
// ====================================================================

freefem_path = get_absolute_file_path('builder_gateway_c.sce');

TABLE = ['ff_exec',     'int_ff_exec'; ...
         'ff_problem',  'int_ff_problem'; ...
         'ff_end',      'int_ff_end'; ...
         'getffResult', 'int_get_ff_result'; ...
         'buildMesh',   'int_put_mesh'; ...
         'UpdateMesh',  'int_put_new_mesh'; ...
         'getMatrix',   'int_get_matrix'];

FILES = ['freefemfi.c'];

if getos()=='Windows' then
  CFLAGS  = '-D__USE_DEPRECATED_STACK_FUNCTIONS__';
  LDFLAGS = '';
else
  CFLAGS  = '-ggdb -D__USE_DEPRECATED_STACK_FUNCTIONS__';
  LDFLAGS = '';
end
  
tbx_build_gateway('freefem_c', TABLE, FILES, freefem_path, ['../../src/c/libdelete_lb_lobj_'], LDFLAGS, CFLAGS);

clear tbx_build_gateway;
