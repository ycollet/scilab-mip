// ====================================================================
// Copyright 2010 - Yann COLLETTE
// This file is released into the public domain
// ====================================================================
libpath = get_absolute_file_path('cleanmacros.sce');

binfiles = ls(libpath+'/*.bin');
for i = 1:size(binfiles,'*')
  mdelete(binfiles(i));
end

mdelete(libpath+'/names');
mdelete(libpath+'/lib');

// ====================================================================
