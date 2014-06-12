// ====================================================================
// Yann COLLETTE
// Copyright 2009-2010
// This file is released into the public domain
// ====================================================================

if MSDOS then
  rep = unix_g('nmake /Y /nologo /f makefile.vc');
else
  rep = unix_g('make -f makefile.u');
end

