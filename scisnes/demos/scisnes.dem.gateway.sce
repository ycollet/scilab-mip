// ====================================================================
// Copyright 2009
// Yann COLLETTE
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path("scisnes.dem.gateway.sce");

subdemolist = ["demo snes","fsolve_snes.sce"i];

subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
