// ====================================================================
// Copyright 2009
// Yann COLLETTE
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path("scinox.dem.gateway.sce");

subdemolist = ["demo nox","fsolve_nox.sce"i];

subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
