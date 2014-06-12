// ====================================================================
// Copyright 2010
// Yann COLLETTE
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path("gammatest.dem.gateway.sce");

subdemolist = ["demo Gamma Test","demo_gammatest.sce"]

subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
