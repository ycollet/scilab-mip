// ====================================================================
// Yann COLLETTE
// Copyright 2009
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path("dace_scilab.dem.gateway.sce");

subdemolist = ["demo Dace Scilab","DACEdemo.sce"];

subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
