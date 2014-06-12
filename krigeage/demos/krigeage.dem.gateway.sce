// ====================================================================
// Yann COLLETTE
// Copyright 2009
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path("krigeage.dem.gateway.sce");

subdemolist = ["demo Krigeage","krigeage.dem.sce"];

subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
