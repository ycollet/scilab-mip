// ====================================================================
// Yann COLLETTE
// Copyright 2009
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path("lolimot.dem.gateway.sce");

subdemolist = ["demo lolimot" ,"testLolimot.sce"];

subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
