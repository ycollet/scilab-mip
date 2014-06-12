// ====================================================================
// Copyright 2009
// Yann COLLETTE
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path("tsp_toolbox.dem.gateway.sce");

subdemolist = ["demo TSP and SA"             ,"SAdemo.sce";]; ..

subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
