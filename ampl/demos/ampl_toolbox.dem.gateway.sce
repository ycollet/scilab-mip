// ====================================================================
// Yann COLLETTE
// Copyright 2009
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path("ampl_toolbox.dem.gateway.sce");

subdemolist = ["demo AMPL Toolbox", "demo_nl.sce"];

subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
