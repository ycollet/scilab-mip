// ====================================================================
// Yann COLLETTE
// Copyright 2009
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path("simplex.dem.gateway.sce");

subdemolist = ["demo optimization via Nelder and Mead"  ,"nm_optim.sce"; ...
               "demo step by step via Nelder and Mead"  ,"nm_step.sce"];

subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
