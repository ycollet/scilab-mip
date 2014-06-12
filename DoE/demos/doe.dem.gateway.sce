// ====================================================================
// Yann COLLETTE
// Copyright 2009
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path("doe.dem.gateway.sce");

subdemolist = ["Quasi random demos",    "QuasiRandomdemo.sce"; ...
               "Classic DoE demos",     "classic_doe_demo.sce"; ...
               "D optimality demo",     "d_opti_demo.sce"; ...
               "Model selection demo",  "model_select_demo.sce"; ...
               "LARS demo",             "lars_test.sce"];

subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
