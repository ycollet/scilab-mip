// Copyright (C) 2010 - Yann COLLETTE
//
// This file is released into the public domain

demopath = get_absolute_file_path("sciglpk.dem.gateway.sce");

subdemolist = ["demo GLPK sparse",  "glpksparse.sce"; ...
               "demo GLPK test 1",  "glpktest1.sce"; ...
               "demo GLPK test 2",  "glpktest2.sce"; ...
               "demo MPS and GLPK", "mps_glpk_test.sce"];

subdemolist(:,2) = demopath + subdemolist(:,2);
