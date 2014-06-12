// Copyright (C) 2010 - Yann COLLETTE
//
// This file is released into the public domain

demopath = get_absolute_file_path("scilpsolve.dem.gateway.sce");

subdemolist = ["demo lpsolve",         "lpsolve_test.sce"; ...
               "demo SOS lpsolve",     "lpsolve_sos.sce"; ...
               "demo MPS and lpsolve", "mps_lpsolve_test.sce"];

subdemolist(:,2) = demopath + subdemolist(:,2);
