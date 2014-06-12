// ====================================================================
// Copyright Yann COLLETTE 2009
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path("scicoinor.dem.gateway.sce");

subdemolist = ["demo of clp",       "mps_clp_test.sce"; ..
               "demo of cbc",       "mps_cbc_test.sce"; ..
               "demo of osi",       "mps_osi_test.sce"; ..
               "demo of bonmin",    "bonmin_demo.sce"; ..
               "demo of ipopt",     "ipopt_demo.sce"; ..
               "demo of optim_slp", "SLPdemo.sce"];

subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
