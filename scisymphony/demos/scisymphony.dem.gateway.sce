// ====================================================================
// Copyright Yann COLLETTE 2011
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path("scisymphony.dem.gateway.sce");

subdemolist = ["demo of symphony", "mps_symphony_test.sce"];

subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
