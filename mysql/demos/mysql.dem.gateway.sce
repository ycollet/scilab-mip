// ====================================================================
// Yann COLLETTE
// Copyright 2009
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path("mysql.dem.gateway.sce");

subdemolist = ["demo MySQL 1","sql_1.sce"; ...
               "demo MySQL 2","sql_2.sce"; ...
               "demo MySQL 3","sql_3.sce"; ...
               "demo MySQL 4","sql_4.sce"];

subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
