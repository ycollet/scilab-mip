// ====================================================================
// Yann COLLETTE
// Copyright 2009
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path("conmin.dem.gateway.sce");

subdemolist = ["demo conmin example 1" ,"example1.sce"; ..
               "demo conmin example 2" ,"example2.sce"     ; ..
               "demo conmin example 3" ,"example3.sce" ; ..
               "demo conmin example 4" ,"example4.sce" ; ];

subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
