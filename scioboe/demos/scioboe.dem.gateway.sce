// Copyright (C) 2010 - Yann COLLETTE
//
// This file is released into the public domain

demopath = get_absolute_file_path("scioboe.dem.gateway.sce");

subdemolist = ["demo OBOE", "oboe_demo.sce"];

subdemolist(:,2) = demopath + subdemolist(:,2);
