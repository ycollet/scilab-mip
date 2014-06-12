// ====================================================================
// Copyright 2009
// Yann COLLETTE
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path("newuoa.dem.gateway.sce");

subdemolist = ["demo parabola",  "parabola.sce"; ...
               "demo chebyquad", "chebyquad.sce"];

subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
