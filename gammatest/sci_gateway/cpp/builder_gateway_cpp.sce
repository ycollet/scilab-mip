// ====================================================================
// Copyright 2010 - Yann COLLETTE
// This file is released into the public domain
// ====================================================================

Files = ['Filter.cpp','GammaTest.cpp','kdtree_2.cpp', ...
         'kdtree_static.cpp','sci_gammatest.cpp','helper.cpp'];

Table = ['gammatest','sci_gammatest'];

gammatest_path = get_absolute_file_path('builder_gateway_cpp.sce');

libs    = [];
ldflags = [];
cflags  = '-ggdb -I' + gammatest_path;

tbx_build_gateway('gammatest_cpp', Table, Files, gammatest_path, libs, ldflags, cflags);

clear tbx_build_gateway;
