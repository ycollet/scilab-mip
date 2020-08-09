// ====================================================================
// Yann COLLETTE
// Copyright 2012
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path('scifreefem.dem.gateway.sce');

subdemolist = ['FreeFem to Scilab : PDE''s',                'ffi_ex1.sce'; ..
               'FreeFem to Scilab : Variational equations', 'ffi_ex2.sce'; ..
               'Scilab to Freefem : Moving border',         'ffi_ex3.sce'; ];

subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
