// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) 2012 - Yann COLLETTE
//
// This file is released under the 3-clause BSD license. See COPYING-BSD.

function subdemolist = demo_gateway()
  demopath = get_absolute_file_path("scifilterSD.dem.gateway.sce");

  subdemolist = ["The dense version of the HS106 problem"   ,"hs106_dense.sce"; ..
                 "The sparse version of the hs106 problem"  ,"hs106_sparse.sce"; ..
                 "The dense version of the hs72 problem"    ,"hs72_dense.sce" ; ..
                 "The sparse version of the hs72 problem"   ,"hs72_sparse.sce" ];

  subdemolist(:,2) = demopath + subdemolist(:,2);
  
endfunction

subdemolist = demo_gateway();
clear demo_gateway; // remove demo_gateway on stack
