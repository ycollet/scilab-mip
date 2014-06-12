// Copyright (C) 2010 - Yann COLLETTE
//
// This file is released into the public domain

demopath = get_absolute_file_path("scicoinutils.dem.gateway.sce");

subdemolist = ["demo read GMS", "read_gms.sce"; ...
               "demo read MPS", "read_mps.sce"; ...
               "demo write parameters", "write_data_parameters.sce"];

subdemolist(:,2) = demopath + subdemolist(:,2);
