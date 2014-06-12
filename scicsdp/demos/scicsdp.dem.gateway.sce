// Copyright (C) 2010 - Yann COLLETTE
//
// This file is released into the public domain

demopath = get_absolute_file_path("scicsdp.dem.gateway.sce");

subdemolist = ["demo CSDP",            "csdp_demo.sce"; ...
	       "Read / Write SDPA",    "sdpa_read_write_demo.sce"; ...
	       "Read SDPA test files", "demo_sdpa_read.sce"];

subdemolist(:,2) = demopath + subdemolist(:,2);
