// ====================================================================
// Copyright DIGITEO 2010
// Yann COLLETTE
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path("scinlopt.dem.gateway.sce");

subdemolist = ["demo nlopt unconstrained",             "NLOpt_unconstr.sce"; ...
	       "demo nlopt constrained",               "NLOpt_constr.sce"; ...
               "demo nlopt C functions unconstrained", "NLOpt_c_demos.sce"];
subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
