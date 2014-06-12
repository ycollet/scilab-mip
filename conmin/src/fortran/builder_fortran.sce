// ====================================================================
// Yann COLLETTE
// Copyright 2009
// This file is released into the public domain
// ====================================================================
libs = '../c/libcommon_conmin';

names = ['get_ct','get_obj','get_igoto','get_nac',...
         'set_delfun','set_dabfun','set_fdch','set_fdchm','set_ct','set_ctmin','set_ctl','set_ctlmin', ...
	 'set_alphax','set_abobj1','set_theta','set_obj','set_ndv','set_ncon','set_nside','set_iprint', ...
	 'set_nfdg','set_nscal','set_linobj','set_itmax','set_itrm','set_icndir','set_igoto','set_nac', ...
	 'set_info','set_infog','set_iter',...
	 'conmin'];

conmin_path = get_absolute_file_path('builder_fortran.sce')

files = 'conmin_scilab.f';

if getos()=='Windows' then
  cflags  = '/I' + SCI + '/modules/core/includes';
  ldflags = '';

  files = [files, 'lnblnk.f'];
else
  cflags  = '-I' + SCI + '/modules/core/includes';
  ldflags = '';
end
	
tbx_build_src(names, files, 'f', conmin_path, [], ldflags, cflags);
clear tbx_build_src;
