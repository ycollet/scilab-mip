// ====================================================================
// Yann COLLETTE
// Copyright 2012
// This file is released into the public domain
// ====================================================================
mode(-1);
lines(0);
try
 v = getversion('scilab');
catch
 error(gettext('Scilab 5.3 or more is required.'));  
end

if v(2) < 3 then
  // new API in scilab 5.3
  error(gettext('Scilab 5.3 or more is required.'));
end

// ====================================================================
if ~with_module('development_tools') then
  error(msprintf(gettext('%s module not installed."),'development_tools'));
end
// ====================================================================
TOOLBOX_NAME = 'SciFreeFEM';
TOOLBOX_TITLE = 'SciFreeFEM';
// ====================================================================
toolbox_dir = get_absolute_file_path('builder.sce');

tbx_builder_macros(toolbox_dir);
tbx_builder_src(toolbox_dir);
tbx_builder_gateway(toolbox_dir);
tbx_builder_help(toolbox_dir);
tbx_build_loader(TOOLBOX_NAME, toolbox_dir);

clear toolbox_dir TOOLBOX_NAME TOOLBOX_TITLE;
// ====================================================================
