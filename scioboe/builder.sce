// ====================================================================
// Copyright Yann COLLETTE 2010
// This file is released into the public domain
// ====================================================================
mode(-1);
lines(0);
try
 getversion('scilab');
catch
 error(gettext('Scilab 5.0 or more is required.'));  
end;
// Uncomment to make a Debug version
//setenv("DEBUG_SCILAB_DYNAMIC_LINK","YES")
// ====================================================================
if ~with_module('development_tools') then
  error(msprintf(gettext('%s module not installed."),'development_tools'));
end
// ====================================================================
TOOLBOX_NAME = 'scioboe';
TOOLBOX_TITLE = 'SciOBOE';
// ====================================================================
toolbox_dir = get_absolute_file_path('builder.sce');
// ====================================================================

tbx_builder_gateway(toolbox_dir);
tbx_builder_help(toolbox_dir);
tbx_build_loader(TOOLBOX_NAME, toolbox_dir);
tbx_build_cleaner(TOOLBOX_NAME, toolbox_dir);

clear toolbox_dir TOOLBOX_NAME TOOLBOX_TITLE;

// ====================================================================
