// ====================================================================
// Yann COLLETTE
// Copyright 2012
// This file is released into the public domain
// ====================================================================

mode(-1);
lines(0);

function main_builder()

  TOOLBOX_NAME = 'doe';
  TOOLBOX_TITLE = 'DoE';
  toolbox_dir   = get_absolute_file_path("builder.sce");

  // Check Scilab's version
  // =============================================================================

  try
    v = getversion("scilab");
  catch
    error(gettext("Scilab 5.3 or more is required."));
  end

  if v(2) < 3 then
    // new API in scilab 5.3
    error(gettext('Scilab 5.3 or more is required.'));
  end

  // Check modules_manager module availability
  // =============================================================================

  if ~isdef('tbx_build_loader') then
    error(msprintf(gettext('%s module not installed."), 'modules_manager'));
  end
  
  // Action
  // =============================================================================
  
  tbx_builder_macros(toolbox_dir);
  //tbx_builder_src(toolbox_dir);
  //tbx_builder_gateway(toolbox_dir);
  tbx_builder_help(toolbox_dir);
  tbx_build_loader(TOOLBOX_NAME, toolbox_dir);
  tbx_build_cleaner(TOOLBOX_NAME, toolbox_dir);
  
endfunction

// =============================================================================
main_builder();
clear main_builder; // remove main_builder on stack
// =============================================================================
