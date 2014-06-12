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
TOOLBOX_NAME = 'scilpsolve';
TOOLBOX_TITLE = 'SciLPSolve';
// ====================================================================
toolbox_dir = get_absolute_file_path('builder.sce');
// ====================================================================
// Build thirdparty first

UseProxy = %t;
VERSION = '5.5.2.0';
PROXY = '10.38.22.2:8080';
// http://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.0.15/lp_solve_5.5.0.15_dev_win64.zip/download

if getos()=='Linux' then
  cd('thirdparty');
  printf('Downloading LPSolve ... please wait\n\n');
  unix_w('wget http://sourceforge.net/projects/lpsolve/files/lpsolve/' + VERSION + '/lp_solve_' + VERSION + '_source.tar.gz/download');
  unix_w('mv download lp_solve_' + VERSION + '_source.tar.gz');
  unix_w('tar xvfz lp_solve_' + VERSION + '_source.tar.gz');
  cd('lp_solve_5.5');
  if Build_64Bits then
    CFLAGS='-fPIC';
  end
  printf('\n\nBuilding LPSolve ... please wait\n\n');
  //cd('lp_solve');
  //unix_w('sh ccc');
  cd('lpsolve55');
  if Build_64Bits then
    unix('sed -e ""s/INVERSE_LUSOL/INVERSE_LUSOL -fPIC/g"" < ccc > ccc.scilab');
  else
    unix('cp ccc ccc.scilab');
  end
  unix_w('sh ccc.scilab');
  cd('../../..');
elseif getos()=='Darwin' then
  cd('thirdparty');
  printf('Downloading LPSolve ... please wait\n\n');
  unix_w('wget http://sourceforge.net/projects/lpsolve/files/lpsolve/' + VERSION + '/lp_solve_' + VERSION + '_source.tar.gz/download');
  unix_w('tar xvfz lp_solve_' + VERSION + '_source.tar.gz');
  cd('lp_solve_5.5');
  printf('\n\nBuilding LPSolve ... please wait\n\n');
  cd('lpsolve55');
  unix_w('sh ccc.osx');
  cd('../../..');
else
  cd('thirdparty/win');
  if win64() then
    printf('\n\nDownloading LPSolve ... please wait\n\n');
    if UseProxy then
      unix_w(SCI + '\tools\curl\curl.exe -L -x ' + PROXY + ' -o lp_solve_' + VERSION + '_dev_win64.zip http://sourceforge.net/projects/lpsolve/files/lpsolve/' + VERSION + '/lp_solve_' + VERSION + '_dev_win64.zip/download');
    else
      unix_w(SCI + '\tools\curl\curl.exe -L -o lp_solve_' + VERSION + '_dev_win64.zip http://sourceforge.net/projects/lpsolve/files/lpsolve/' + VERSION + '/lp_solve_' + VERSION + '_dev_win64.zip/download');
    end

    printf('\n\nUncompressing LPSolve ... please wait\n\n');
    unix_w(SCI + '\tools\zip\unzip.exe -o lp_solve_' + VERSION + '_dev_win64.zip');
    unix_w('del lp_solve_' + VERSION + '_dev_win64.zip');
  else
    printf('\n\nDownloading LPSolve ... please wait\n\n');
    if UseProxy then
      unix_w(SCI + '\tools\curl\curl.exe -L -x ' + PROXY + ' -o lp_solve_' + VERSION + '_dev_win32.zip http://sourceforge.net/projects/lpsolve/files/lpsolve/' + VERSION + '/lp_solve_' + VERSION + '_dev_win32.zip/download');
    else
      unix_w(SCI + '\tools\curl\curl.exe -L -o lp_solve_' + VERSION + '_dev_win32.zip http://sourceforge.net/projects/lpsolve/files/lpsolve/' + VERSION + '/lp_solve_' + VERSION + '_dev_win32.zip/download');
    end
    printf('\n\nUncompressing LPSolve ... please wait\n\n');
    unix_w(SCI + '\tools\zip\unzip.exe -o lp_solve_' + VERSION + '_dev_win32.zip');
    unix_w('del lp_solve_' + VERSION + '_dev_win32.zip');
  end
  cd('../..');
end

// ====================================================================

tbx_builder_macros(toolbox_dir);
tbx_builder_gateway(toolbox_dir);
tbx_builder_help(toolbox_dir);
tbx_build_loader(TOOLBOX_NAME, toolbox_dir);
tbx_build_cleaner(TOOLBOX_NAME, toolbox_dir);

clear toolbox_dir TOOLBOX_NAME TOOLBOX_TITLE;

// ====================================================================
// Remove the ipopt build tree + install tree
if getos()~='Windows' then
  printf('\n\nCleaning LPSolve ... please wait\n\n');
  cd('thirdparty');
  unix_w('rm -rf lp_solve_5.5*');
  cd('../');
end
