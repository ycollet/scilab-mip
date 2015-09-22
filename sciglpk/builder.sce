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
TOOLBOX_NAME = 'sciglpk';
TOOLBOX_TITLE = 'SciGLPK';
// ====================================================================
toolbox_dir = get_absolute_file_path('builder.sce');
// ====================================================================
// Build thirdparty first

Build_64Bits = %t;
UseProxy = %f;
PROXY = '66.66.66.66:8080';
VERSION = '4.55';

if getos()~='Windows' then
  mkdir('thirdparty');
  cd('thirdparty');
  printf('Downloading GLPK ... please wait\n\n');
  unix_w('wget http://ftp.gnu.org/gnu/glpk/glpk-' + VERSION + '.tar.gz');
  unix_w('tar xvfz glpk-' + VERSION + '.tar.gz');
  cd('glpk-' + VERSION);
  printf('\n\nConfiguring GLPK ... please wait\n\n');
  if Build_64Bits then
    FLAGS='--with-pic';
  else
    FLAGS='';
  end
  unix_w('./configure --prefix=' + toolbox_dir + '/thirdparty/glpk ' + FLAGS + ' --enable-static --disable-shared --enable-debug');

  printf('\n\nBuilding GLPK ... please wait\n\n');
  unix_w('make');
  printf('\n\nInstalling GLPK ... please wait\n\n');
  unix_w('make install');
  cd('../..');
else
  mkdir('thirdparty');
  mkdir('thirdparty/win');
  cd('thirdparty/win');
  printf('Downloading GLPK ... please wait\n\n');
  if UseProxy then
    unix_w(SCI + '\tools\curl\curl.exe -L -x ' + PROXY + ' -o glpk-' + VERSION + '.zip http://sourceforge.net/projects/winglpk/files/winglpk/GLPK-' + VERSION + '/winglpk-' + VERSION + '.zip/download');
  else
    unix_w(SCI + '\tools\curl\curl.exe -L -o glpk-' + VERSION + '.zip http://sourceforge.net/projects/winglpk/files/winglpk/GLPK-' + VERSION + '/winglpk-' + VERSION + '.zip/download');
  end
  
  printf('\n\nUncompressing GLPK ... please wait\n\n');
  unix_w(SCI + '\tools\zip\unzip.exe -o glpk-' + VERSION + '.zip');
  
  cd('glpk-' + VERSION);
  
  if win64() then
    cd('w64');
  else
    cd('w32');
  end
  unix_w('copy config_VC config.h');
  unix_w('copy makefile_VC makefile.mak');
  G_make('','');
  cd('../..');
  
  clear tmp_version;
  
  cd('../..');
end

// ====================================================================
// Store the GLPK version number into a macro
// ====================================================================

fd = mopen('macros/glpk_version.sci','w');
mputl(['function glpk_ver = glpk_version()', ...
       '  glpk_ver = ''' + VERSION + ''';', ...
       'endfunction'], fd);
mclose(fd);
clear fd;

// ====================================================================

tbx_builder_macros(toolbox_dir);
tbx_builder_gateway(toolbox_dir);
tbx_builder_help(toolbox_dir);
tbx_build_loader(TOOLBOX_NAME, toolbox_dir);
tbx_build_cleaner(TOOLBOX_NAME, toolbox_dir);

clear toolbox_dir TOOLBOX_NAME TOOLBOX_TITLE;

// ====================================================================
// Remove the build tree + install tree
if getos()~='Windows' then
  printf('\n\nCleaning GLPK ... please wait\n\n');
  cd('thirdparty');
  unix_w('rm -rf glpk-' + VERSION + '*');
  cd('../');
end
