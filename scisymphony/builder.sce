// ====================================================================
// Copyright Yann COLLETTE 2011
// This file is released into the public domain
// ====================================================================
mode(-1);
lines(0);
try
 getversion('scilab');
catch
 error(gettext('Scilab 5.0 or more is required.'));  
end;
// ====================================================================
if ~with_module('development_tools') then
  error(msprintf(gettext('%s module not installed."),'development_tools'));
end
// ====================================================================
TOOLBOX_NAME = 'scisymphony';
TOOLBOX_TITLE = 'SciSymphony';
// ====================================================================
toolbox_dir = get_absolute_file_path('builder.sce');
// ====================================================================

// Build thirdparty first

VERSION  = '5.4.3';
Build_64Bits = %t;
Build_Symphony = %t;
UseProxy = %f;
proxy    = "10.38.2.10:8080";
MSV_PATH = '""c:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\IDE\""';

if Build_Symphony then
  if getos()~='Windows' then
    cd('thirdparty');
    printf('Downloading SYMPHONY ... please wait\n\n');
    if UseProxy then
      unix_w('wget --proxy=' + proxy + 'http://www.coin-or.org/download/source/SYMPHONY/SYMPHONY-' + VERSION + '.tgz');
    else
      unix_w('wget http://www.coin-or.org/download/source/SYMPHONY/SYMPHONY-' + VERSION + '.tgz');
    end
    unix_w('tar xvfz SYMPHONY-' + VERSION + '.tgz');
    cd('SYMPHONY-' + VERSION);
    printf('\n\nConfiguring SYMPHONY ... please wait\n\n');
    if Build_64Bits then
      FLAGS='--with-pic';
    else
      FLAGS='';
    end
    unix_w('./configure --prefix=' + toolbox_dir + '/thirdparty/symphony ' + FLAGS + ' --enable-static --disable-shared --enable-debug');
  
    printf('\n\nBuilding SYMPHONY ... please wait\n\n');
    unix_w('make');
  
    printf('\n\nInstalling SYMPHONY ... please wait\n\n');
    unix_w('make install');
    cd('../..');
  else
    cd('thirdparty/win');
  
    printf('\n\nDownloading SYMPHONY ... please wait\n\n');
    if UseProxy then
      unix_w(SCI + '\tools\curl\curl.exe -x ' + proxy + ' -o SYMPHONY-' + VERSION + '.zip http://www.coin-or.org/download/source/SYMPHONY/SYMPHONY-' + VERSION + '.zip');
    else
      unix_w(SCI + '\tools\curl\curl.exe -o SYMPHONY-' + VERSION + '.zip http://www.coin-or.org/download/source/SYMPHONY/SYMPHONY-' + VERSION + '.zip');
    end
  
    printf('\n\nUncompressing SYMPHONY ... please wait\n\n');
    unix_w(SCI + '\tools\zip\unzip.exe -o SYMPHONY-' + VERSION + '.zip');
  
    printf('\n\nBuilding the solution ... please wait\n\n')
    cd('SYMPHONY-'+VERSION+'\SYMPHONY\MSVisualStudio\v10');
    
    [output, err] = dos(MSV_PATH + '..\\..\\VC\\vcvarsall.bat');
    disp(output);
    [output, err] = dos(MSV_PATH + 'devenv.exe symphony.sln /build release');
    disp(output);
    cd('../../../../../..');
  end
end

tbx_builder_macros(toolbox_dir);
tbx_builder_gateway(toolbox_dir);
tbx_builder_help(toolbox_dir);
tbx_build_loader(TOOLBOX_NAME, toolbox_dir);
tbx_build_cleaner(TOOLBOX_NAME, toolbox_dir);

clear toolbox_dir TOOLBOX_NAME TOOLBOX_TITLE;
// ====================================================================
// Remove the build tree + install tree
if getos()~='Windows' then
  printf('\n\nCleaning SYMPHONY ... please wait\n\n');
  cd('thirdparty');
  unix_w('rm -rf symphony');
  unix_w('rm -rf SYMPHONY*');
  cd('../');
end
