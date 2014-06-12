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
TOOLBOX_NAME = 'scicoinutils';
TOOLBOX_TITLE = 'SciCoinUtils';
// ====================================================================
toolbox_dir = get_absolute_file_path('builder.sce');
// ====================================================================
// Build thirdparty first

Build        = %t;
Build_64Bits = %f;
VERSION  = '2.8.2';
UseProxy = %t;
proxy    = "on"; // For some old wget, you must set this to on and add a http_proxy= line into .wgetrc
//proxy    = "10.38.2.10:8080";
MSV_PATH = '""c:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\IDE\""';

current_path = pwd();

if Build then
    if getos()~='Windows' then
      cd('thirdparty');
      printf('Downloading CoinUtils ... please wait\n\n');
      if (UseProxy) then
        unix_w('wget --proxy=" + proxy + " http://www.coin-or.org/download/source/CoinUtils/CoinUtils-' + VERSION + '.tgz');
      else
        unix_w('wget http://www.coin-or.org/download/source/CoinUtils/CoinUtils-' + VERSION + '.tgz');
      end
      unix_w('tar xvfz CoinUtils-' + VERSION + '.tgz');
      cd('CoinUtils-' + VERSION);
      printf('\n\nConfiguring CoinUtils ... please wait\n\n');
      if Build_64Bits then
        FLAGS='--with-pic';
      else
        FLAGS='';
      end
      
      unix_w('./configure --prefix=' + toolbox_dir + '/thirdparty/coinutils ' + FLAGS + ' --enable-static --disable-shared --enable-debug');
    
      printf('\n\nBuilding CoinUtils ... please wait\n\n');
      unix_w('make');
    
      printf('\n\nInstalling CoinUtils ... please wait\n\n');
      unix_w('make install');
      cd('../..');
    else
      cd('thirdparty/win');
    
      printf('\n\nDownloading CoinUtils ... please wait\n\n');
      if (UseProxy) then
        unix_w(SCI + '\tools\curl\curl.exe -x ' + proxy + ' http://www.coin-or.org/download/source/CoinUtils/CoinUtils-' + VERSION + '.zip -o CoinUtils-' + VERSION + '.zip');
      else
        unix_w(SCI + '\tools\curl\curl.exe -o CoinUtils-' + VERSION + '.zip http://www.coin-or.org/download/source/CoinUtils/CoinUtils-' + VERSION + '.zip');
      end
    
      printf('\n\nUncompressing CoinUtils ... please wait\n\n');
      unix_w(SCI + '\tools\zip\unzip.exe -o CoinUtils-' + VERSION + '.zip');

      printf('\n\nBuilding the solution ... please wait\n\n')
      cd('CoinUtils-' + VERSION + '\CoinUtils\MSVisualStudio\v10');
      dos(MSV_PATH + '..\..\VC\vcvarsall.bat');
      dos(MSV_PATH + 'devenv.exe CoinUtils.sln /build Release');
    end
    cd(current_path);
end


// ====================================================================

tbx_builder_gateway(toolbox_dir);
tbx_builder_help(toolbox_dir);
tbx_build_loader(TOOLBOX_NAME, toolbox_dir);
tbx_build_cleaner(TOOLBOX_NAME, toolbox_dir);

clear toolbox_dir TOOLBOX_NAME TOOLBOX_TITLE;

// ====================================================================
// Remove the build tree + install tree
if getos()~='Windows' then
  printf('\n\nCleaning CoinUtils ... please wait\n\n');
  cd('thirdparty');
  unix_w('rm -rf coinutils');
  unix_w('rm -rf CoinUtils*');
  cd('../');
end
