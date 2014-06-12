// ====================================================================
// Yann COLLETTE
// Copyright 2009
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
TOOLBOX_NAME = 'fann_toolbox';
TOOLBOX_TITLE = 'FANN Toolbox';
// ====================================================================
toolbox_dir = get_absolute_file_path('builder.sce');
// ====================================================================
// Build thirdparty first

// You can donwload and install the FANN srcs manually.
// uncompress the archive of the 2.2.0 version in thirdparty.
// go into FANN-2.2.0, make a 'build' directory, go into that directory.
// Enter the cmake command:
// - For 64 bits linux:
// cmake -DCMAKE_INSTALL_PREFIX=../../../fann-bin -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=OFF -DCMAKE_C_FLAGS_RELEASE=-fPIC ..
// - For 32 bits linux:
// cmake -DCMAKE_INSTALL_PREFIX=../../../fann-bin -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=OFF ..
// Now make and install:
// make
// make install

// Options:

Build_64bits = %t; // Must be set to %t if your OS is a 64 bits one.
Download_Srcs = %f; // If %f, use the sources installed in thirdparty/fann-bin.
Build_Srcs = %f; // If sources are already installed in thirdparty dir, we can set this to %f.
UseProxy = %f; // If you are behind a proxy, set this to %t and set the proxy address below.
VERSION = '2.2.0'; // version of fann to be used.
PROXY = '10.38.22.2:8080';
MSV_PATH = '""c:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\IDE\""';
// http://sourceforge.net/projects/fann/files/fann/2.2.0/FANN-2.2.0-Source.zip/download

if (Build_Srcs) then
  if getos()=='Linux' | getos()=='Darwin' then
    cd('thirdparty');
    if (Download_Srcs) then
      printf('Downloading FANN ... please wait\n\n');
      unix_w('wget http://sourceforge.net/projects/fann/files/fann/' + VERSION + '/FANN-' + VERSION + '-Source.tar.gz/download');
      unix_w('mv download FANN-' + VERSION + '-Source.tar.gz');
      unix_w('tar xvfz FANN-' + VERSION + '-Source.tar.gz');
      unix_w('rm FANN-' + VERSION + '-Source.tar.gz');
    end
    // Overwrite the current CMakeLists.txt to allow static compilation
    unix_w('cp CMakeLists.txt FANN-' + VERSION + '-Source/src');
    cd('FANN-' + VERSION + '-Source');

    printf('\n\nBuilding FANN ... please wait\n\n');
    unix_w('mkdir build');
    unix_w('cd build');
    build_cmd = 'cmake -DCMAKE_INSTALL_PREFIX=' + toolbox_dir + '/thirdparty/fann-bin -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=OFF';
    if (Build_64bits) then
      build_cmd = build_cmd + ' -DCMAKE_C_FLAGS_RELEASE=-fPIC ..';
    else
      build_cmd = build_cmd + ' ..';
    end
    unix_w(build_cmd);
    unix_w('make');
    unix_w('make install');
    cd('..');
  else
    cd('thirdparty');
    if (Download_Srcs) then
      printf('\n\nDownloading FANN ... please wait\n\n');
      if UseProxy then
        unix_w(SCI + '\tools\curl\curl.exe -L -x ' + PROXY + ' -o FANN-' + VERSION + '-Source.zip http://sourceforge.net/projects/fann/files/fann/' + VERSION + '/FANN-' + VERSION + '-Source.zip/download');
      else
        unix_w(SCI + '\tools\curl\curl.exe -L -o FANN-' + VERSION + '-Source.zip http://sourceforge.net/projects/fann/files/fann/' + VERSION + '/FANN-' + VERSION + '-Source.zip/download');
      end
      
      printf('\n\nUncompressing FANN ... please wait\n\n');
      unix_w(SCI + '\tools\zip\unzip.exe -o FANN-' + VERSION + '-Source.zip');
      unix_w('del FANN-' + VERSION + '-Source.zip');
    end

    printf('\n\nBuilding the solution ... please wait\n\n')
    cd('FANN-' + VERSION + '-Source\VS2010');
    dos(MSV_PATH + '..\..\VC\vcvarsall.bat');
    dos(MSV_PATH + 'devenv.exe fann.sln /build Release');
    cd('../..');
  end
end

// ====================================================================

//tbx_builder_macros(toolbox_dir);
//tbx_builder_src(toolbox_dir);
tbx_builder_gateway(toolbox_dir);
tbx_builder_help(toolbox_dir);
tbx_build_loader(TOOLBOX_NAME, toolbox_dir);

clear toolbox_dir TOOLBOX_NAME TOOLBOX_TITLE build_cmd;
// ====================================================================
