// Copyright (C) 2010 - DIGITEO

// This file is released into the public domain

mode(-1);
lines(0);

TOOLBOX_NAME  = "scinlopt";
TOOLBOX_TITLE = "SciNLOpt";
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

// Check development_tools module avaibility
// =============================================================================

if ~with_module('development_tools') then
  error(msprintf(gettext('%s module not installed."),'development_tools'));
end

// ====================================================================
// Build thirdparty first

Build_64Bits = %t;
Download = %t;
Version  = "2.2.3";
UseProxy = %t;
proxy    = "on"; // For some old wget, you must set this to on and add a http_proxy= line into .wgetrc
//proxy    = "10.38.2.10:8080";

if getos()~='Windows' then
  cd('thirdparty');
  if (Download) then
    printf('Downloading NLOPT ... please wait\n\n');
    if getos()=='Darwin' then
      if (UseProxy) then
        unix_w('curl -o nlopt-' + Version + '.tar.gz http://ab-initio.mit.edu/nlopt/nlopt-' + Version + '.tar.gz');
      else
        unix_w('curl -x ' + proxy + ' -o nlopt-' + Version + '.tar.gz http://ab-initio.mit.edu/nlopt/nlopt-' + Version + '.tar.gz');
      end
    else
      if (UseProxy) then
        unix_w('wget --proxy=' + proxy + ' http://ab-initio.mit.edu/nlopt/nlopt-' + Version + '.tar.gz');
      else
        unix_w('wget http://ab-initio.mit.edu/nlopt/nlopt-' + Version + '.tar.gz');
      end
    end
  end
  unix_w('tar xvfz nlopt-' + Version + '.tar.gz');
  cd('nlopt-' + Version);
  printf('\n\nConfiguring NLOPT ... please wait\n\n');
  if Build_64Bits then
    FLAGS='--with-pic';
  else
    FLAGS='';
  end
  unix_w('./configure --prefix=' + toolbox_dir + '/thirdparty/nlopt ' + FLAGS + ' --enable-static --disable-shared --enable-debug --with-cxx');

  printf('\n\nBuilding NLOPT ... please wait\n\n');
  unix_w('make');
  printf('\n\nInstalling NLOPT ... please wait\n\n');
  unix_w('make install');
  cd('../..');
else
  cd('thirdparty/win');
  if (Download) then
    printf('\n\nDownload NLOPT ... please wait\n\n');

    if (UseProxy) then
      unix_w(SCI + '\tools\curl\curl.exe -x ' + proxy + ' http://ab-initio.mit.edu/nlopt/nlopt-' + Version + '-dll.zip -o nlopt-' + Version + '-dll.zip');
    else
      unix_w(SCI + '\tools\curl\curl.exe http://ab-initio.mit.edu/nlopt/nlopt-' + Version + '-dll.zip -o nlopt-' + Version + '-dll.zip');
    end
  end
  // If libnlopt-0.lib is not present in the archive, you must
  // Open a visual studio console, go to the thirdparty/win directory
  // and enter the following command:
  // lib /DEF:libnlopt-0.def /OUT:libnlopt-0.lib
  // Now switch the Download variable to %f and start again the
  // builder.sce

  printf('\n\nUncompressing NLOPT ... please wait\n\n');
  unix_w(SCI + '\tools\zip\unzip.exe -o nlopt-' + Version + '-dll.zip');

  cd('../..');
end

// Action
// =============================================================================

tbx_builder_macros(toolbox_dir);
tbx_builder_gateway(toolbox_dir);
tbx_builder_help(toolbox_dir);
tbx_build_loader(TOOLBOX_NAME, toolbox_dir);
tbx_build_cleaner(TOOLBOX_NAME, toolbox_dir);

// Clean variables
// =============================================================================

clear toolbox_dir TOOLBOX_NAME TOOLBOX_TITLE;

// ====================================================================
// Remove the ipopt build tree + install tree
if getos()~='Windows' then
  printf('\n\nCleaning NLOPT ... please wait\n\n');
  cd('thirdparty');
  unix_w('rm -rf nlopt*');
  cd('../');
end
