;##############################################################################################################
; Inno Setup Install script for Toolbox_skeleton
; http://www.jrsoftware.org/isinfo.php
; Yann COLLETTE
; This file is released into the public domain
;##############################################################################################################
; modify this path where is toolbox_skeleton directory
#define BinariesSourcePath "E:\toolboxes\ampli-1.3"
;
#define AMPLToolboxversion "1.3"
#define CurrentYear "2010"
#define AMPLToolboxDirFilename "scilab-ampl-1.3"
;##############################################################################################################
[Setup]
; Debut Données de base à renseigner suivant version
SourceDir={#BinariesSourcePath}
AppName=AMPL Toolbox
AppVerName={#AMPLToolboxversion}
DefaultDirName={pf}/{#AMPLToolboxDirFilename}
InfoAfterfile=readme.txt
LicenseFile=license.txt
WindowVisible=true
AppPublisher=Home
BackColorDirection=lefttoright
AppCopyright=Copyright © {#CurrentYear}
Compression=lzma/max
InternalCompressLevel=normal
SolidCompression=true
VersionInfoVersion={#AMPLToolboxversion}
VersionInfoCompany=Home
;##############################################################################################################
[Files]
; Add here files that you want to add
Source: loader.sce; DestDir: {app}
Source: license.txt; DestDir: {app}
Source: readme.txt; DestDir: {app}
Source: changelog.txt; DestDir: {app}
Source: etc\ampl_toolbox.quit; DestDir: {app}\etc
Source: etc\ampl_toolbox.start; DestDir: {app}\etc
Source: macros\buildmacros.sce; DestDir: {app}\macros
Source: macros\lib; DestDir: {app}\macros
Source: macros\names; DestDir: {app}\macros
Source: macros\*.sci; DestDir: {app}\macros
Source: macros\*.bin; DestDir: {app}\macros
Source: sci_gateway\loader_gateway.sce; DestDir: {app}\sci_gateway
Source: sci_gateway\c\loader.sce; DestDir: {app}\sci_gateway\c
Source: sci_gateway\c\ampl_c.dll; DestDir: {app}\sci_gateway\c
;Source: sci_gateway\cpp\loader.sce; DestDir: {app}\sci_gateway\cpp
;Source: sci_gateway\cpp\skeleton_cpp.dll; DestDir: {app}\sci_gateway\cpp
;Source: sci_gateway\fortran\loader.sce; DestDir: {app}\sci_gateway\fortran
;Source: sci_gateway\fortran\skeleton_fortran.dll; DestDir: {app}\sci_gateway\fortran
;Source: src\c\libcsum.dll; DestDir: {app}\src\c
;Source: src\c\loader.sce; DestDir: {app}\src\c
;Source: src\fortran\libfsum.dll; DestDir: {app}\src\fortran
;Source: src\fortran\loader.sce; DestDir: {app}\src\fortran
;Source: tests\*.*; DestDir: {app}\tests; Flags: recursesubdirs
;Source: includes\*.h; DestDir: {app}\includes; Flags: recursesubdirs
;Source: locales\*.*; DestDir: {app}\locales; Flags: recursesubdirs
Source: demos\*.*; DestDir: {app}\demos
Source: demos\data\*.*; DestDir: {app}\demos\data; Flags: recursesubdirs
Source: jar\*.*; DestDir: {app}\jar
;
;##############################################################################################################
