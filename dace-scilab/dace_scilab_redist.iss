;##############################################################################################################
; Inno Setup Install script for Toolbox_skeleton
; http://www.jrsoftware.org/isinfo.php
; Allan CORNET
; This file is released into the public domain
;##############################################################################################################
; modify this path where is toolbox_skeleton directory
#define BinariesSourcePath "D:\Scilab-5\modules\Dace-scilab-2.0"
;
#define Toolbox_skeleton_version "2.0"
#define CurrentYear "2009"
#define Toolbox_skeletonDirFilename "Dace_scilab"
;##############################################################################################################
[Setup]
; Debut Données de base à renseigner suivant version
SourceDir={#BinariesSourcePath}
AppName=Toolbox skeleton
AppVerName=Dace Scilab version 2.0
DefaultDirName={pf}\{#Toolbox_skeletonDirFilename}
InfoAfterfile=readme.txt
LicenseFile=license.txt
WindowVisible=true
AppPublisher=Home made
BackColorDirection=lefttoright
AppCopyright=Copyright © {#CurrentYear}
Compression=lzma/max
InternalCompressLevel=normal
SolidCompression=true
VersionInfoVersion={#Toolbox_skeleton_version}
VersionInfoCompany=Home made
;##############################################################################################################
[Files]
; Add here files that you want to add
Source: loader.sce; DestDir: {app}
Source: etc\dace_scilab.quit; DestDir: {app}\etc
Source: etc\dace_scilab.start; DestDir: {app}\etc
Source: macros\buildmacros.sce; DestDir: {app}\macros
Source: macros\lib; DestDir: {app}\macros
Source: macros\names; DestDir: {app}\macros
Source: macros\*.sci; DestDir: {app}\macros
Source: macros\*.bin; DestDir: {app}\macros
Source: tests\*.*; DestDir: {app}\tests; Flags: recursesubdirs
Source: demos\*.*; DestDir: {app}\locales; Flags: recursesubdirs
;
;##############################################################################################################
