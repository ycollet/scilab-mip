;##############################################################################################################
; Inno Setup Install script for Toolbox_skeleton
; http://www.jrsoftware.org/isinfo.php
; Allan CORNET
; This file is released into the public domain
;##############################################################################################################
; modify this path where is toolbox_skeleton directory
#define BinariesSourcePath "D:\Scilab5\modules\Lolimot-2.0"
;
#define Toolbox_skeleton_version "2.0"
#define CurrentYear "2009"
#define Toolbox_skeletonDirFilename "Lolimot"
;##############################################################################################################
[Setup]
; Debut Données de base à renseigner suivant version
SourceDir={#BinariesSourcePath}
AppName=Toolbox skeleton
AppVerName=Lolimot version 2.0
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
VersionInfoCompany=Your Company
;##############################################################################################################
[Files]
; Add here files that you want to add
Source: loader.sce; DestDir: {app}
Source: etc\lolimot.quit; DestDir: {app}\etc
Source: etc\lolimot.start; DestDir: {app}\etc
Source: macros\buildmacros.sce; DestDir: {app}\macros
Source: macros\lib; DestDir: {app}\macros
Source: macros\names; DestDir: {app}\macros
Source: macros\*.sci; DestDir: {app}\macros
Source: macros\*.bin; DestDir: {app}\macros
Source: tests\*.*; DestDir: {app}\tests; Flags: recursesubdirs
Source: demos\*.*; DestDir: {app}\locales; Flags: recursesubdirs
Source: help\*.*; DestDir: {app}\help; Flags: recursesubdirs
;
;##############################################################################################################
