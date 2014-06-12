;##############################################################################################################
; Inno Setup Install script for SciGLPK Module
; http://www.jrsoftware.org/isinfo.php
; Yann COLLETTE
; This file is released into the public domain
;##############################################################################################################
; modify this path where is sciglpk directory
#define BinariesSourcePath "Z:\toolboxes_google\sciglpk"
;
#define SCIGLPK_Module_version "1.0"
#define CurrentYear "2010"
#define SCIGLPK_ModuleDirFilename "sciglpk-1.0"
;##############################################################################################################
[Setup]
; Debut Données de base à renseigner suivant version
SourceDir={#BinariesSourcePath}
AppName=SciGLPK Module
AppVerName=SciGLPK Module version 1.0
DefaultDirName={pf}/{#SCIGLPK_ModuleDirFilename}
InfoAfterfile=readme.txt
LicenseFile=license.txt
WindowVisible=true
AppPublisher=Home
BackColorDirection=lefttoright
AppCopyright=Copyright © {#CurrentYear}
Compression=lzma/max
InternalCompressLevel=normal
SolidCompression=true
VersionInfoVersion={#SCIGLPK_Module_version}
VersionInfoCompany=Home
;##############################################################################################################
[Files]
; Add here files that you want to add
Source: loader.sce; DestDir: {app}
Source: license.txt; DestDir: {app}
Source: readme.txt; DestDir: {app}
Source: etc\sciglpk.quit; DestDir: {app}\etc
Source: etc\sciglpk.start; DestDir: {app}\etc
Source: macros\buildmacros.sce; DestDir: {app}\macros
Source: macros\lib; DestDir: {app}\macros
Source: macros\names; DestDir: {app}\macros
Source: macros\*.sci; DestDir: {app}\macros
Source: macros\*.bin; DestDir: {app}\macros
Source: sci_gateway\loader_gateway.sce; DestDir: {app}\sci_gateway
Source: sci_gateway\cpp\loader.sce; DestDir: {app}\sci_gateway\cpp
Source: sci_gateway\cpp\*.dll; DestDir: {app}\sci_gateway\cpp
Source: demos\*.*; DestDir: {app}\demos;
Source: jar\*.*; DestDir: {app}\jar
Source: thirdparty\win\*.*; DestDir: {app}\thirdparty\win Flags: recursesubdirs
;
;##############################################################################################################
