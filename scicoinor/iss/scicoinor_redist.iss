;##############################################################################################################
; Inno Setup Install script for Scicoinor Toolbox
; http://www.jrsoftware.org/isinfo.php
; Yann COLLETTE
; This file is released into the public domain
;##############################################################################################################
; modify this path where is toolbox_skeleton directory
#define BinariesSourcePath "Z:\toolboxes_google\scicoinor"
;
#define SCICOINOR_Toolbox_version "1.2"
#define CurrentYear "2009"
#define SCICOINOR_ToolboxDirFilename "scicoinor-1.2"
;##############################################################################################################
[Setup]
; Debut Données de base à renseigner suivant version
SourceDir={#BinariesSourcePath}
AppName=SciCoinOR Toolbox
AppVerName=SciCoinOR Toolbox version 1.2
DefaultDirName={pf}/{#SCICOINOR_ToolboxDirFilename}
InfoAfterfile=readme.txt
LicenseFile=license.txt
WindowVisible=true
AppPublisher=Home
BackColorDirection=lefttoright
AppCopyright=Copyright © {#CurrentYear}
Compression=lzma/max
InternalCompressLevel=normal
SolidCompression=true
VersionInfoVersion={#SCICOINOR_Toolbox_version}
VersionInfoCompany=Home
;##############################################################################################################
[Files]
; Add here files that you want to add
Source: loader.sce; DestDir: {app}
Source: license.txt; DestDir: {app}
Source: readme.txt; DestDir: {app}
Source: changelog.txt; DestDir: {app}
Source: etc\scicoinor.quit; DestDir: {app}\etc
Source: etc\scicoinor.start; DestDir: {app}\etc
Source: macros\buildmacros.sce; DestDir: {app}\macros
Source: macros\lib; DestDir: {app}\macros
Source: macros\names; DestDir: {app}\macros
Source: macros\*.sci; DestDir: {app}\macros
Source: macros\*.bin; DestDir: {app}\macros
Source: sci_gateway\loader_gateway.sce; DestDir: {app}\sci_gateway
Source: sci_gateway\cpp\loader.sce; DestDir: {app}\sci_gateway\cpp
Source: sci_gateway\cpp\sci_coinor.dll; DestDir: {app}\sci_gateway\cpp
Source: src\win32\lapackpp\lib\lapackpp.dll; DestDir: {app}\src\win32\lapackpp\lib
Source: demos\*.*; DestDir: {app}\demos;
Source: demos\data\*.*; DestDir: {app}\demos\data; Flags: recursesubdirs
Source: jar\*.*; DestDir: {app}\jar
;
;##############################################################################################################
