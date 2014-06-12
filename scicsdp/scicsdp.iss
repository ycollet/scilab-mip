;##############################################################################################################
; Inno Setup Install script for SciCSDP Module
; http://www.jrsoftware.org/isinfo.php
; Yann COLLETTE
; This file is released into the public domain
;##############################################################################################################
; modify this path where is sciglpk directory
#define BinariesSourcePath "E:\Scilab\toolboxes\scicsdp"
;
#define SCICSDP_Module_version "0.1"
#define CurrentYear "2012"
#define SCICSDP_ModuleDirFilename "scicsdp-0.1"
;##############################################################################################################
[Setup]
; Debut Données de base à renseigner suivant version
SourceDir={#BinariesSourcePath}
AppName=SciCSDP Module
AppVerName=SciCSDP Module version 0.1
DefaultDirName={pf}/{#SCICSDP_ModuleDirFilename}
InfoAfterfile=readme.txt
LicenseFile=license.txt
WindowVisible=true
AppPublisher=Home
BackColorDirection=lefttoright
AppCopyright=Copyright © {#CurrentYear}
Compression=lzma/max
InternalCompressLevel=normal
SolidCompression=true
VersionInfoVersion={#SCICSDP_Module_version}
VersionInfoCompany=Home
;##############################################################################################################
[Files]
; Add here files that you want to add
Source: loader.sce; DestDir: {app}
;Source: unloader.sce; DestDir: {app}
Source: license.txt; DestDir: {app}
Source: readme.txt; DestDir: {app}
Source: etc\scicsdp.quit; DestDir: {app}\etc
Source: etc\scicsdp.start; DestDir: {app}\etc
Source: sci_gateway\loader_gateway.sce; DestDir: {app}\sci_gateway
Source: sci_gateway\cpp\loader.sce; DestDir: {app}\sci_gateway\cpp
Source: sci_gateway\cpp\*.dll; DestDir: {app}\sci_gateway\cpp
Source: src\c\loader.sce; DestDir: {app}\src\c
Source: src\c\libsdp.dll; DestDir: {app}\src\c
Source: demos\*.*; DestDir: {app}\demos;
Source: demos\data\sdpa\*.*; DestDir: {app}\demos\data\sdpa;
Source: jar\*.*; DestDir: {app}\jar
;
;##############################################################################################################
