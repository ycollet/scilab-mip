;##############################################################################################################
; Inno Setup Install script for SciNLopt Module
; http://www.jrsoftware.org/isinfo.php
; Yann COLLETTE
; This file is released into the public domain
;##############################################################################################################
; modify this path where is scinlopt directory
#define BinariesSourcePath "E:\Scilab\toolboxes\scinlopt"
;
#define SCINLOPT_Module_version "1.0"
#define CurrentYear "2011"
#define SCINLOPT_ModuleDirFilename "scinlopt-1.0"
;##############################################################################################################
[Setup]
; Debut Données de base à renseigner suivant version
SourceDir={#BinariesSourcePath}
AppName=SciNLOpt Module
AppVerName=SciNLOpt Module version 1.0
DefaultDirName={pf}/{#SCINLOPT_ModuleDirFilename}
InfoAfterfile=readme.txt
LicenseFile=license.txt
WindowVisible=true
AppPublisher=Home
BackColorDirection=lefttoright
AppCopyright=Copyright © {#CurrentYear}
Compression=lzma/max
InternalCompressLevel=normal
SolidCompression=true
VersionInfoVersion={#SCINLOPT_Module_version}
VersionInfoCompany=Home
;##############################################################################################################
[Files]
; Add here files that you want to add
Source: loader.sce; DestDir: {app}
Source: license.txt; DestDir: {app}
Source: readme.txt; DestDir: {app}
Source: etc\scinlopt.quit; DestDir: {app}\etc
Source: etc\scinlopt.start; DestDir: {app}\etc
Source: sci_gateway\loader_gateway.sce; DestDir: {app}\sci_gateway
Source: sci_gateway\c\loader.sce; DestDir: {app}\sci_gateway\c
Source: sci_gateway\c\*.dll; DestDir: {app}\sci_gateway\c
Source: thirdparty\win\*.dll; DestDir: {app}\thirdparty\win
Source: demos\*.*; DestDir: {app}\demos;
Source: jar\*.*; DestDir: {app}\jar
;
;##############################################################################################################
