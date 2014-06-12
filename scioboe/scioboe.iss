;##############################################################################################################
; Inno Setup Install script for SciOBOE Module
; http://www.jrsoftware.org/isinfo.php
; Yann COLLETTE
; This file is released into the public domain
;##############################################################################################################
; modify this path where is scioboe directory
#define BinariesSourcePath "Z:\toolboxes_google\scioboe"
;
#define SCIOBOE_Module_version "1.0"
#define CurrentYear "2010"
#define SCIOBOE_ModuleDirFilename "scioboe-1.0"
;##############################################################################################################
[Setup]
; Debut Données de base à renseigner suivant version
SourceDir={#BinariesSourcePath}
AppName=SciOBOE Module
AppVerName=SciOBOE Module version 1.0
DefaultDirName={pf}/{#SCIOBOE_ModuleDirFilename}
InfoAfterfile=readme.txt
LicenseFile=license.txt
WindowVisible=true
AppPublisher=Home
BackColorDirection=lefttoright
AppCopyright=Copyright © {#CurrentYear}
Compression=lzma/max
InternalCompressLevel=normal
SolidCompression=true
VersionInfoVersion={#SCIOBOE_Module_version}
VersionInfoCompany=Home
;##############################################################################################################
[Files]
; Add here files that you want to add
Source: loader.sce; DestDir: {app}
Source: license.txt; DestDir: {app}
Source: readme.txt; DestDir: {app}
Source: etc\scioboe.quit; DestDir: {app}\etc
Source: etc\scioboe.start; DestDir: {app}\etc
Source: sci_gateway\loader_gateway.sce; DestDir: {app}\sci_gateway
Source: sci_gateway\cpp\loader.sce; DestDir: {app}\sci_gateway\cpp
Source: sci_gateway\cpp\*.dll; DestDir: {app}\sci_gateway\cpp
Source: demos\*.*; DestDir: {app}\demos;
Source: jar\*.*; DestDir: {app}\jar
Source: help\*.*; DestDir: {app}\help Flags: recursesubdirs
;
;##############################################################################################################
