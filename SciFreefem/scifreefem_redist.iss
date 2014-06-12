;##############################################################################################################
; Inno Setup Install script for Toolbox_skeleton
; http://www.jrsoftware.org/isinfo.php
; Y. Collette
; This file is released into the public domain
;##############################################################################################################
; modify this path where is toolbox_skeleton directory
#define BinariesSourcePath "E:\Scilab\toolboxes\SciFreeFEM"
;
#define SCIFREEFEM_Toolbox_version "1.2"
#define CurrentYear "2012"
#define SCIFREEFEM_ToolboxDirFilename "scifreefem-1.2"
;##############################################################################################################
[Setup]
; Debut Données de base à renseigner suivant version
SourceDir={#BinariesSourcePath}
AppName=SciFreeFEM Toolbox
AppVerName=SciFreeFEM Toolbox version 1.2
DefaultDirName={pf}/{#SCIFREEFEM_ToolboxDirFilename}
InfoAfterfile=readme.txt
LicenseFile=license.txt	
WindowVisible=true
AppPublisher=Home
BackColorDirection=lefttoright
AppCopyright=Copyright © {#CurrentYear}
Compression=lzma/max
InternalCompressLevel=normal
SolidCompression=true
VersionInfoVersion={#SCIFREEFEM_Toolbox_version}
VersionInfoCompany=Home
;##############################################################################################################
[Files]
; Add here files that you want to add
Source: loader.sce; DestDir: {app}
Source: license.txt; DestDir: {app}
Source: changelog.txt; DestDir: {app}
Source: readme.txt; DestDir: {app}
Source: etc\SciFreeFEM.quit; DestDir: {app}\etc
Source: etc\SciFreeFEM.start; DestDir: {app}\etc
Source: macros\buildmacros.sce; DestDir: {app}\macros
Source: macros\lib; DestDir: {app}\macros
Source: macros\names; DestDir: {app}\macros
Source: macros\*.sci; DestDir: {app}\macros
Source: macros\*.bin; DestDir: {app}\macros
Source: sci_gateway\loader_gateway.sce; DestDir: {app}\sci_gateway
Source: sci_gateway\c\loader.sce; DestDir: {app}\sci_gateway\c
Source: sci_gateway\c\freefem_c.dll; DestDir: {app}\sci_gateway\c
Source: src\c\libdelete_lb_obj_.dll; DestDir: {app}\src\c
Source: src\c\loader.sce; DestDir: {app}\src\c
Source: help\pdf\*.*; DestDir: {app}\help\pdf;
Source: demos\*.*; DestDir: {app}\demos;
Source: demos\EDP\*.*; DestDir: {app}\demos\EDP;
Source: jar\*.*; DestDir: {app}\jar
;
;##############################################################################################################
