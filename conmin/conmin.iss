;##############################################################################################################
; Inno Setup Install script for Toolbox_skeleton
; http://www.jrsoftware.org/isinfo.php
; Allan CORNET
; This file is released into the public domain
;##############################################################################################################
; modify this path where is toolbox_skeleton directory
#define BinariesSourcePath "E:\Scilab\toolboxes\conmin"
;
#define Conmin_version "2.1"
#define CurrentYear "2011"
#define ConminDirFilename "conmin"
;##############################################################################################################
[Setup]
; Debut Données de base à renseigner suivant version
SourceDir={#BinariesSourcePath}
AppName=Conmin
AppVerName=Conmin version 2.1
DefaultDirName={pf}\{#ConminDirFilename}
InfoAfterfile=readme.txt
LicenseFile=license.txt
WindowVisible=true
AppPublisher=Home made
BackColorDirection=lefttoright
AppCopyright=Copyright © {#CurrentYear}
Compression=lzma/max
InternalCompressLevel=normal
SolidCompression=true
VersionInfoVersion={#Conmin_version}
VersionInfoCompany=Home made
;##############################################################################################################
[Files]
; Add here files that you want to add
Source: loader.sce; DestDir: {app}
Source: *.txt; DestDir: {app}
Source: etc\conmin.quit; DestDir: {app}\etc
Source: etc\conmin.start; DestDir: {app}\etc
Source: help\builder_help.sce; DestDir: {app}\help
Source: help\en_US\*.xml; DestDir: {app}\help\en_US
Source: help\en_US\*.sce; DestDir: {app}\help\en_US
Source: jar\*; DestDir: {app}\jar
Source: macros\*; DestDir: {app}\macros
Source: sci_gateway\loader_gateway.sce; DestDir: {app}\sci_gateway
Source: sci_gateway\c\conmin_c.dll; DestDir: {app}\sci_gateway\c
Source: sci_gateway\c\loader.sce; DestDir: {app}\sci_gateway\c
Source: src\fortran\libget_ct.dll; DestDir: {app}\src\fortran
Source: src\fortran\loader.sce; DestDir: {app}\src\fortran
Source: demos\*; DestDir: {app}\demos; Flags: recursesubdirs
;
;##############################################################################################################
