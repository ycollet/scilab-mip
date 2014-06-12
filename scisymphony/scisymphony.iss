;##############################################################################################################
; Inno Setup Install script for Toolbox_skeleton
; http://www.jrsoftware.org/isinfo.php
; Yann COLLETTE
; This file is released into the public domain
;##############################################################################################################
; modify this path where is toolbox_skeleton directory
#define BinariesSourcePath "E:\Scilab\toolboxes\scisymphony"
;
#define Scisymphony_version "1.0"
#define CurrentYear "2011"
#define ScisymphonyDirFilename "scisymphony"
;##############################################################################################################
[Setup]
; Debut Données de base à renseigner suivant version
SourceDir={#BinariesSourcePath}
AppName=Scysimphony
AppVerName=Scisymphony version 1.0
DefaultDirName={pf}\{#ScisymphonyDirFilename}
InfoAfterfile=readme.txt
LicenseFile=license.txt
WindowVisible=true
AppPublisher=Home made
BackColorDirection=lefttoright
AppCopyright=Copyright © {#CurrentYear}
Compression=lzma/max
InternalCompressLevel=normal
SolidCompression=true
VersionInfoVersion={#Scisymphony_version}
VersionInfoCompany=Home made
;##############################################################################################################
[Files]
; Add here files that you want to add
Source: loader.sce; DestDir: {app}
Source: *.txt; DestDir: {app}
Source: etc\*.quit; DestDir: {app}\etc
Source: etc\*.start; DestDir: {app}\etc
Source: help\builder_help.sce; DestDir: {app}\help
Source: help\en_US\*.xml; DestDir: {app}\help\en_US
Source: help\en_US\*.sce; DestDir: {app}\help\en_US
Source: jar\*; DestDir: {app}\jar
Source: macros\*; DestDir: {app}\macros
Source: sci_gateway\loader_gateway.sce; DestDir: {app}\sci_gateway
Source: sci_gateway\cpp\loader.sce; DestDir: {app}\sci_gateway\cpp
Source: sci_gateway\cpp\sci_symphony.dll; DestDir: {app}\sci_gateway\cpp
Source: demos\*; DestDir: {app}\demos; Flags: recursesubdirs
;
;##############################################################################################################
