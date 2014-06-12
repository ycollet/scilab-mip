;##############################################################################################################
; Inno Setup Install script for MySQL toolbox
; http://www.jrsoftware.org/isinfo.php
; Yann COLLETTE
; This file is released into the public domain
;##############################################################################################################
; modify this path where is toolbox_skeleton directory
#define BinariesSourcePath "Z:\toolboxes_google\mysql"
;
#define Toolbox_mysql_version "2.0"
#define CurrentYear "2010"
#define Toolbox_mysqlDirFilename "mysql"
;##############################################################################################################
[Setup]
; Debut Données de base à renseigner suivant version
SourceDir={#BinariesSourcePath}
AppName=MySQL
AppVerName=MySQL toolbox version 2.0
DefaultDirName={pf}\{#Toolbox_mysqlDirFilename}
InfoAfterfile=readme.txt
LicenseFile=license.txt
WindowVisible=true
AppPublisher=Home
BackColorDirection=lefttoright
AppCopyright=Copyright © {#CurrentYear}
Compression=lzma/max
InternalCompressLevel=normal
SolidCompression=true
VersionInfoVersion={#Toolbox_mysql_version}
VersionInfoCompany=Home
;##############################################################################################################
[Files]
; Add here files that you want to add
Source: loader.sce; DestDir: {app}
Source: etc\mysql.quit; DestDir: {app}\etc
Source: etc\mysql.start; DestDir: {app}\etc
;Source: macros\buildmacros.sce; DestDir: {app}\macros
;Source: macros\lib; DestDir: {app}\macros
;Source: macros\names; DestDir: {app}\macros
;Source: macros\*.sci; DestDir: {app}\macros
;Source: macros\*.bin; DestDir: {app}\macros
Source: help\*.*; DestDir: {app}\help
Source: jar\*.*; DestDir: {app}\jar
Source: sci_gateway\loader_gateway.sce; DestDir: {app}\sci_gateway
Source: sci_gateway\cpp\loader.sce; DestDir: {app}\sci_gateway\cpp
Source: sci_gateway\cpp\mysql_cpp.dll; DestDir: {app}\sci_gateway\cpp
;Source: tests\*.*; DestDir: {app}\tests; Flags: recursesubdirs
Source: demos\*.*; DestDir: {app}\locales; Flags: recursesubdirs
;
;##############################################################################################################
