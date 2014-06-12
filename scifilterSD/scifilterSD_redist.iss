;##############################################################################################################
; Inno Setup Install script for scifilterSD
; http://www.jrsoftware.org/isinfo.php
; Yann COLLETTE
; This file is released under the 3-clause BSD license. See COPYING-BSD.
;##############################################################################################################
; modify this path where is toolbox_skeleton directory
#define BinariesSourcePath "E:\Scilab\Toolboxes\scifilterSD"
;
#define scifilterSD_version "0.1"
#define CurrentYear "2012"
#define scifilterSDDirFilename "scifilterSD"
;##############################################################################################################
[Setup]
; Debut Données de base à renseigner suivant version
SourceDir={#BinariesSourcePath}
AppName=scifilterSD
AppVerName=scifilterSD 0.1
DefaultDirName={pf}\{#scifilterSDDirFilename}
InfoAfterfile=readme.txt
LicenseFile=license.txt
WindowVisible=true
AppPublisher=Home
BackColorDirection=lefttoright
AppCopyright=Copyright © {#CurrentYear}
Compression=lzma/max
InternalCompressLevel=normal
SolidCompression=true
VersionInfoVersion={#scifilterSD_version}
VersionInfoCompany=Home
;##############################################################################################################
[Files]
; Add here files that you want to add
Source: loader.sce; DestDir: {app}
Source: etc\scifilterSD.quit; DestDir: {app}\etc
Source: etc\scifilterSD.start; DestDir: {app}\etc
;Source: macros\buildmacros.sce; DestDir: {app}\macros
;Source: macros\lib; DestDir: {app}\macros
;Source: macros\names; DestDir: {app}\macros
;Source: macros\*.sci; DestDir: {app}\macros
;Source: macros\*.bin; DestDir: {app}\macros
Source: sci_gateway\loader_gateway.sce; DestDir: {app}\sci_gateway
Source: sci_gateway\dense\loader.sce; DestDir: {app}\sci_gateway\dense
Source: sci_gateway\dense\filtersd_dense_c.dll; DestDir: {app}\sci_gateway\dense
Source: sci_gateway\sparse\loader.sce; DestDir: {app}\sci_gateway\sparse
Source: sci_gateway\sparse\filtersd_sparse_c.dll; DestDir: {app}\sci_gateway\sparse
;Source: tests\*.*; DestDir: {app}\tests; Flags: recursesubdirs
Source: demos\*.*; DestDir: {app}\demos; Flags: recursesubdirs
Source: jar\*.*; DestDir: {app}\jar; Flags: recursesubdirs
Source: help\*.*; DestDir: {app}\help; Flags: recursesubdirs
;
;##############################################################################################################
