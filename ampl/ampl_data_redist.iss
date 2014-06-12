;##############################################################################################################
; Inno Setup Install script for Toolbox_skeleton
; http://www.jrsoftware.org/isinfo.php
; Allan CORNET
; This file is released into the public domain
;##############################################################################################################
; modify this path where is toolbox_skeleton directory
#define BinariesSourcePath "Z:\toolboxes_google\ampl"
;
#define AMPL_Toolbox_version "1.0"
#define CurrentYear "2009"
#define AMPL_ToolboxDirFilename "scilab-ampl-1.0"
;##############################################################################################################
[Setup]
; Debut Données de base à renseigner suivant version
SourceDir={#BinariesSourcePath}
AppName=AMPL Toolbox
AppVerName=AMPL Toolbox version 1.0
DefaultDirName={pf}/{#AMPL_ToolboxDirFilename}
InfoAfterfile=readme.txt
LicenseFile=license.txt
WindowVisible=true
AppPublisher=Home
BackColorDirection=lefttoright
AppCopyright=Copyright © {#CurrentYear}
Compression=lzma/max
InternalCompressLevel=normal
SolidCompression=true
VersionInfoVersion={#AMPL_Toolbox_version}
VersionInfoCompany=Home
;##############################################################################################################
[Files]
; Add here files that you want to add
Source: demos\data2\*.*; DestDir: {app}\demos\data2; Flags: recursesubdirs
;
;##############################################################################################################

