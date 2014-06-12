;##############################################################################################################
; Inno Setup Install script for Scicoinor Data Set CORAL
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
AppName=SciCoinOR Toolbox Data Set CORAL
AppVerName=SciCoinOR Toolbox version 1.2
DefaultDirName={pf}/{#SCICOINOR_ToolboxDirFilename}
InfoAfterfile=readme_data2.txt
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
Source: demos\data2\Coral\*.*; DestDir: {app}\demos\data2\Coral; Flags: recursesubdirs
;
;##############################################################################################################
