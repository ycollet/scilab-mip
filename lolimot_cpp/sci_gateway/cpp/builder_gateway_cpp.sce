// ====================================================================
// Allan CORNET
// DIGITEO 2008
// This file is released into the public domain
// ====================================================================

src_path = get_absolute_file_path('builder_gateway_cpp.sce');

FILES = ['Analysis.cpp', ...
         'Validate.cpp', ...
	 'ExportC.cpp','ExportMatlab.cpp','Learning.cpp', ...
	 'LL_Lolimot.cpp','Post.cpp', ...
	 'DefaultParam.cpp','ExportCommonC_Cpp.cpp', ...
	 'LL_Cut.cpp','LL_Mesure.cpp', ...
	 'TrainLolimotStruct.cpp','DisplayData.cpp', ...
	 'ExportCpp.cpp','LL_Dimension.cpp', ...
	 'LL_Partition.cpp','sci_learn_lolimot.cpp'];

TABLE = ['lolimot','sci_lolimot'];

LIBS = [];

LDFLAGS = [];
CFLAGS  = '-I' + src_path;

tbx_build_gateway('lolimot_cpp', TABLE, FILES, src_path, LIBS,LDFLAGS,CFLAGS);

clear tbx_build_gateway;
