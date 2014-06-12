// ====================================================================
// Y. Collette
// Copyright 2012
// This file is released into the public domain
// ====================================================================

if (getos()=='Linux' | getos()=='Darwin' | getos()=='SunOS') then
  CFLAGS = '-fPIC -g -I' + get_absolute_file_path('builder_c.sce') + ' -DNOMAIN  -DFREEFEM -DF77_UNDERSCORE';
else
  CFLAGS = '-I' + get_absolute_file_path('builder_c.sce') + ' -DNOMAIN -DFREEFEM -DF77_UNDERSCORE';
end

UseGraph = %f;

SRC = ['Mesh2.cpp', 'MeshDraw.cpp', 'MeshGeom.cpp', 'MeshQuad.cpp', 'MeshRead.cpp', 'MeshWrite.cpp', 'Meshio.cpp', ...
       'Metric.cpp', 'QuadTree.cpp', 'R2.cpp', 'SetOfE4.cpp', 'BamgFreeFem.cpp', 'analyse.cpp', ...
       'fonction.cpp', 'gestchar.cpp', 'gibbs.cpp', 'ivarsol.cpp', 'cglapl.cpp', 'convec.cpp', ...
       'grid.cpp', 'vect.cpp', 'fem2.cpp', 'list.cpp', 'graph.cpp'];

if UseGraph then
  if (getos()=='Linux' | getos()=='SunOS') then
    SRC = [SRC, 'Xrgraph.cpp'];
  end
  
  if (getos()=='Windows') then
    SRC = [SRC, 'pcrgraph.cpp'];
  end
  
  if (getos()=='Darwin') then
    SRC = [SRC, 'macrgraph.cpp'];
  end
else
  SRC = [SRC, 'txtgraph.cpp'];
end

NAMES = ['delete_lb_lobj_', 'delete_lexp_', 'add_lobj_', 'delete_scilabana_', ...
         'freefem_code_', 'ff_problem_', 'get_ff_result_', 'put_scilab_mesh_', ...
	 'put_scilab_border_1_', 'put_scilab_border_2_', 'build_scilab_mesh_', ...
	 'get_matrix_'];

LIBRARY = 'libfreefem';

tbx_build_src(NAMES, SRC, "c", get_absolute_file_path('builder_c.sce'),"","",CFLAGS);

clear tbx_build_src;
