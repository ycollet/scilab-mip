// Standard formulation problem, Equations are described in FreeFeem language.

// User can enable / disable graphics by setting Graphics variable to %t or %f
if ~exists('Graphics') then
  Graphics = %t;
end

if Graphics then
  xset('wdim',700,700);
  xset('colormap',hotcolormap(254));
end

// Definition of the border

bord = tlist(['border';'a';'b';'c';'d';'e';'f'],...
             list('x = t; y = 0',0,1,6,1),...
	     list('x = 1; y = t',0,0.5,4, 1),...
	     list('x = 1 - t; y = 0.5',0,0.5,4,1),...
	     list('x = 0.5; y = t',0.5,1,4,1),...
	     list('x = 1 - t; y = 1',0.5,1,4,1),...
	     list('x = 0; y = 1 - t',0,1,6,1));

buildMesh(bord,'th'); // Mesh building

[noeul,trianl] = getffResult(); // Get mesh description in Scilab

erreur = 0.1;
coef   = 0.1^(1./5.); // Error level will be divided by 10 every 5 iterations.

for i=1:4;
  // Define and solve the problem
  ff_problem('solve(u) {pde(u) laplace(u) = 1; on(a,b,c,d,e,f) u = 0;};',1);
  
  // Get the result as a Scilab variable
  [tFunc] = getffResult('u');

  // Get the matrix of the linear problem
  [Mat,jlow,jhigh,SizeBloc] = getMatrix();

  // Show the mesh
  if Graphics then
    clf(); xsetech([0,0,1,0.5]);
    meshvisu(); 

    // Show the result
    xsetech([0,0.5,1,0.5]);
    resultvisu(noeul,trianl,tFunc); 
  end

  if messagebox('Adapt mesh','Information','question',['Yes','No'],'modal')==2 then break,end

  erreur = erreur*coef; // Adapt error coefficient

  // Adapt mesh

  ff_adaptmesh('mesh th = adaptmesh  (th,u)',verbosity = 5,...
               abserror = 1, nbjacoby = 2, err = erreur, nbvx = 5000, ...
               omega = 1.8 ,ratio = 1.8, nbsmooth = 3, splitpbedge = 1, ...
               maxsubdiv = 5, rescaling = 1);
  if (i < 4) then
    [noeul,trianl] = getffResult(); // Get the new mesh
  end
end

ff_problem('solve(u) {pde(u) laplace(u) = 1; on(a,b,c,d,e,f) u=0;};');

// Get mesh named  'th' and result in Scilab environment
[noeul,trianl,tFunc] = getffResult('th,u');

if Graphics then
  clf();
  xsetech([0,0,0.5,0.5]);
  meshvisu(100); // Show mesh

  xsetech([0.5,0,0.5,0.5]);
  meshvisu(160,[0.4 0.4 0.55 0.55]); 

  xsetech([0,0.5,1,0.5]);
  resultvisu(noeul,trianl,tFunc);
end

ff_end(); // Destroy Fem interpretor
