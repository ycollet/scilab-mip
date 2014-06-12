// Variational formulation, Equations are described in Scilab language

// User can enable / disable graphics by setting Graphics variable to %t or %f
if ~exists('Graphics') then
  Graphics = %t;
end

// Definition of the border
deff('[x,y] = f1(t)','x = 0; y = 1 - t');

bord = tlist(['border';'a';'b';'c';'d';'e';'f'],...
             list('x = t; y = 0',0,1,6,1),...
	     list('x = 1; y = t',0,0.5,4, 1),...
	     list('x = 1 - t; y = 0.5',0,0.5,4,1),...
	     list('x = 0.5; y = t',0.5,1,4,1),...
	     list('x = 1 - t; y = 1',0.5,1,4,1),...
	     list(valf(f1, 0:1/6:1), 1));

buildMesh(bord,'th');

funcprot(0); // Avoid a warning for the redefinition of f1(t)
u = defvar('u'); // Set the unknown
f = defvar('1'); // Set the right hand side of the equation

if Graphics then
  xset('colormap',hotcolormap(254));
  xset('wdim',700,700);
end

for i=1:3;
  // Define and solve the problem

  ff_problem(strcat(['solve(u)  {',...
	             pde_sol(u,div(grad(u)),f),...
                     'on(a,b,c,d,e,f) u = 0; };']));

  // Get the result as a Scilab variables
  [noeul,trianl,tFunc] = getffResult('u');

  if Graphics then
    // Show the mesh
    clf(); xsetech([0,0,1,0.5]);
    meshvisu();
    
    // Show the result
    xsetech([0,0.5,1,0.5]);
    resultvisu(noeul,trianl,tFunc);
  end

  if messagebox('Update mesh','Information','question',['Yes','No'],'modal')==2 then break,end

  // Update mesh
  deff('[x,y] = f1(t)','x = (0.5/(4-i))*cos(t); y = 0.5+0.5*sin(t)');
  bord('f')(1) = valf(f1, %pi/2:%pi/(5*i):3*%pi/2);
  UpdateMesh(bord,'th'); 
end

ff_end(); // Destroy Fem interpretor
