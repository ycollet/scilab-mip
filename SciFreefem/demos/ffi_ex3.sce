// Variational formulation, moving border,Equations are described in Scilab language

// User can enable / disable graphics by setting Graphics variable to %t or %f
if ~exists('Graphics') then
  Graphics = %t;
end

// Deformation of a truss under its own weight

if Graphics then
  xset('wdim',700,700);
  xset('colormap',hotcolormap(256));
end

bord = tlist(['border';'a';'b';'c';'d'],...
              list('x = t; y = 0',0,10,40,1),...
	      list('x = 10; y = t',0,2,20, 1),...
	      list('x = 10 - t; y = 2',0,10,40,1),...
	      list('x = 0; y = t',2,0,20,1));

buildMesh(bord,'th');

// Constants
E      = 21.5;
sigma  = 0.29;
mu     = E/(2*(1+sigma));
lambda = E*sigma/((1+sigma)*(1-2*sigma));
ff_var('gravity',0.05); // Define variable 'gravity' in FreeFem

u = defvar(['u1';'u2']);
v = defvar(['v1';'v2']);
g = defvar(['0';'gravity']);

eps_u = (grad(u)+grad(u)')/2;
eps_v = (grad(v)+grad(v)')/2;

// Previous 2 lines are equivalent to
// eps_u=defvar(['dx(u1)' '(dy(u1)+dx(u2))/2';...
//               '(dy(u1)+dx(u2))/2' 'dy(u2)']);
// eps_v=defvar(['dx(v1)' '(dy(v1)+dx(v2))/2';...
//               '(dy(v1)+dx(v2))/2' 'dy(v2)']);

// Define the problem and solve it

ff_problem(strcat(['varsolve(th,0) bb(u1,v1,u2,v2) with {',...
            pde_varsol('bb','domain',...
            2*mu*(eps_u.*.eps_v)+lambda*tr(eps_u)*tr(eps_v)/2-g*v),...
	    '+ on(b,d)(v1)(u1=0)',...
	    '+ on(b,d)(v2)(u2=0) };']),1);

// Get the result as a Scilab variables
[noeul,trianl,tFunc] = getffResult('u1');

// [Mat,jlow,jhigh,SizeBloc] = getMatrix();	

if Graphics then
  clf();
  xsetech([0,0,1,1]);
  meshvisu();
  //xclick();

  clf();
  xsetech([0,0,1,0.5]);
  xtitle('Displacement versus x');
  resultvisu(noeul,trianl,tFunc);
end

[tFunc] = getffResult('u2');

if Graphics then
  xsetech([0,0.5,1,0.5]);
  xtitle('Displacement versus y');
  resultvisu(noeul,trianl,tFunc);
end

if messagebox('Show final mesh ?','Information','question',['Yes','No'],'modal')==1 then
  // Plot the mesh after deformation
  ff_exec('mesh th1 = movemesh(th, x-u1, y-u2);');
  [noeul,trianl] = getffResult();
  if Graphics then
    clf();
    xsetech([0,0,1,1]);
    meshvisu();
  end
end

ff_end();
