SciCoinOR Toolbox

This toolbox allows you to use the (mixed integer) linear programming tools and the nonlinear constrained optimization tools under scilab.

With this installation, a some .mps files are shipped. These files allow to test the tools on linear programs and mixed integer linear programs. A MPS file is a format dedicated to handle mixed integer linear programs.

For Windows:

To use this toolbox under scilab: go into the scicoinor toolbox installation directory and do a:

exec loader.sce;

For Linux:

unpack the source code directory, read the wiki page on 

http://code.google.com/p/scilab-mip/wiki/CompilationCoinORLinux

and install GLPK, LPSOLVE and some of the CoinOR tools.
Don't forget to change some paths in sci_gateway/cpp/builder_gateway_cpp.sce file.
Once this is done, you can build the toolbox:

exec builder.sce;

And now, you can load the toolbox:

exec loader.sce;


Once the toolbox has been loaded, close and reopen the help browser to access to the documentation of the functions.
You have also access to a demo via the ?->Demonstrations of Scilab button from the GUI.

In the demos/data,  .mps files are shipped. A description of these files is available in the demos/listoffiles.html file.
A lot more .mps data files are available separately via the other installer.

Yann COLLETTE (ycollette dot nospam at free dot fr)
