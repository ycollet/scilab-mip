Description 
=========== 
Scilab Interface to FreeFem
Author
======
Emmanuel Geay (updated by Y. Collette)


Installation
============
To Install this toolbox: (Scilab-5.3 and above)

We Suppose here that <PATH> stands for  the path of the directory
containing this README file.

- On Unix/Linux / Windows systems
     * Administrator
        Has to execute, once and for all, the following instruction 
        within Scilab:
        exec <PATH>/builder.sce
        This operation requires a C++ and C compiler and
        permission to write in
	- <PATH>/macros  to generate *.bin, names and lib files 
        - <PATH>/src/c to compile FreeFem routines
        - <PATH>/sci_gateway/c to compile FreeFem interface

     * User
  	Should execute the following instruction within Scilab:
	exec <PATH>/loader.sce
	before using the toolbox, he  can also put it  in his
        .scilab startup file for automatic loading.

Contents
========
README             : this file
builder.sce        : script for buliding library
src                : directory of C++ routines
                     *.cpp *.h FreeFem library files
sci_gateway        : directory of the C freefemfi.c interface with Scilab
macros             : directory of Scilab functions
     *.sci         : source versions
help               : directory for help.
