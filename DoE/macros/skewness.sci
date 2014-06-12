// Copyright (C) 1996, 1997 John W. Eaton
//
// This file is part of Octave.
//
// Octave is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// Octave is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Octave; see the file COPYING.  If not, write to the Free
// Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.

// -*- texinfo -*-
// @deftypefn {Function File} {} skewness (@var{x}, @var{dim})
// If @var{x} is a vector of length @math{n}, return the skewness
// @iftex
// @tex
// $$
// {\rm skewness} (x) = {1\over N \sigma(x)^3} \sum_{i=1}^N (x_i-\bar{x})^3
// $$
// @end tex
// @end iftex
// @ifinfo
//
// @example
// skewness (x) = N^(-1) std(x)^(-3) sum ((x - mean(x)).^3)
// @end example
// @end ifinfo
//
// @noindent
// of @var{x}.  If @var{x} is a matrix, return the skewness along the
// first non-singleton dimension of the matrix. If the optional
// @var{dim} argument is given, operate along this dimension.
// @end deftypefn

// Author: KH <Kurt.Hornik@wu-wien.ac.at>
// Created: 29 July 1994
// Adapted-By: jwe

function retval = skewness (x)

[nargout, nargin] = argn();

if (nargin ~= 1)
  error ("skewness (x)");
end

if (size(x,1)==1) then
  x = x';
end

sz = size (x,1);

c = size(x,1);
x = x - ones(sz, 1) * mean(x);
retval = 0;
s = stdev (x);
ind = find (s > 0);
x = sum (x .^ 3);
retval = x ./ (c * s ^ 3);
  
endfunction
