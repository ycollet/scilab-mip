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
// @deftypefn {Function File} {} kurtosis (@var{x}, @var{dim})
// If @var{x} is a vector of length @math{N}, return the kurtosis
// @iftex
// @tex
// $$
//  {\rm kurtosis} (x) = {1\over N \sigma(x)^4} \sum_{i=1}^N (x_i-\bar{x})^4 - 3
// $$
// @end tex
// @end iftex
// @ifinfo
//
// @example
// kurtosis (x) = N^(-1) std(x)^(-4) sum ((x - mean(x)).^4) - 3
// @end example
// @end ifinfo
//
// @noindent
// of @var{x}.  If @var{x} is a matrix, return the kurtosis over the
// first non-singleton dimension. The optional argument @var{dim}
// can be given to force the kurtosis to be given over that 
// dimension.
// @end deftypefn

// Author: KH <Kurt.Hornik@wu-wien.ac.at>
// Created: 29 July 1994
// Adapted-By: jwe

function retval = kurtosis (x)

[nargout, nargin] = argn();

if (nargin ~= 1) then
  error("kurtosis (x)");
end

if (size(x,1)==1) then
  x = x';
end

sz = size (x,1);

c = sz;
// x = x - repmat (mean (x, dim), idx);
x = x - ones(sz, 1) * mean (x);
retval = 0;
s = stdev(x);
x = sum(x.^4);
ind = find (s > 0);
retval = x ./ (c * s ^ 4) - 3;

endfunction

