// Copyright (C) 1998 John W. Eaton
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
// @deftypefn {Function File} {} randperm (@var{n})
// Return a row vector containing a random permutation of the
// integers from 1 to @var{n}.
// @end deftypefn

// Author: "James R. Van Zandt" <jrv@vanzandt.mv.com>
// Adapted-By: jwe

function retval = randperm(n)

[nargout, nargin] = argn();

if ((nargin == 1) & ((type(n)==1) & (size(n,1)==1) & (size(n,2)==1)) & (floor (n) == n)) then
  if (n >= 0) then
    [junk, retval] = sort(rand (1, n));
  else
    error ("randperm: argument must be non-negative");
  end
else
  error("randperm (n)");
end
endfunction
