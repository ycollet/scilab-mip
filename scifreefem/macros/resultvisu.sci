function [] = resultvisu(node,triangle,func,varargin)
  [lhs,rhs] = argn(0);
  if (rhs < 3) then error('incorrect number of arguments'); end;
  fec(node(:,2),node(:,3),triangle,func,varargin(:));
endfunction
