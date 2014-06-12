function []=meshvisu(col,rect)
  // Mesh visualisation
  // uses global variables

  [lhs,rhs] = argn(0);

  if rhs<=0 then col=1;end;
  if rhs<=1 then rect=[min(noeul(:,2)),min(noeul(:,3)),max(noeul(:,2)),max(noeul(:,3))];end;
  if rhs<=2 then iso='1';end;

  plot2d(1,1,[1],'031',' ',rect);
  xset('clipgrf');
  xx = trianl(:,2:4);
  xx = matrix(xx,prod(size(xx)),1);
  x  = noeul(xx,2);
  triang = size(x,'*')/3;
  x = matrix(x,triang,3);
  y = noeul(xx,3);
  y = matrix(y,triang,3);
  x = [x,x(:,1)]';
  y = [y,y(:,1)]';
  xpolys(x,y,col*ones(1,triang));
  xset('clipoff');
endfunction
