function []=nvisu(rect)
  // Visualisation des noeuds

  [lhs,rhs] = argn(0);
  if rhs==0 then rect=[min(noeul(:,2)),min(noeul(:,3)),max(noeul(:,2)),max(noeul(:,3))]; end;
  plot2d(1,1,[1],'031',' ',rect);
  xset('clipgrf');
  bords = noeul(find(noeul(:,4)>0),:);
  [no,ign] = size(bords);
  for i=1:no
    xstring(bords(i,2),bords(i,3),string(bords(i,4)));
  end
  xset('clipoff');
endfunction
