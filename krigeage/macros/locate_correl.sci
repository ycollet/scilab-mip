function list_correl = locate_correl()
 list_aux = who('local'); 
 list_aux = list_aux(grep(list_aux,'correl_'));
 if ~isempty(list_aux) then
   list_aux = list_aux(find(part(list_aux,1)=='c'));
   list_correl = stripblanks(part(list_aux,8:max(length(list_aux))));
 else
   list_correl = [];
 end
 list_correl = [list_correl;'cubic';'exp';'expg';'gauss';'lin';'materm';'sinus';'spherical';'spline'];
 list_correl = unique(list_correl);
endfunction
