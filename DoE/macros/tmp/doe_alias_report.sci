function doe_alias_report(alias_list)
for i=1:size(alias_list)
  printf('Alias generator : %s - monomials aliased : ', alias_list(i)(1));
  if (size(alias_list(i))==1) then
    printf('no monomials.\n');
  elseif (size(alias_list(i))==2) then
    printf('%s.\n', alias_list(i)($));
  else
    for j=2:size(alias_list(i))-1
      printf('%s, ', alias_list(i)(j));
    end
    printf('%s.\n', alias_list(i)($));
  end
end
endfunction
