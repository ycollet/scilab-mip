function dr = derive(d,ch)
b = strindex(ch,'dx');
c = strindex(ch,'dy');
a = 0;
res = ch;

if (b <> []) then
  res = strsubst(res,'dx(','d'+d+'x(');
  a   = 1;
end

if (c <> []) then
  res = strsubst(res,'dy(','d'+d+'y(');
  a   = 1;
end

if (a == 0) then
  res = 'd'+d+'('+res+')';
end

dr = res;
endfunction
