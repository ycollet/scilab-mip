function [TSP_Opt] = read_tsp_opt(Filename)
[fd, err] = mopen(Filename,'r');
if (err) then
  printf('error while reading file %s\n',Filename);
  abort;
end

TSP_Opt = [];
Line    = [];

while(Line~='TOUR_SECTION')
  Line = mgetl(fd,1);
  Line = stripblanks(Line);
end

Index = 1;
while(Line~='-1')
  Line = mgetl(fd,1);
  Line = stripblanks(Line);
  if Line=='-1' then break; end
  TSP_Opt(Index) = eval(Line);
  Index = Index + 1
end
mclose(fd);
endfunction
