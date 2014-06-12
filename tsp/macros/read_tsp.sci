function [TSP_Town_List, TSPSize, TSPComment] = read_tsp(Filename)
[fd, err] = mopen(Filename,'r');
if (err) then
  printf('error while reading file %s\n',Filename);
  abort;
end

Line = mgetl(fd,1);
Line = stripblanks(tokens(Line,':'));

while(Line(1)~='EOF')
  select(Line(1))
  case 'NAME' then
    TSPName = Line(2);
  case 'COMMENT' then
    TSPComment = Line(2);
  case 'TYPE' then
    if (Line(2)~='TSP') then
      error('the instance to read must be of ''TSP'' type');
    end
  case 'EDGE_WEIGHT_TYPE' then
    if (Line(2)~='EUC_2D') then
      error('the edges of this instance to read must be of ''EUC_2D'' type');
    end
  case 'NODE_COORD_SECTION' then
    TspDataType = Line(1);
    if (TspDataType~='NODE_COORD_SECTION') then
      error('this function is only able to read data of the type NODE_COORD_SECTION');
    end
  case 'DIMENSION' then
    TSPSize = eval(Line(2));
  else
    Line = stripblanks(tokens(Line,' '));
    Index = 0;
    while(Line(1)~='EOF')
      Index = Index + 1;
      for i=1:size(Line,1)
        TSP_Town_List(Index,i) = eval(Line(i));
      end
      Line = mgetl(fd,1);
      Line = stripblanks(tokens(Line,' '));
    end // While
  end
  if (Line(1)~='EOF') then
    Line = mgetl(fd,1);
    Line = stripblanks(tokens(Line,':'));
  end
end
mclose(fd);
endfunction
