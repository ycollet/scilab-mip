function plot_tsp(TSPTownList, TSPVar, TSPComment)
if (size(TSPVar,1)==1) then
  TSPVar = TSPVar';
end
maxX = max(TSPTownList(:,2));
minX = min(TSPTownList(:,2));
maxY = max(TSPTownList(:,3));
minY = min(TSPTownList(:,3));
scf();
a = gca();
a.data_bounds = [minX minY;maxX maxY];
drawlater;
for i=1:size(TSPVar,1)-1
  plot([TSPTownList(TSPVar(i),2);TSPTownList(TSPVar(i+1),2)],[TSPTownList(TSPVar(i),3);TSPTownList(TSPVar(i+1),3)], 'k+-');
end
plot([TSPTownList(TSPVar(1),2);TSPTownList(TSPVar($),2)],[TSPTownList(TSPVar(1),3);TSPTownList(TSPVar($),3)], 'k+-');
xtitle(TSPComment,'X','Y');
drawnow;
endfunction
