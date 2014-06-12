//////////////////////////////////////////////////
// plot_lolimot_part                            //
//                                              //
// Plot the partitions of a lolimot model in 2D //
//////////////////////////////////////////////////

function plot_lolimot_part(lolModel,lol_title)
if typeof(lolModel)~='lolimot' then
  error('plot_lolimot_part: parameter must be of type ''lolimot''\n');
end
if lolModel('nbdim')~=2 then
  return
end

listofrects = [];

for i=1:size(lolModel('listofcutinf'),1)
  // (x,y): point in the upper left corner
  x = lolModel('listofcutinf')(i,1);
  y = lolModel('listofcutplus')(i,2);
  w = lolModel('listofcutplus')(i,1) - lolModel('listofcutinf')(i,1);
  h = lolModel('listofcutplus')(i,2) - lolModel('listofcutinf')(i,2);
  listofrects = [listofrects;[x,y,w,h]];
end

xrects(listofrects');

h = gcf();
h.children.data_bounds = [min(lolModel('listofcutinf'),'r'); max(lolModel('listofcutplus'),'r')];

if isdef('lol_title','local') then
  xtitle(lol_title,'x','y');
else
  xtitle('Cuts','x','y');
end
endfunction
