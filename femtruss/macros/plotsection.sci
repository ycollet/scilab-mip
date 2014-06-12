function plotsection(T, P, S)
//---------------------------------------------------
// plotsection(T, P, S, NodeLabels, EltLabels)
// 
//  T the connectivity between points
//  P the points
//  S sections of the bars (same number of lines as T)

Net = size(T,1);
_3D_problem = (size(P,2)==3);

h = gcf();
h.color_map = graycolormap(256);

drawlater;

delta_x = max(P(:,1))-min(P(:,1));
delta_y = max(P(:,2))-min(P(:,2));
x_max = max(P(:,1));
x_min = min(P(:,1));
y_max = max(P(:,2));
y_min = min(P(:,2));

if _3D_problem then
  delta_z = max(P(:,2))-min(P(:,2));
  z_max = max(P(:,2));
  z_min = min(P(:,2));

  f.data_bounds = [x_min-0.1*delta_x, y_min-0.1*delta_y, z_min-0.1*delta_z; x_max+0.1*delta_x, y_max+0.1*delta_y, z_max+0.1*delta_z];
else
  f.data_bounds = [x_min-0.1*delta_x, y_min-0.1*delta_y; x_max+0.1*delta_x, y_max+0.1*delta_y];
end

_min_section = min(S);
_max_section = max(S);
nb_color_max = 200;
nb_color_min = 20;

colorbar(_min_section,_max_section,[nb_color_min, nb_color_max]);

for ie = 1:Net
  Section = S(ie);
  _color = floor((nb_color_max - nb_color_min)*(_max_section - Section)/max(_max_section - _min_section,%eps)) + nb_color_min;
  index_begin = T(ie,1);
  index_end   = T(ie,2);
  plot2d([P(index_begin,1) P(index_end,1)],[P(index_begin,2) P(index_end,2)],color(ceil(255*h.color_map(_color,1)),ceil(255*h.color_map(_color,2)),ceil(255*h.color_map(_color,3))));
  e = gce();
  if _3D_problem then
    e.children.data = [P(index_begin,1) P(index_end,2) P(index_end,3)],[P(index_begin,1) P(index_end,2) P(index_end,3)];
  end
end 

if _3D_problem then
  h.children(2).view = '3d';
end
h.children(2).axes_visible = ['off','off','off'];

drawnow;
endfunction
