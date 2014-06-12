function plot_pipe_graph(g_in,_title,_xlabel,_ylabel,params)
if ~isdef('_title','local')  then _title  = ''; end
if ~isdef('_xlabel','local') then _xlabel = ''; end
if ~isdef('_ylabel','local') then _ylabel = ''; end
if ~isdef('params','local')  then
  params = [];
end

[node_label, err]    = get_param(params,'node_label', %T);
[edge_label, err]    = get_param(params,'edge_label', %T);
[edge_flow, err]     = get_param(params,'edge_flow',  %T);
[node_pressure, err] = get_param(params,'node_pressure', %T);
[node_angle, err]    = get_param(params,'node_angle', 0);
[edge_angle, err]    = get_param(params,'edge_angle', 0);

x_offset = 0.01*(max(g_in('node_x')) - min(g_in('node_x')));
y_offset = 0.01*(max(g_in('node_y')) - min(g_in('node_y')));

drawlater;
for i=1:g_in('node_number')
  plot(g_in('node_x')(i),g_in('node_y')(i),'ro');
  n_label = [];
  if node_label    then n_label = n_label + sprintf('%d',i); end
  if node_pressure then n_label = n_label + sprintf(' (p=%.2f)',g_in.nodes.data.pressure(i)); end
  if ~isempty(n_label) then 
    xstring(g_in('node_x')(i),g_in('node_y')(i),n_label); 
    e = gce();
    e.font_angle = node_angle;
  end 
end

for i=1:length(g_in('tail'))
  x_start = g_in('node_x')(g_in('tail')(i));
  x_end   = g_in('node_x')(g_in('head')(i));
  y_start = g_in('node_y')(g_in('tail')(i));
  y_end   = g_in('node_y')(g_in('head')(i));
  plot([x_start x_end],[y_start y_end],'k-');
  e_label = [];
  if edge_label then e_label = e_label + sprintf('[%d %d]',g_in('tail')(i),g_in('head')(i)); end
  if edge_flow  then e_label = e_label + sprintf(' (f=%.2f)',g_in.edges.data.flow(i)); end
  if ~isempty(e_label) then 
    xstring((x_start+x_end)/2,(y_start+y_end)/2,e_label); 
    e = gce();
    e.font_angle = edge_angle;
  end
end

a=gca();
a.axes_visible = ['off','off','off'];
a.box = 'off';
bounds  = a.data_bounds;
delta_x = 0.025*(max(g_in('node_x')) - min(g_in('node_x')));
delta_y = 0.025*(max(g_in('node_y')) - min(g_in('node_y')));
a.data_bounds = a.data_bounds + [-delta_x,-delta_y;delta_x,delta_y];
xtitle(_title,_xlabel,_ylabel);
drawnow;
endfunction

