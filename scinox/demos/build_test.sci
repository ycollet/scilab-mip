function [g,I_source,I_leaves] = build_test(_type,do_plot)
if ~isdef('do_plot','local') then
  do_plot = %F;
end

g      = [];
I_head = [];
I_tail = [];

if _type=='simple' then
  // Node parameters
  // n°                 1   2   3   4   5   6   7
  node_x           = [250 150 350 100 200 300 400];
  node_y           = [300 200 200 100 100 100 100];
  node_z           = [  0  10  10  20  20  20  20];
  consumption_node = [  0   0   0   1   2   3   4];
  node_pressure    = [  0   0   0   0   0   0   0];
  
  // Edge parameters
  // n°         1  2  3  4  5  6
  I_start   = [ 1  1  2  2  3  3];
  I_end     = [ 2  3  4  5  6  7];
  diameter  = [50 50 25 25 25 25]/1000; // Meters
  flow_pipe = [ 0  0  0  0  0  0];
  
  InstanceName    = 'simple';
  I_source        = [1];         // The source nodes
  source_pressure = [10]; // 10 bars
  I_leaves        = [4 5 6 7];   // The leaves nodes
elseif _type=='simple_loop_1' then
  // Node parameters
  // n°                 1   2   3   4   5   6   7
  node_x           = [250 150 350 100 200 300 400];
  node_y           = [300 200 200 100 100 100 100];
  node_z           = [  0  10  10  20  20  20  20];
  consumption_node = [  0   0   0   1   2   3   4];
  node_pressure    = [  0   0   0   0   0   0   0];
  
  // Pipe parameters
  // n°          1  2  3  4  5  6  7
  diameter  = [100 50 50 25 25 25 25]/1000; // Meters
  I_start   = [  1  1  2  2  3  3  2];
  I_end     = [  2  3  4  5  6  7  3];
  flow_pipe = [  0  0  0  0  0  0  0];
  
  InstanceName    = 'simple loop';
  I_source        = [1];
  source_pressure = [10]; // 10 bars
  I_leaves        = [4 5 6 7];
elseif _type=='simple_loop_2' then
  L = 100;
  // Node parameters
  // n°               1   2   3   4   5   6   7   8
  node_x           = [1   1   2   1   2   1   0   0]*L;
  node_y           = [3   2   2   1   1   0   2   1]*L;
  node_z           = [0   0   0   0   0   0   0   0];
  consumption_node = [0   0   0   0   0 100   0   0];
  node_pressure    = [0   0   0   0   0   1   0   0];
  
  // Edge parameters
  // n°          1  2  3  4  5  6  7  8  9
  I_start   = [  1  2  3  5  2  4  2  7  8];
  I_end     = [  2  3  5  4  4  6  7  8  4];
  diameter  = [100 50 50 50 50 50 50 50 50]/1000; // Meters
  flow_pipe = [  0  0  0  0  0  0  0  0  0];
  
  InstanceName    = 'simple loop 2';
  I_source        = [1];
  source_pressure = [10]; // 10 bars
  I_leaves        = [6];
elseif _type=='simple_loop_3' then
  L = 100;
  // Node parameters
  // n°               1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
  node_x           = [0   0   1   2   3   0   1   2   3   0   1   2   3   0   1   2   3   0   3   2.5]*L;
  node_y           = [5   4   4   4   4   3   3   3   3   2   2   2   2   1   1   1   1   0   1.5 1]*L;
  node_z           = [0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0];
  consumption_node = [0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 100   0   0];
  node_pressure    = [0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0];
  
  // Edge parameters
  // n°           1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28
  I_start     = [ 1  2  3  4  2  3  4  5  6  7  8  6  7  8  9 10 11 12 10 11 12 13 19 14 15 16 20 14];
  I_end       = [ 2  3  4  5  6  7  8  9  7  8  9 10 11 12 13 11 12 13 14 15 16 19 17 15 16 20 17 18];
  diameter    = [50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50]/1000; // Meters
  flow_pipe   = [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0];
  
  InstanceName    = 'simple loop 3';
  I_source        = [1];
  source_pressure = [10]; // 10 bars
  I_leaves        = [18];
elseif _type=='test_hc_1' then
  // Test from:
  // Hardy Cross, "Analysis of flow in networks of conduits or conductors", 
  // University of Illinois Engineering Experiment Station, Bulletin n° 286, november 1936

  L = 100;
  // Node parameters
  // n°               1   2   3   4   5   6
  node_x           = [0   0   1   0   1   0]*L;
  node_y           = [3   2   2   1   1   0]*L;
  node_z           = [0   0   0   0   0   0];
  consumption_node = [0   0   0   0   0 100];
  node_pressure    = [0   0   0   0   0   0];
  
  // Edge parameters
  // n°            1   2  3  4  5  6
  I_start     = [  1   2  3  5  2  4];
  I_end       = [  2   3  5  4  4  6];
  diameter    = [100  50 50 50 50 50]/1000; // Meters
  flow_pipe   = [  0   0  0  0  0  0];
  
  InstanceName    = 'Test Hardy Cross 1';
  I_source        = [1];
  source_pressure = [10]; // 10 bars
  I_leaves        = [6];
elseif _type=='test_hc_1_bis' then
  // Test from:
  // Hardy Cross, "Analysis of flow in networks of conduits or conductors", 
  // University of Illinois Engineering Experiment Station, Bulletin n° 286, november 1936

  L = 100;
  // Node parameters
  // n°               1   2   3   4   5   6
  node_x           = [0   0   1   0   1   0]*L;
  node_y           = [3   2   2   1   1   0]*L;
  node_z           = [0   0   0   0   0   0];
  consumption_node = [0   0   0   0   0 100];
  node_pressure    = [0   0   0   0   0   0];
  
  // Edge parameters
  // n°            1   2  3  4  5  6  7
  I_start     = [  1   2  3  5  2  4  5];
  I_end       = [  2   3  5  4  4  6  6];
  diameter    = [100  50 50 50 50 50 50]/1000; // Meters
  flow_pipe   = [  0   0  0  0  0  0  0];
  
  InstanceName    = 'Test Hardy Cross 1 bis';
  I_source        = [1];
  source_pressure = [10]; // 10 bars
  I_leaves        = [6];
elseif _type=='test_hc_2' then
  // Test from:
  // Hardy Cross, "Analysis of flow in networks of conduits or conductors", 
  // University of Illinois Engineering Experiment Station, Bulletin n° 286, november 1936

  L = 100;
  // Node parameters
  // n°               1   2   3   4   5   6   7   8   9   10  11
  node_x           = [1   0   1   2   0   1   2   0   1   2   1]*L;
  node_y           = [4   3   3   3   2   2   2   1   1   1   0]*L;
  node_z           = [0   0   0   0   0   0   0   0   0   0   0];
  consumption_node = [0   0   0   0   0   0   0   0   0   0 100];
  node_pressure    = [0   0   0   0   0   0   0   0   0   0   0];
  
  // Edge parameters
  // n°         1  2  3  4  5  6  7  8  9 10 11 12 13 14
  I_start   = [ 1  3  3  2  3  4  6  6  5  6  7  9  9  9];
  I_end     = [ 3  2  4  5  6  7  5  7  8  9 10  8 10 11];
  diameter  = [50 50 50 50 50 50 50 50 50 50 50 50 50 50]/1000; // Meters
  flow_pipe = [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0];
  
  InstanceName    = 'Test Hardy Cross 2';
  I_source        = [1];
  source_pressure = [10]; // 10 bars
  I_leaves        = [11];
elseif _type=='test_hc_3' then
  // Test from:
  // Hardy Cross, "Analysis of flow in networks of conduits or conductors", 
  // University of Illinois Engineering Experiment Station, Bulletin n° 286, november 1936

  L = 100;
  // Node parameters
  // Node n°          1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  15  17  18
  node_x           = [0   1   2   3   2   4   4   3   3   2   1   0   0   2   3 -0.5 4.5 2]*L;
  node_y           = [3   3   3   3   4   4   2.5 2.5 2   2   2   2   1   1   1  3.5 4.5 0.5]*L;
  node_z           = [0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  0   0   0];
  consumption_node = [0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  0   0 100];
  node_pressure    = [0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  0   0   0];
  
  // Edge parameters
  // Edge n°    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
  I_start   = [ 1  2  3 16  5  5 17  6  6  7  4  8  9 10 11  1  2  9 10 12 13 14 14];
  I_end     = [ 2  3  4  1  3  6  6  4  7  8  8  9 10 11 12 12 11 15 14 13 14 15 18];
  diameter  = [50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50]/1000;
  flow_pipe = [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0];
  
  InstanceName    = 'Test Hardy Cross 3';
  I_source        = [16 17];
  source_pressure = [10 15]; // 10 and 15 bars
  I_leaves        = [18];
elseif _type=='test_pr1' then
  // Node parameters
  // n°                  1    2    3    4      5   6   7   8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23
  node_x           = [-250 -250 -250    0    250 250 250 500  500 -500 -500 -500 -500    0  500    0 -750 -750 -750 -750 -250 -250    0];
  node_y           = [ 375  125 -125 -312.5 -125 125 375 375 -125  375  125 -125 -500 -500 -500 -750 -750 -500 -125  375  375  400 -800];
  node_z           = [  20   20   20   20     20  20  20 -20  -20  -20  -20  -20  -20  -20  -20    0    0    0    0    0  200   20    0];
  consumption_node = [   0    2    6    3      2   3   2   3    2    2    2    2    2    2    2    2    2    2    2    2    0    0    0];
  node_pressure    = [   0    0    0    0      0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0];
  
  // Edge parameters
  // n°         1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19
  //           20    21    22    23    24    25    26    27    28    29    30    31
  I_start   = [ 1     3     3     4     6     6     7    15     2     7     5     8     1    10    11    12    12    13     4 ...
               14     9    14    16    17    18    19    19    20     1    22    16];
  I_end     = [ 2     2     4     5     5     7     1    16     6     8     9     9    10    11    12     3    13    14    14 ...
               15    15    16    17    18    19    12    20    10    21     1    23];
  diameter  = [ 0.15  0.1   0.225 0.225 0.1   0.15  0.15  0.155 0.225 0.155 0.155 0.155 0.155 0.155 0.155 0.155 0.166 0.166 0.166 ...
                0.155 0.155 0.155 0.155 0.155 0.155 0.155 0.155 0.155 0.45  0.45  0.45]*1; // Meters
  flow_pipe = [ 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0  ...
                0     0     0     0     0     0     0     0     0     0     0     0];
                
  InstanceName    = 'Test problem 1';
  I_source        = [22];
  source_pressure = [10]; // 10 bars
  I_leaves        = [23];
elseif _type=='test_pr2' then
  // YC: ajouter un noeud source et un noeud client
  // Node parameters
  // n°                 1       2       3       4       5       6       7       8       9      10      11       12      13      14
  //                   15      16      17      18      19      20      21      22      23
  node_x           = [275     -25    -250     -25    -562.5   500    -562.5   500    -562.5   812.5  -562.5   -562.5   812.5   312.5   ... 
                      812.5   812.5  812.5     62.5  -250     275    -250      62.5    62.5];
  node_y           = [750     750     500     500     750     500     500     750     250     750    -312.5   -562.5   500    -562.5   ...
                      250     -62.5  -562.5  -312.5  -312.5   750     -62.5   -62.5   250];
  node_z           = [213.4   212.9   212.7   212.4   213.9   212.2   212.7   212.45  212.1   212.45  211.6    251.5   212     241.65  ...
                      215.75  221.7   231.7   231.7   211.65  511.7   511.9   271.85  261.95];
  consumption_node = [  5.207   8.525   5.602   7.097   3.418   6.908   6.835   3.418   6.568   1.3     9.875    3.327   3.875   1.143 ...
                       28.42    0       6.532   5.012   0       6.532   6.55    7.135  11.55];
  node_pressure    = [  0       0       0       0       0       0       0       0       0       0       0        0       0       0     ...
                        0       0       0       0       0       0       0       0       0];
  
  // Edge parameters
  // n°        1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18
  //          19    20    21    22    23    24    25    26    27    28    29    30    31    32
  I_start   = [1     2     5     7     9    11    12    14    18    19    19    18    21    22    23    23    23     7 ...
               3     4     6     8     2     8     6    14    13    13    15    16    16    20];
  I_end     = [2     5     7     9    11    12    14    18    19    11    21    22    22    23     9    15     3     3 ...
               4     6     8     1     4    10    13    17    10    15    16    17    22    1];
  diameter  = [0.155 0.155 0.155 0.155 0.175 0.185 0.185 0.185 0.175 0.175 0.155 0.175 0.155 0.155 0.155 0.165 0.155 0.155 ...
               0.155 0.155 0.155 0.155 0.155 0.155 0.165 0.2   0.155 0.175 0.185 0.2   0.2   0.255]*1; // Meters
  flow_pipe = [0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 ...
               0     0     0     0     0     0     0     0     0     0     0     0     0     0];
               
  InstanceName    = 'Test problem 2';
  I_source        = [1];
  source_pressure = [10]; // 10 bars
  I_leaves        = [11];
elseif _type=='test_pr3' then
  // YC: ajouter un noeud source et un noeud client
  // Node parameters
  // n°                1  2  3  4  5  6  7  8  9
  node_x           = [ 0  1  2  3  4  5  6  7  8]*100;
  node_y           = [ 0  2  0  0  2  0  0  2  0]*100;
  node_z           = [ 0  0  0  0  0  0  0  0  0];
  consumption_node = [10  0 10 10  0 10 10  0 10];
  node_pressure    = [ 0 10  0  0 10  0  0 10  0];
  
  // Edge parameters
  // n°        
  I_start   = [2 2 5 5 8 8];
  I_end     = [1 3 4 6 7 9];
  diameter  = [1 1 1 1 1 1]*0.155; // Meters
  flow_pipe = [0 0 0 0 0 0];
               
  InstanceName    = 'Test problem 3';
  I_source        = [2 5 8];
  source_pressure = [10 10 10]; // 10 bars
  I_leaves        = [1 3 4 6 7 9];
else
  error('wrong instance name');
end

// Set the pressure of the source
node_pressure(I_source) = source_pressure;

//
// Set the type of pipes:
//
// - 0 normal
// - 1 regulator

//
// Set the type of nodes:
//
// - 0 normal - always consum
// - 1 disconnectable - if P<P_threshold then we deconnect the consumer

NodeNum = length(node_x);

// Set the graph
g = make_graph(InstanceName,1,NodeNum,I_start,I_end);
g('node_x') = node_x;
g('node_y') = node_y;
// Node parameters
g = add_node_data(g,'node_z',node_z);
g = add_node_data(g,'flow',consumption_node);
g = add_node_data(g,'pressure',node_pressure);
// Pipe parameters
g = add_edge_data(g,'diameter',diameter);
g = add_edge_data(g,'flow',flow_pipe);

g('node_color') = [3 ones(size(node_x,1),size(node_x,2)-1)];
g('name') = InstanceName;
g('default_node_diam') = (max(g('node_x') - min(g('node_x')) + max(g('node_y') - min(g('node_y')))))/500;

if do_plot then
  params = init_param();
  params = add_param(params,'node_label',   %T);
  params = add_param(params,'node_pressure',%T);
  params = add_param(params,'node_angle',   45);
  params = add_param(params,'edge_label', %F);
  params = add_param(params,'edge_flow',  %T);
  params = add_param(params,'edge_angle', 45);

  plot_pipe_graph(g,InstanceName,params=params);
end
endfunction

