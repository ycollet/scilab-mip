function [x0, lower, upper] = pipe_init(g_in, I_sources,params)

if ~isdef('params','local') then
  params = [];
end

[Qmax,err]   = get_param(params,'Qmax', 250);
[Pmax,err]   = get_param(params,'Pmax', 30);
[Pmin,err]   = get_param(params,'Pmin', 0);
[k_pipe,err] = get_param(params, 'k_pipe', 5/100);
[Debug,err]  = get_param(params, 'Debug', %F);

n_nodes  = g_in.node_number;
n_pipes  = length(g_in.tail);

//////////////////////////
// Initialisation of x0 //
//////////////////////////

///////////////////////
// Variable ordering //
///////////////////////
// 1           to n_nodes           - the pressure at node i
// n_nodes + 1 to n_nodes + n_pipes - the flow in pipe i

// Random initialisation
x0 = [(Pmax - Pmin).*rand(1,n_nodes) + Pmin 2*Qmax*rand(1,n_pipes) - Qmax]';

// Set the pressure
for i=1:length(I_sources)
  x0(I_sources(i)) = g_in.nodes.data.pressure(I_sources(i));
end

// Set the consumer flow
for i=1:n_pipes
  x0(n_nodes + i) = g_in.nodes.data.flow(i);
end

//////////////////////////////////////////////////
// Definition of the upper and lower boundaries //
//////////////////////////////////////////////////

lower = [ 0*ones(1,n_nodes)    -Qmax*ones(1,n_pipes)]';
upper = [ Pmax*ones(1,n_nodes)  Qmax*ones(1,n_pipes)]';
endfunction

function g_out = pipe_update_graph(g_in, x_opt)

n_nodes  = g_in.node_number;
n_pipes  = length(g_in.tail);

///////////////////////
// Variable ordering //
///////////////////////
// 1           to n_nodes           - the pressure at node i
// n_nodes + 1 to n_nodes + n_pipes - the flow in pipe i

// We put the results in the output graph
g_out = g_in;
for i=1:n_nodes
  g_out.nodes.data.pressure(i) = x_opt(i);
end
// Update the flow consummed
for i=1:n_nodes
  g_out.nodes.data.flow(i) = x_opt(i);
end
for i=1:n_pipes
  g_out.edges.data.flow(i) = x_opt(n_nodes+i);
end
endfunction

//////////////////////////////
// The equality constraints //
//////////////////////////////

function h = compute_pipe_eq(x,g_in,I_sources,params)

if ~isdef('params','local') then
  params = [];
end

[Qmax,err]   = get_param(params,'Qmax', 250);
[Pmax,err]   = get_param(params,'Pmax', 30);
[Pmin,err]   = get_param(params,'Pmin', 0);
[k_pipe,err] = get_param(params, 'k_pipe', 5/100);
[Debug,err]  = get_param(params, 'Debug', %F);

n_nodes  = g_in.node_number;
n_pipes  = length(g_in.tail);

///////////////////////
// Variable ordering //
///////////////////////
// 1            to n_nodes           - the pressure at node i
// n_nodes + 1  to n_nodes + n_pipes - the flow in pipe i

// Generate the constraints
h  = [];
Index_constraint = 1;

if Debug then
  printf('Equality constraints\n');
end

///////////////////////////////////////////////////////////////////////
// Set 1: We run across all the nodes and produce the node equations //
///////////////////////////////////////////////////////////////////////
for i=1:n_nodes
  // We skip the source node
  if ~isempty(find(I_sources==i)) then continue; end
  
  Index_begin = find(g_in.tail==i);
  Index_end   = find(g_in.head==i);

  ///////////////////////////////////////////
  // Set 1-b : Produce the node constraint //
  ///////////////////////////////////////////
  h(Index_constraint) = -g_in.nodes.data.flow(i);
  if Debug then
    printf('normal consumer - %d - (node %d) - %f ', Index_constraint, i, -g_in.nodes.data.flow(i));
  end

  ////////////////////////////////////////////
  // Set 1-c : continue the node constraint //
  ////////////////////////////////////////////
  if ~isempty(Index_begin) then h(Index_constraint) = h(Index_constraint) - sum(x(n_nodes + Index_begin)); end
  if ~isempty(Index_end)   then h(Index_constraint) = h(Index_constraint) + sum(x(n_nodes + Index_end)); end

  Index_constraint = Index_constraint + 1;
  if Debug then    
    if ~isempty(Index_begin) then
      for ii=1:length(Index_begin)
        printf('- pipe[%d %d] ', g_in.tail(Index_begin(ii)), g_in.head(Index_begin(ii)));
      end
    end
    if ~isempty(Index_end) then
      for ii=1:length(Index_end)
        printf('+ pipe[%d %d] ', g_in.tail(Index_end(ii)), g_in.head(Index_end(ii)));
      end
    end
    printf('\n');
  end
end

////////////////////////////////////////////////////////////////////////
// Set 2: We run across all the pipes and produce the pipes equations //
////////////////////////////////////////////////////////////////////////
for i=1:n_pipes
  // Get the flow in pipe i
  Index_begin = g_in.tail(i);
  Index_end   = g_in.head(i);

  Q1 = x(n_nodes+i); // Flow in pipe
  P1 = x(Index_begin);
  P2 = x(Index_end);

  L  = sqrt((g_in.node_x(Index_begin) - g_in.node_x(Index_end))^2 + (g_in.node_y(Index_begin) - g_in.node_y(Index_end))^2);
  D  = g_in.edges.data.diameter(i)*1000;
  
  // Produce the pipe equation
  h(Index_constraint) = P1^2 - P2^2 - k_pipe * L / D^5 * Q1 * abs(Q1);
  if Debug then
    printf('Case 1 - constraint %d: pressure node %d ^2 - pressure node %d ^2 - flow pipe [%d %d] = 0\n', Index_constraint, Index_begin, Index_end, Index_begin, Index_end);
  end
  Index_constraint = Index_constraint + 1;
end

//////////////////////////////////////////////////////////////////////////////
// Set 3: Deal with the source node: at source node, the pressure is fixed. //
//////////////////////////////////////////////////////////////////////////////
for i=1:length(I_sources)
  h(Index_constraint) = x(I_sources(i)) - g_in.nodes.data.pressure(I_sources(i));
  if Debug then
    printf('constraint %d: pressure node %d - %f (%f) = 0\n', Index_constraint, I_sources(i), x(I_sources(i)), g_in.nodes.data.pressure(I_sources(i)));
  end
  Index_constraint = Index_constraint + 1;
end
endfunction

function dh = compute_pipe_deq(x,g_in,I_sources,params)

if ~isdef('params','local') then
  params = [];
end

[Qmax,err]   = get_param(params,'Qmax', 250);
[Pmax,err]   = get_param(params,'Pmax', 30);
[Pmin,err]   = get_param(params,'Pmin', 0);
[k_pipe,err] = get_param(params, 'k_pipe', 5/100);
[Debug,err]  = get_param(params, 'Debug', %F);

n_nodes  = g_in.node_number;
n_pipes  = length(g_in.tail);

///////////////////////
// Variable ordering //
///////////////////////
// 1           to n_nodes           - the pressure at node i
// n_nodes + 1 to n_nodes + n_pipes - the flow in pipe i

// Generate the constraints
dh = [];
Index_constraint = 1;

///////////////////////////////////////////////////////////////////////
// Set 1: We run across all the nodes and produce the node equations //
///////////////////////////////////////////////////////////////////////
for i=1:n_nodes
  // We skip the source node
  if ~isempty(find(I_sources==i)) then continue; end

  Index_begin = find(g_in.tail==i);
  Index_end   = find(g_in.head==i);

  dh(Index_constraint,:) = spzeros(size(x,1),size(x,2))';

  ///////////////////////////////////////////
  // Set 1-b : Produce the node constraint //
  ///////////////////////////////////////////
  dh(Index_constraint,:) = spzeros(size(x,1),size(x,2))';

  ////////////////////////////////////////////
  // Set 1-c : continue the node constraint //
  ////////////////////////////////////////////
  if ~isempty(Index_begin) then 
    dh(Index_constraint,n_nodes + Index_begin) = dh(Index_constraint,n_nodes + Index_begin) - ones(1,length(Index_begin)); 
  end
  if ~isempty(Index_end)   then 
    dh(Index_constraint,n_nodes + Index_end) = dh(Index_constraint,n_nodes + Index_end) + ones(1,length(Index_end));   
  end
  Index_constraint = Index_constraint + 1;
end

////////////////////////////////////////////////////////////////////////
// Set 2: We run across all the pipes and produce the pipes equations //
////////////////////////////////////////////////////////////////////////
for i=1:n_pipes
  // Get the flow in pipe i
  Index_begin = g_in.tail(i);
  Index_end   = g_in.head(i);
  
  Q1 = x(n_nodes+i); // Flow in pipe
  P1 = x(Index_begin);
  P2 = x(Index_end);
  
  L  = sqrt((g_in.node_x(Index_begin) - g_in.node_x(Index_end))^2 + (g_in.node_y(Index_begin) - g_in.node_y(Index_end))^2);
  D  = g_in.edges.data.diameter(i)*1000;

  // Produce the pipe equation
  dh(Index_constraint,:)            = spzeros(size(x,1),size(x,2))';
  dh(Index_constraint, Index_begin) =  2 * P1;
  dh(Index_constraint, Index_end)   = -2 * P2;
  dh(Index_constraint, n_nodes + i) = k_pipe * L / D^5 * (-abs(Q1)^2-Q1^2)/max(%eps,abs(Q1));
  Index_constraint = Index_constraint + 1;
end

//////////////////////////////////////
// Set 3: Deal with the source node //
//////////////////////////////////////
dh(Index_constraint,:) = spzeros(size(x,1),size(x,2))';
for i=1:length(I_sources)
  dh(Index_constraint, I_sources(i)) = 1;
  Index_constraint = Index_constraint + 1;
end

dh = sparse(dh');
endfunction

