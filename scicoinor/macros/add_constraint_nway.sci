function param = add_constraint_nway(param,which)
if ~isdef('param','local') then
  error('add_constraint_nway: param must be a clist');
end
if (typeof(param)~='clist') then
  error('add_constraint_nway: param must be a clist');
end  

param('which')($+1)       = [];
param('weight')($+1)      = [];
param('type')($+1)        = 0;
param('clique_type')($+1) = 0;
param('column')($+1)      = 0;
param('range')($+1)       = 0;
param('id_obj')($+1)      = 0;
param('length')($+1)      = 0;
param('id')($+1)          = length(param('length'))-1;

param('which')($)  = which;
param('length')($) = length(which);
// 0 -> CbcSOS, 1-> CbcClique, 2->CbcLotsize, 3->CbcNWay, 4->CbcLink
param('id_obj')($) = 3; // CbcNWay
endfunction

