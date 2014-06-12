function param = add_constraint_sos(param,_type,which,weight,constr_no)
if ~isdef('param','local') then
  error('add_constraint_sos: param must be defined');
end
if (typeof(param)~='clist') then
  error('add_constraint_sos: param must be a clist');
end  
if ~isdef('constr_no','local') then
  constr_no = 0;
end

param('which')($+1)       = [];
param('weight')($+1)      = [];
param('type')($+1)        = 0;
param('clique_type')($+1) = 0;
param('column')($+1)      = 0;
param('range')($+1)       = 0;
param('id_obj')($+1)      = 0;
param('length')($+1)      = 0;
param('id')($+1)          = constr_no;

if isempty('weight') then
  weight = 1:length(which);
end

if ((_type<=0) |(_type>=3)) then
  error(sprintf('_type must be equal to 1 or 2 only\n','add_constraint_sos'));
end

//YC: add a new parameter (lower + upper for example) so as to check if the indexes given in which are correct.

param('type')($)   = _type;
param('which')($)  = which;
param('length')($) = length(which);
param('weight')($) = weight;

if (length(weight)~=length(which)) then
  error('add_constraint_sos: which and length params must have the same size');
end
// 0 -> CbcSOS, 1-> CbcClique, 2->CbcLotsize, 3->CbcNWay, 4->CbcLink
param('id_obj')($) = 0; // CbcSOS
endfunction

