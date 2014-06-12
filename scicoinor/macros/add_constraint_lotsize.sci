function param = add_constraint_lotsize(param,column,_range,_min,_max,lower,upper)
if ~isdef('param','local') then
  error('add_constraint_lotsize: param must be a clist');
end
if ~isdef('_max','local') then
  _max = [];
end
if (typeof(param)~='clist') then
  error('add_constraint_lotsize: param must be a clist');
end  
if and(lower>upper) then
  error('add_constraint_lotsize: upper must be higher than lower');
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

points = zeros(1,length(_min)+length(_max));
if isempty(_max) | ~_range then
  points = _min;
else
  points(1:2:$) = _min;
  points(2:2:$) = _max;
end

if (min(points)>lower(column)) | (max(points)<upper(column)) then
  error(sprintf('%s: error, you must fine tune the upper and lower bounds of column %d so as to match the given range\n','add_constraint_lotsize',column));
end

param('column')($) = column;
param('weight')($) = points;
param('length')($) = length(points);
param('range')($)  = _range;
if range==1 then
  res = and(_min<=_max);
  if ~res then
    error('add_constraint_lotsize: min is not less or equal to max');
  end
end
// 0 -> CbcSOS, 1-> CbcClique, 2->CbcLotsize, 3->CbcNWay, 4->CbcLink
param('id_obj')($) = 2; // CbcLotsize
endfunction

