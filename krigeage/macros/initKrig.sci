function kModelOut = initKrig(Data,model,tendency,var)

if ~isdef('Data','local') then
  error('initKrig: Data is mandatory');
end
if ~isdef('model','local') then
  model = 'exp';
end
if ~isdef('tendency','local') then
  tendency = '1';
end
if ~isdef('var','local') then
  var = [];
  for i=1:size(Data,2)-1
    var = var + ' x'+string(i);
  end
end

kModelOut = mlist(['krig','model','tendency','var','data','K','param'],['','','',[],[],[]]);
// Model Name
kModelOut('model') = model;
// Tendency
kModelOut('tendency') = tendency;
// var
kModelOut('var') = var;
// Data
kModelOut('data') = Data;
endfunction
