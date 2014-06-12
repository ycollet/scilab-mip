//////////////////////////////////////////////////////////////////////////////////////////////
// learn_valid_lolimot                                                                      //
//                                                                                          //
// This function build a Lolimot model for a given data set and a given validation data set //
//////////////////////////////////////////////////////////////////////////////////////////////

function [modelOut,stat] = learn_valid_lolimot(data_learn, data_valid,sigma,nbpart,maximp,nbCut,vec,Log,modelIn,pLinear)
if ~isdef('sigma','local') then
  sigma = 0.33;
end
if ~isdef('nbpart','local') then
  nbpart = 10;
end
if ~isdef('maximp','local') then
  maximp = 0;
end
if ~isdef('nbCut','local') then
  nbCut = 2;
end
if ~isdef('vec','local') then
  vec = %F;
end
if ~isdef('Log','local') then
  Log = %F;
end
if ~isdef('data_learn','local') then
  error('learn_valid_lolimot: data_learn is required\n');
end
if ~isdef('data_valid','local') then
  error('learn_valid_lolimot: data_valid is required\n');
end
if size(data_valid,2)~=size(data_learn,2) then
  error('learn_valid_lolimot: data_learn and data_valid must have the same number of dimensions\n');
end
if ~isdef('pLinear','local') then
  pLinear = %F;
end

if nbCut==1 then nbCut = 2; end

[nargout,nargin] = argn();

stat_is_defined = (nargout==2);

if stat_is_defined then
  stat = list();
end

if isdef('modelIn','local') then
  // Learn a model with 1 partition
  if stat_is_defined then
    [modelOut,stat_aux] = learn_lolimot(data_learn,sigma,2,maximp,nbCut,vec,Log,modelIn,pLinear);
    stat = lstcat(stat,stat_aux);
  else
    modelOut = learn_lolimot(data_learn,sigma,2,maximp,nbCut,vec,Log,modelIn,pLinear);
  end
else
  // Learn a model with 1 partition
  if stat_is_defined then
    [modelOut,stat_aux] = learn_lolimot(data_learn,sigma,2,maximp,nbCut,vec,Log,pLinear=pLinear);
    stat = lstcat(stat,stat_aux);
  else
    modelOut = learn_lolimot(data_learn,sigma,2,maximp,nbCut,vec,Log,pLinear=pLinear);
  end
end
// Perform the validation
Res = 0;
if vec then
  for j=1:size(data_valid,1)
    Res = Res + (estim_vec_lolimot(data_valid(j,1:$-1),modelOut) - data_valid(j,$))^2;
  end
else
  for j=1:size(data_valid,1)
    Res = Res + (estim_lolimot(data_valid(j,1:$-1),modelOut) - data_valid(j,$))^2;
  end
end

Min_Valid = Res;
if Log then
  printf('learn_valid_lolimot: Partition %d / %d - validation residual = %.2f\n',1,nbpart,Min_Valid);
end

for i=2:nbpart
  old_model = modelOut;
  // Learn a model with part+1 partition
  if stat_is_defined then
    [modelOut,stat_aux] = learn_lolimot(data_learn,sigma,2,maximp,nbCut,vec,Log,modelOut,pLinear);
    stat = lstcat(stat,stat_aux);
  else
    modelOut = learn_lolimot(data_learn,sigma,2,maximp,nbCut,vec,Log,modelOut,pLinear);
  end
  // Perform the validation
  Res = 0;
  if Vec then
    for j=1:size(data_valid,1)
      Res = Res + (estim_vec_lolimot(data_valid(j,1:$-1),modelOut) - data_valid(j,$))^2;
    end
  else
    for j=1:size(data_valid,1)
      Res = Res + (estim_lolimot(data_valid(j,1:$-1),modelOut) - data_valid(j,$))^2;
    end
  end
  
  if Res<Min_Valid then
    Min_Valid = Res;
    if Log then
      printf('learn_valid_lolimot: Partition %d / %d - validation residual = %.2f\n',i,nbpart,Min_Valid);
    end
  else
    if Log then
      printf('learn_valid_lolimot: Partition %d / %d - validation residual = %.2f - best residual %.2f - learning phase stopped\n',i,nbpart,Res,Min_Valid);
    end
    // We switch back to the previous model which was better wrt the validation residual
    modelOut = old_model;
    break;
  end
end
endfunction

