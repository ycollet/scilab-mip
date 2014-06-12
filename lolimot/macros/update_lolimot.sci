//////////////////////////////////////////////////////////////////////////
// update_lolimot                                                       //
//                                                                      //
// We recompute the linear model of each partition in the Lolimot model //
//////////////////////////////////////////////////////////////////////////

function modelOut = update_lolimot(data,lolModel,vec)

modelOut = lolModel;

if ~isdef('vec','local') then
  vec = %F;
end

// Update the linear model of each partition
for i=1:size(modelOut('listofcutinf'),1)
  modelOut = _learn_model(data,i,modelOut);
end

// Compute the global residual
tmp_res = 0;
if vec then
  for i=1:size(data,1)
    tmp_res = tmp_res + (data(i,$) - estim_vec_lolimot(data(i,1:$-1),modelOut))^2;
  end
else
  for i=1:size(data,1)
    tmp_res = tmp_res + (data(i,$) - estim_lolimot(data(i,1:$-1),modelOut))^2;
  end
end
modelOut('residual') = tmp_res;
endfunction

