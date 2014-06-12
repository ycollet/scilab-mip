/////////////////////////////////////////////////////////
// _learn_model                                        //
//                                                     //
// Learn the linear model corresponding to a partition //
/////////////////////////////////////////////////////////
function modelOut = _learn_model(data,nopart,lolModel,pLinear)
  modelOut = lolModel;
  // We build the weighted data set for a partition
  if pLinear then
    for i=1:size(data,1)
      phi_coeff  = _phi_norm_data_linear(data(i,1:($-1)),nopart,modelOut);
      Input(i,:) = phi_coeff*[1 data(i,1:($-1))];
      Output(i)  = phi_coeff*data(i,$);
    end
  else
    for i=1:size(data,1)
      phi_coeff  = _phi_norm_data(data(i,1:($-1)),nopart,modelOut);
      Input(i,:) = phi_coeff*[1 data(i,1:($-1))];
      Output(i)  = phi_coeff*data(i,$);
    end
  end
  lin_model = Input \ Output;
  tmp_res = 0;
  for i=1:size(Input,1)
    tmp_res = tmp_res + (lin_model' * Input(i,:)' - Output(i))^2;
  end
  // Updating lolModel
  modelOut('listofmod')(nopart,:)    = lin_model';
  modelOut('listofresidual')(nopart) = tmp_res;
endfunction
