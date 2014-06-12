function saveKrig(Filename, kModel, Plot)
  DataMeas = kModel('data');
  parm_sol = kModel('param');

  if ((~isempty(kModel('data')))&(~isempty(kModel('param')))) then
    fid = mopen(VarName,'w');

    // Writing model name
    mfprintf(fid,'%s\n',kModel('model'));
    // Writing tendency 
    mfprintf(fid,'%s\n',kModel('tendency'));
    // Writing var list
    mfprintf(fid,'%s\n',kModel('var'));
    // Writing DataMeas
    mfprintf(fid,'%d\n',size(kModel('data'),1));
    mfprintf(fid,'%d\n',size(kModel('data'),2));
    for i=1:size(kModel('data'),1)
      for j=1:size(kModel('data'),2)
        mfprintf(fid,'%f ',kModel('data')(i,j));
      end
      mfprintf(fid,'\n');
    end
    // Writing K
    mfprintf(fid,'%d\n',size(kModel('K'),1));
    mfprintf(fid,'%d\n',size(kModel('K'),2));
    for i=1:size(kModel('K'),1)
      for j=1:size(kModel('K'),2)
        mfprintf(fid,'%f ',kModel('K')(i,j));
      end
      mfprintf(fid,'\n');
    end
    // Writing parm_sol
    mfprintf(fid,'%d\n',length(kModel('param')));
    for i=1:length(kModel('param'))
      mfprintf(fid,'%f ',kModel('param')(i));
    end
    mclose(fid);
  else
    if (isempty(kModel('data'))) then
      if (Plot) then
        warningMessage('No dataset found');
      else
        printf('No dataset found\n');
      end
    end
    if (isempty(kModel('param'))) then
      if (Plot) then
        warningMessage('No model parameters found');
      else
        printf('No model parameters found\n');
      end
    end
  end
endfunction
