//////////////////////////////////////////////////////////////
// learn_lolimot                                            //
//                                                          // 
// This function build a Lolimot model for a given data set //
//////////////////////////////////////////////////////////////

function [modelOut,stat] = learn_lolimot(data,sigma,nbpart,maximp,nbCut,vec,Log,modelIn,pLinear)
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
if ~isdef('pLinear','local') then
  pLinear = %F;
end
if nbCut==1 then nbCut = 2; end

[nargout,nargin] = argn();

stat_is_defined = (nargout==2);

if isdef('modelIn','local') then
  if typeof(modelIn)~='lolimot' then
    error('learn_lolimot: error - modelIn argument must be a ''lolimot'' structure\n');
  end

  modelOut = modelIn;

  if modelOut('nbdim')~=size(data,2)-1 then
    error('learn_lolimot: nbdim of the model doesn''t correspond to the numbre of dimensions in the data set\n');
  end

  modelOut('sigma') = sigma;
  
  modelOut = _learn_model(data,1,modelOut,pLinear);
  modelOut = update_lolimot(data,modelOut,vec);
  
  old_residual = modelOut('residual');  
else
  // Initialisation of the Lolimot structure
  modelOut = mlist(['lolimot','type','nbdim','sigma','listofmod','listofcutinf','listofcutplus','listofresidual','residual', 'min', 'max', ...
                              [],    [],     [],     [],         [],           [],              [],              []          []     []]);

  // Initialisation of the first partition
  modelOut('nbdim')         = size(data,2)-1;
  modelOut('listofcutplus') = max(data(:,1:$-1),'r');
  modelOut('listofcutinf')  = min(data(:,1:$-1),'r');
  modelOut('sigma')         = sigma;
  modelOut('min')           = modelOut('listofcutinf');
  modelOut('max')           = modelOut('listofcutplus');
  if pLinear==%T then
    modelOut('type')        = 'lin';
  else
    modelOut('type')        = 'exp';  
  end

  modelOut = _learn_model(data,1,modelOut,pLinear);
  modelOut = update_lolimot(data,modelOut,vec);

  old_residual = modelOut('residual');
end

if stat_is_defined then
  stat = list();
end

for i=1:nbpart-1
  t_start = getdate();
  
  if Log then
    printf('learn_lolimot: partition %d / %d - residual = %.2f', i, nbpart,modelOut('residual'));
  end
  
  if stat_is_defined then
    stat($+1) = [];
    stat($)($+1) = modelOut('residual'); // Residual before adding a partition
  end
  
  // Find the partition with the maximum residual
  [tmp,Index] = max(modelOut('listofresidual')); Index = Index(1);
  // We copy the selected cut at the end of the list
  modelOut('listofmod')     = [modelOut('listofmod');     modelOut('listofmod')(Index,:)];
  modelOut('listofcutinf')  = [modelOut('listofcutinf');  modelOut('listofcutinf')(Index,:)];
  modelOut('listofcutplus') = [modelOut('listofcutplus'); modelOut('listofcutplus')(Index,:)];
  // Test some cuttings
  Index_ToCut = 0;
  Max_ToCut   = %inf;
  for j=1:modelOut('nbdim')
    for k=2:nbCut
      mid_pos = (modelOut('listofcutplus')(Index,j) - modelOut('listofcutinf')(Index,j)) * (k-1) / nbCut + modelOut('listofcutinf')(Index,j);
      modelOut('listofcutinf')($,j)      = mid_pos;
      modelOut('listofcutplus')(Index,j) = mid_pos;
      modelOut = update_lolimot(data,modelOut,vec);
      if modelOut('residual')<Max_ToCut then
        Max_ToCut   = modelOut('residual');
        Index_ToCut = j;
        CutPos      = k;
      end
      // Undo the change
      modelOut('listofcutinf')($,j)      = modelOut('listofcutinf')(Index,j);
      modelOut('listofcutplus')(Index,j) = modelOut('listofcutplus')($,j);
    end
  end
  // Perform again the best cutting
  mid_pos = (modelOut('listofcutplus')(Index,Index_ToCut) - modelOut('listofcutinf')(Index,Index_ToCut)) * (CutPos-1) / nbCut + modelOut('listofcutinf')(Index,Index_ToCut);
  modelOut('listofcutinf')($,Index_ToCut)      = mid_pos;
  modelOut('listofcutplus')(Index,Index_ToCut) = mid_pos;
  modelOut = update_lolimot(data,modelOut,vec);

  t_end = getdate();

  if Log then
    printf(' / %.2f - improvement = %3.2f %% - elapsed time = %.2f s\n', ...
           modelOut('residual'),100*abs(old_residual - modelOut('residual')) / max([modelOut('residual') %eps]), etime(t_end,t_start));
  end

  if stat_is_defined then
    stat($)($+1) = modelOut('residual'); // Residual after adding a partition
    stat($)($+1) = etime(t_end,t_start); // time required for learning
    stat($)($+1) = Index; // Which partition has been cut
    stat($)($+1) = Index_ToCut; // The cut dimension
  end

  if abs(old_residual - modelOut('residual')) / max([modelOut('residual') %eps]) < maximp then
    if Log then
      printf('learn_lolimot: %f improvement reached - learning phase stopped\n', maximp);
    end
    break;
  end
  
  old_residual = modelOut('residual');
end
endfunction

