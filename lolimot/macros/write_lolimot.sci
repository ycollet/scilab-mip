///////////////////////////////////
// write_lolimot                 //
//                               //
// Store a learnt model on drive //
///////////////////////////////////

function err = write_lolimot(Filename, lolModel)
[fid, err] = mopen(Filename,'w');

if err then return; end;
mfprintf(fid,'%d\n',lolModel('nbdim'));
mfprintf(fid,'%f\n',lolModel('sigma'));
mfprintf(fid,'%f\n',lolModel('residual'));
mfprintf(fid,'%d\n',size(lolModel('listofmod'),1)); // The number of partitions
// We save the list of models
for i=1:size(lolModel('listofmod'),1)
  for j=1:size(lolModel('listofmod'),2)
    mfprintf(fid,'%f ',lolModel('listofmod')(i,j));
  end
  mfprintf(fid,'\n');
end
// We save the list of inf cuts
for i=1:size(lolModel('listofcutinf'),1)
  for j=1:size(lolModel('listofcutinf'),2)
    mfprintf(fid,'%f ',lolModel('listofcutinf')(i,j));
  end
  mfprintf(fid,'\n');
end
// We save the list of inf plus
for i=1:size(lolModel('listofcutplus'),1)
  for j=1:size(lolModel('listofcutplus'),2)
    mfprintf(fid,'%f ',lolModel('listofcutplus')(i,j));
  end
  mfprintf(fid,'\n');
end
// We save the list of residuals
for i=1:size(lolModel('listofresidual'),1)
  mfprintf(fid,'%f ',lolModel('listofresidual')(i));
end
mfprintf(fid,'\n');
// We save the min and max bounds
for i=1:length(lolModel('min'))
  mfprintf(fid,'%f ',lolModel('min')(i));
end
mfprintf(fid,'\n');
for i=1:length(lolModel('max'))
  mfprintf(fid,'%f ',lolModel('max')(i));
end
mfprintf(fid,'\n');
mfprintf(fid,'%s\n',lolModel('type'));
mclose(fid);
endfunction

