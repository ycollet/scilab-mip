///////////////////////////////////////
// export_model                      //
//                                   //
// export a model to scilab / matlab //
///////////////////////////////////////

function err = export_model(modelname, lolModel, file_type)

[nargout,nargin] = argn();

if nargin<2 then
  error('export_model: modelname and lolModel are mandatory. file_type is optional');
end

if ~isdef('file_type','local') then
  file_type = 'scilab';
end

if file_type=='scilab' then
  comment   = '//';
  extension = '.sci';
else
  comment = '%';
  extension = '.m';
end

[fid, err] = mopen(modelname+extension,'w');

if err then return; end

mfprintf(fid,'function y = %s(x)\n',modelname);
mfprintf(fid,'SumPhi = 0; y = 0;\n');

// Center of the partitions
mfprintf(fid,'%s Center of the partitions\n',comment);
for i=1:size(lolModel('listofcutinf'),1)
  mfprintf(fid,'xc(%d,:) = [',i);
  for j=1:size(lolModel('listofcutinf'),2)
    mfprintf(fid,' %f ', (lolModel('listofcutinf')(i,j) + lolModel('listofcutplus')(i,j))/2);
  end
  mfprintf(fid,'];\n');
end
mfprintf(fid,'\n');

// Size of the partitions
mfprintf(fid,'%s Size of the partitions\n',comment);
for i=1:size(lolModel('listofcutinf'),1)
  mfprintf(fid,'dx(%d,:) = [',i);
  for j=1:size(lolModel('listofcutinf'),2)
    mfprintf(fid,' %f ', (lolModel('listofcutplus')(i,j) - lolModel('listofcutinf')(i,j)));
  end
  mfprintf(fid,'];\n');
end
mfprintf(fid,'\n');

// Models of the partitions
mfprintf(fid,'%s Models of the partitions\n',comment);
for i=1:size(lolModel('listofmod'),1)
  mfprintf(fid,'models(%d,:) = [',i);
  for j=1:size(lolModel('listofmod'),2)
    mfprintf(fid,' %f ', lolModel('listofmod')(i,j));
  end
  mfprintf(fid,'];\n');
end
mfprintf(fid,'\n');

if lolModel('type')=='exp' then
  // Generation of the parameters related to the partitions - exponential membership function
  mfprintf(fid,'%s Computing the membership function values\n',comment);
  mfprintf(fid,'for i=1:%d\n',size(lolModel('listofcutinf'),1));
  mfprintf(fid,'  Phi(i) = 0;\n');
  mfprintf(fid,'  for j=1:%d\n',size(lolModel('listofcutinf'),2));
  mfprintf(fid,'    Phi(i) = Phi(i) + (x(j) - xc(i,j))^2/(sqrt(2)*%f*dx(i,j))^2;\n',lolModel('sigma'));
  mfprintf(fid,'  end\n');
  mfprintf(fid,'  Phi(i) = exp(-Phi(i));\n');
  mfprintf(fid,'  SumPhi = SumPhi + Phi(i)\n');
  mfprintf(fid,'end\n');
else
  // Generation of the parameters related to the partitions - piecewise linear membership function
  mfprintf(fid,'\nk_overlap = 0.5/%f\n',lolModel('sigma'));
  
  mfprintf(fid,'for i=1:%d\n',size(lolModel('listofcutinf'),1));
  mfprintf(fid,'  Phi(i) = 1.0;\n');
  mfprintf(fid,'  for j=1:%d\n',size(lolModel('listofcutinf'),2));
  mfprintf(fid,'    if (xc(i,j) - %f*dx(i,j)/2 < x(j)) & (x(j) < xc(i,j) + %f*dx(i,j)/2) then\n',lolModel('sigma'),lolModel('sigma'));
  mfprintf(fid,'      tmp = 1.0;\n');
  mfprintf(fid,'    elseif (xc(i,j) - dx(i,j) - k_overlap * %f*dx(i,j)/2 <= x(j)) & (x(j) <= xc(i,j) - %f*dx(i,j)/2) then\n',lolModel('sigma'),lolModel('sigma'));
  mfprintf(fid,'      tmp = (x(j) - xc(i,j) + dx(i,j) + k_overlap * %f*dx(i,j)/2) / (xc(i,j) - %f*dx(i,j)/2 - xc(i,j) + dx(i,j) + k_overlap * %f*dx(i,j)/2);\n',lolModel('sigma'),lolModel('sigma'),lolModel('sigma'));
  mfprintf(fid,'    elseif (xc(i,j) + %f*dx(i,j)/2 <= x(j)) & (x(j) <= xc(i,j) + dx(i,j) + k_overlap * %f*dx(i,j)/2) then\n', lolModel('sigma'),lolModel('sigma'));
  mfprintf(fid,'      tmp = (xc(i,j) + dx(i,j) + k_overlap * %f*dx(i,j)/2 - x(j)) / (xc(i,j) + dx(i,j) + k_overlap * %f*dx(i,j)/2 - xc(i,j) + %f*dx(i,j)/2);\n',lolModel('sigma'),lolModel('sigma'),lolModel('sigma'));
  mfprintf(fid,'    else\n');
  mfprintf(fid,'      tmp = 0;\n');
  mfprintf(fid,'    end\n');
  mfprintf(fid,'    Phi(i) = min(Phi(i), tmp);\n');
  mfprintf(fid,'  end\n');
  mfprintf(fid,'  SumPhi = SumPhi + Phi(i)\n');
  mfprintf(fid,'end\n');
end

if file_type=='scilab' then
  mfprintf(fid,'\nSumPhi = max(SumPhi,\%eps);\n\n');
else
  mfprintf(fid,'\nSumPhi = max(SumPhi,eps);\n\n');
end

// Generation of the Lolimot model
mfprintf(fid,'%s Computing the Lolimot model\n',comment);
mfprintf(fid,'for i=1:%d\n',size(lolModel('listofcutinf'),1));
mfprintf(fid,'  tmp = models(i,1);\n');
mfprintf(fid,'  for j=1:%d\n',size(lolModel('listofcutinf'),2));
mfprintf(fid,'    tmp = tmp + x(j)*models(i,j+1);\n');
mfprintf(fid,'  end\n');
mfprintf(fid,'  y = y + tmp*Phi(i)/SumPhi;\n');
mfprintf(fid,'end\n');

if file_type=='scilab' then
  mfprintf(fid,'endfunction\n');
end

mclose(fid);

res = err;
endfunction

