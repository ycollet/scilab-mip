//////////////////////////////////////////////////
// export_der_model                             //
//                                              //
// export a derivative model to scilab / matlab //
//////////////////////////////////////////////////

function err = export_der_model(modelname, lolModel, file_type)

[nargout,nargin] = argn();

if nargin<2 then
  error('export_der_model: modelname and lolModel are mandatory. file_type is optional');
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
    mfprintf(fid,' %f ', (lolModel('listofcutplus')(i,j) - lolModel('listofcutinf')(i,j))*lolModel('sigma'));
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
  // Exponential model
  mfprintf(fid,'%s Computation of the derivative of the Lolimot model with exponential membership function\n',comment);
  mfprintf(fid,'u = 0; v = 0;\n');
  mfprintf(fid,'for i=1:%d\n',size(lolModel('listofcutplus'),1));
  mfprintf(fid,'  Li = (models(i,1) + sum(x(:)'' .* models(i,2:$)));\n');
  mfprintf(fid,'  u  = u + Li * exp(-sum((x(:)'' - xc(i,:)'')./dx(i,:)'').^2);\n');
  mfprintf(fid,'  v  = v + exp(-sum((x(:)'' - xc(i,:)'')./dx(i,:)'').^2);\n');
  mfprintf(fid,'end\n');

  mfprintf(fid,'for i=1:%d\n',size(lolModel('listofcutplus'),2));
  mfprintf(fid,'  du = 0; dv = 0;\n');
  mfprintf(fid,'  for j=1:%d\n',size(lolModel('listofcutplus'),1));
  mfprintf(fid,'    Lj = (models(j,1) + sum(x(:)'' .* models(j,2:$)));\n');
  mfprintf(fid,'    dv = dv - 2 * ((x(i) - xc(j,i)) / dx(j,i)^2) * exp(-sum((x(:)'' - xc(j,:)'')./dx(j,:)'').^2);\n');
  mfprintf(fid,'    du = du + (models(j,i+1) - 2 * ((x(i) - xc(j,i)) / dx(j,i)^2) * Lj) * exp(-sum((x(:)'' - xc(j,:)'')./dx(j,:)'').^2);\n');
  mfprintf(fid,'  end\n');
  if file_type=='scilab' then
    mfprintf(fid,'  y(i,1) = (du * v - dv * u) / max(v^2,\%eps);\n');
  else
    mfprintf(fid,'  y(i,1) = (du * v - dv * u) / max(v^2,eps);\n');
  end
  mfprintf(fid,'end\n');
else
  mfprintf(fid,'%s Computation of the derivative of the Lolimot model with piecewise linear membership function\n',comment);
  mfprintf(fid,'k_overlap = 0.5 / %f\n',lolModel('sigma'));
  mfprintf(fid,'L = []; Mf = []; dMf = []; U = 0; dU = []; V = 0; dV = [];\n');
  mfprintf(fid,'for i=1:%d\n',size(lolModel('listofcutplus'),1));
  mfprintf(fid,'  for j=1:%d\n',size(lolModel('listofcutplus'),2));
  mfprintf(fid,'    if ((xc(i,j) - dx(i,j)/2 < x(j))&(x(j) < xc(i,j) + dx(i,j)/2)) then\n');
  mfprintf(fid,'      phi_aux = 1.0;\n');
  mfprintf(fid,'    elseif ((xc(i,j) - dx(i,j)/%f/2 - k_overlap * dx(i,j)/2 <= x(j))&(x(j) <= xc(i,j) - dx(i,j)/2)) then\n',lolModel('sigma'));
  mfprintf(fid,'      phi_aux = (x(j) - (xc(i,j) - dx(i,j)/%f/2 - k_overlap *dx(i,j)/2)) / ... \n',lolModel('sigma'));
  mfprintf(fid,'                ((xc(i,j) - dx(i,j)/2) - (xc(i,j) - dx(i,j)/%f / 2 - k_overlap * dx(i,j)/2));\n',lolModel('sigma'));
  mfprintf(fid,'    elseif ((xc(i,j) + dx(i,j)/2 <= x(j))&(x(j) <= xc(i,j) + dx(i,j)/%f/2 + k_overlap * dx(i,j)/2)) then\n',lolModel('sigma'));
  mfprintf(fid,'      phi_aux = (x(j) - (xc(i,j) + dx(i,j)/%f/2 + k_overlap * dx(i,j)/2)) / ...\n',lolModel('sigma'));
  mfprintf(fid,'                ((xc(i,j) + dx(i,j)/2) - (xc(i,j) + dx(i,j)/%f/2 + k_overlap * dx(i,j)/2))\n',lolModel('sigma'));
  mfprintf(fid,'    else\n');
  mfprintf(fid,'      phi_aux = 0;\n');
  mfprintf(fid,'    end\n');
  mfprintf(fid,'    Mf(i,j) = phi_aux;\n');
  mfprintf(fid,'  end\n');
  
  mfprintf(fid,'  Li = models(i,1) + sum(x(:)'' .* models(i,2:$));\n');
    
  mfprintf(fid,'  [_tmp,Index] = min(Mf(i,:)); Index = Index(1);\n');

  mfprintf(fid,'  dMf_tmp = 0;\n');
  mfprintf(fid,'  if ((xc(i,Index) - dx(i,Index)/%f/2 - k_overlap * dx(i,Index)/2 <= x(Index))&(x(Index) <= xc(i,Index) - dx(i,Index)/2)) then\n',lolModel('sigma'));
  mfprintf(fid,'    dMf_tmp = 1 / ((xc(i,Index) - dx(i,Index)/2) - (xc(i,Index) - dx(i,Index)/%f/2 - k_overlap * dx(i,Index)/2));\n',lolModel('sigma'));
  mfprintf(fid,'  end\n');
  mfprintf(fid,'  if ((xc(i,Index) + dx(i,Index)/2 <= x(Index))&(x(Index) <= xc(i,Index) + dx(i,Index)/%f/2 + k_overlap * dx(i,Index)/2)) then\n',lolModel('sigma'));
  mfprintf(fid,'    dMf_tmp = 1 / ((xc(i,Index) + dx(i,Index)/2) - (xc(i,Index) + dx(i,Index)/%f/2 + k_overlap * dx(i,Index)/2));\n',lolModel('sigma'));
  mfprintf(fid,'  end\n');
  mfprintf(fid,'  dMf(i,:)     = zeros(1,length(x));\n');
  mfprintf(fid,'  dMf(i,Index) = dMf_tmp;\n');
  mfprintf(fid,'  U = U + Mf(i,Index)*Li;\n');i
  mfprintf(fid,'  V = V + Mf(i,Index);\n');
  mfprintf(fid,'  dU = dU + models(i,2:$) * Mf(i,Index) + dMf(i,:) * Li;\n');
  mfprintf(fid,'  dV = dV + dMf(i,:);\n');
  mfprintf(fid,'end\n');
  if file_type=='scilab' then
    mfprintf(fid,'y = (dU*V-U*dV)/max(V^2,\%eps);\n');
  else
    mfprintf(fid,'y = (dU*V-U*dV)/max(V^2,eps);\n');
  end
end

if file_type=='scilab' then
  mfprintf(fid,'endfunction\n');
end

mclose(fid);

res = err;
endfunction

