lines(0);

stacksize('max');

write_html = %F;
mps_type   = 0;

mps_miplib = list();
mps_miplib(1) = 'data/misc07.mps';
mps_miplib(2) = 'data/reactor.mps';

if write_html then
  fid = mopen('listoffiles.html','w+');

  mfprintf(fid, '<TABLE BORDER=""1"">\n');
  mfprintf(fid, '<CAPTION> Characteritics of MPS files </CAPTION>\n');
  mfprintf(fid, '<TR>\n');
  mfprintf(fid, ' <TH> filename <H>\n');
  mfprintf(fid, ' <TH> index <H>\n');
  mfprintf(fid, ' <TH> nb_var <H>\n');
  mfprintf(fid, ' <TH> nb_constr <H>\n');
  mfprintf(fid, ' <TH> nb_coeff <H>\n');
  mfprintf(fid, ' <TH> nb_var_int <H>\n');
  mfprintf(fid, '<R>\n');
else
  fid = mopen('listoffiles.txt','w+');

  mfprintf(fid, 'Characteritics of MPS files\n');
  mfprintf(fid, '\n');
  mfprintf(fid, ' || filename || index || nb_var || nb_constr || nb_coeff || nb_var_int ||\n');
end

////////////////
// mps_miplib //
////////////////

for i=1:length(mps_miplib)
  mps_filename = mps_miplib(i);
  
  ///////////////////
  // Read the data //
  ///////////////////

  printf('Reading informations of file |%s|\n', mps_filename);
  mps_file = read_mps_file(mps_filename, mps_type);

  if write_html then
    mfprintf(fid, '<TR>\n');
    mfprintf(fid, ' <TH> %s <H>\n',mps_filename);
    mfprintf(fid, ' <TH> %d <H>\n',i);
    mfprintf(fid, ' <TH> %d <H>\n',mps_file('nb_obj_var'));
    mfprintf(fid, ' <TH> %d <H>\n',mps_file('nb_constr'));
    mfprintf(fid, ' <TH> %d <H>\n',mps_file('nb_val_constr_mat'));
    mfprintf(fid, ' <TH> %d <H>\n',mps_file('nb_int_var'));
    mfprintf(fid, '<R>\n');
  else
    mfprintf(fid, '|| %s ||',mps_filename);
    mfprintf(fid, ' %d ||',i);
    mfprintf(fid, ' %d ||',mps_file('nb_obj_var'));
    mfprintf(fid, ' %d ||',mps_file('nb_constr'));
    mfprintf(fid, ' %d ||',mps_file('nb_val_constr_mat'));
    mfprintf(fid, ' %d ||\n',mps_file('nb_int_var'));
  end
end

if write_html then
  mfprintf(fid,'<TABLE>\n');
end

mclose(fid);

