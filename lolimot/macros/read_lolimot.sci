/////////////////////////////////////////////////////
// read_lolimot                                    //
//                                                 //
// Read a lolimot model previously stored on drive //
/////////////////////////////////////////////////////

function [err,modelOut] = read_lolimot(Filename)
[fid, err] = mopen(Filename,'r');
if err then return; end;

modelOut = mlist(['lolimot','type','nbdim','sigma','listofmod','listofcutinf','listofcutplus','listofresidual','residual', 'min', 'max', ...
                            [],    [],     [],     [],         [],           [],              [],              []          []     []]);

txt = mgetl(fid,1); modelOut('nbdim')    = eval(txt);
txt = mgetl(fid,1); modelOut('sigma')    = eval(txt);
txt = mgetl(fid,1); modelOut('residual') = eval(txt);
txt = mgetl(fid,1); nopart = eval(txt);

// We read the list of models
modelOut('listofmod') = [];
for i=1:nopart
  txt = mgetl(fid,1);
  modelOut('listofmod') = [modelOut('listofmod'); eval(tokens(txt,' '))'];
end
// We read the list of inf cuts
modelOut('listofcutinf') = [];
for i=1:nopart
  txt = mgetl(fid,1);
  modelOut('listofcutinf') = [modelOut('listofcutinf'); eval(tokens(txt,' '))'];
end
// We save the list of inf plus
modelOut('listofcutplus') = [];
for i=1:nopart
  txt = mgetl(fid,1);
  modelOut('listofcutplus') = [modelOut('listofcutplus'); eval(tokens(txt,' '))'];
end
// We save the list of residuals
txt = mgetl(fid,1);
modelOut('listofresidual') = eval(tokens(txt,' '))';
// We read the min bounds
txt = mgetl(fid,1);
modelOut('min') = eval(tokens(txt,' '))';
// We read the max bounds
txt = mgetl(fid,1);
modelOut('max') = eval(tokens(txt,' '))';
txt = mgetl(fid,1);
modelOut('type') = txt;
mclose(fid);
endfunction

