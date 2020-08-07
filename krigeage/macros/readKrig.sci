function kModel = readKrig(Filename)
  // Definition of the kriging structure
  // kModel('model')    : name of the model
  // kModel('tendency') : list of the tendency parameters
  // kModel('var')      : list of the variables
  // kModel('data')     : learning data set
  // kModel('K')        : K matrix
  // kModel('param')    : parameters of the covariance
  
  fid    = mopen(Filename,'r');
  kModel = mlist(['krig','model','tendency','var','data','K','param'],['','','',[],[],[]]);
  // Model
  String = mgetl(fid,1); kModel('model') = String;
  // Tendency
  String = mgetl(fid,1); kModel('tendency') = String;
  // var
  String = mgetl(fid,1); kModel('var') = String;
  // Data
  String = mgetl(fid,1); NbLines = eval(String);
  String = mgetl(fid,1); NbCol   = eval(String);
  DataMeas = zeros(NbLines,NbCol);
  for i=1:NbLines
    String = mgetl(fid,1); ListOfTokens = tokens(String,' ');
    for j=1:NbCol
      DataMeas(i,j) = eval(ListOfTokens(j));
    end
  end
  kModel('data') = DataMeas;
  // K
  String = mgetl(fid,1); NbLines = eval(String);
  String = mgetl(fid,1); NbCol   = eval(String);
  K = zeros(NbLines,NbCol);
  for i=1:NbLines
    String = mgetl(fid,1); ListOfTokens = tokens(String,' ');
    for j=1:NbCol
      K(i,j) = eval(ListOfTokens(j));
    end
  end
  kModel('K') = K
  // Model param
  String = mgetl(fid,1); NbCol = eval(String);
  String = mgetl(fid,1); ListOfTokens = tokens(String,' ');
  parm_sol = zeros(1,NbCol);
  for i=1:NbCol
    parm_sol(i) = eval(ListOfTokens(i));
  end
  kModel('param') = parm_sol;
  mclose(fid);
endfunction

