//////////////////////////////////////////////////////////////////////////////////////
// krig_ui                                                                          //
//                                                                                  //
// This function displays an UI interface to manage the learning of a kriging model //
//////////////////////////////////////////////////////////////////////////////////////

function krig_ui()

global kModel;
global Dimension;
global Dimension_Old;
global DataValid;
global k_var;
global m_text;
global StepDistance;
global listOfModels;
global ExportString_First;
global ExportString;
global Horizon;
global NbRestart;

// Global variables for UI
global hf1;
//global m_read;
//global m_read_learn;
//global m_read_valid;
//global m_read_model;
//global m_krig;
//global m_lean;
//global m_validate;
//global m_param;
//global m_model;
//global m_tend;
//global m_var;
//global m_stepdist;
//global m_horizon;
//global m_restart;
//global m_export_file;
//global m_quit;
global m_text;

kModel = mlist(['krig','model','tendency','var','data','K','param'],['','','',[],[],[]]);
// Definition of the kriging structure
// kModel('model')    : name of the model
// kModel('tendency') : list of the tendency parameters
// kModel('var')      : list of the variables
// kModel('data')     : learning data set
// kModel('K')        : K matrix
// kModel('param')    : parameters of the covariance

// Initialisation of the kriging structure

kModel('model')    = 'exp';
kModel('tendency') = '1';
kModel('var')      = 'x1';
kModel('data')     = [];
kModel('K')        = [];
kModel('param')    = [];

// Build the list of models

listOfModels = locate_correl();

// Default values

Dimension     = 1;
Dimension_Old = 1;
k_var         = 2;
StepDistance  = 0.1;
DataValid     = [];
Horizon       = 1;
NbRestart     = 10;

// Building the UI interface

hf1           = scf();
drawlater;
delmenu(hf1.figure_id,gettext('&Fichier'));
delmenu(hf1.figure_id,gettext('&Outils'));
delmenu(hf1.figure_id,gettext('&Edition'));
delmenu(hf1.figure_id,gettext('&?'));
toolbar(hf1.figure_id,'off');

hf1.position  = [10 10 320 100];
hf1.figure_name = 'Krigeage';

m_read        = uimenu(hf1,'label','Read');
m_read_learn  = uimenu(m_read,'label','Measures','callback','ReadKrigMeasures');
m_read_valid  = uimenu(m_read,'label','Validation','callback','ReadKrigValidation');
m_read_model  = uimenu(m_read,'label','Model','callback','ReadKrigModel');
m_krig        = uimenu(hf1,'label','Krigeage');
m_lean        = uimenu(m_krig,'label','learn','callback','LearnKrig');
m_validate    = uimenu(m_krig,'label','validate','callback','ValidateKrig');
m_param       = uimenu(hf1,'label','Parameters');
m_model       = uimenu(m_param,'label','Set model','callback','SetModel');
m_tend        = uimenu(m_param,'label','Set tendency','callback','SetTendency');
m_var         = uimenu(m_param,'label','Set list of var','callback','SetListOfVar');
m_stepdist    = uimenu(m_param,'label','Set distance step','callback','SetStepDist');
m_horizon     = uimenu(m_param,'label','Set learn horizon length','callback','SetHorizon');
m_restart     = uimenu(m_param,'label','Set number of restart','callback','SetNbRestart');
m_export_file = uimenu(hf1,'label','Export','callback','SetExportFile');
m_quit        = uimenu(hf1,'label','Close','callback','delete(gcf())');
m_text        = uicontrol(hf1,'style','text', ...
                         'position',[20 10 280 80], ...
                         'string','Krigeage', ...
                         'backgroundcolor',[1.0,1.0,1.0]);
drawnow;

endfunction

///////////////////////////////////////////////////
// ReadKrigMeasures                              //
//                                               //
// Display a dialog and read a learning data set //
///////////////////////////////////////////////////

function ReadKrigMeasures()
  global kModel;
  global m_text;
  global Dimension;
  
  set(m_text,'string','Reading measures');
  
  Filename = tk_getfile('*.dat','.',title='Krigeage - Read measurements');
  if (Filename~='') then
    DataMeas = read(Filename, -1, Dimension+1); // Dimension dimensions but Dimension+1 columns (Dimension for input and 1 for output)
    kModel('data') = DataMeas;
  end
  set(m_text,'string','Krigeage');
endfunction

/////////////////////////////////////////////////////
// ReadKrigValidation                              //
//                                                 //
// Display a dialog and read a validation data set //
/////////////////////////////////////////////////////

function ReadKrigValidation()
  global DataValid;
  global m_text;
  global Dimension;
  
  set(m_text,'string','Reading validation');

  Filename  = tk_getfile('*.dat','.',title='Krigeage - Read validation');
  if (Filename~='') then
    DataValid = read(Filename, -1, Dimension+1);
  end
  set(m_text,'string','Krigeage');
endfunction

///////////////////////////////////////////////////////
// ReadKrigModel                                     //
//                                                   //
// Display a dialog and select the correlation model //
///////////////////////////////////////////////////////

function ReadKrigModel()
  global kModel;
  global m_text;
  
  set(m_text,'string','Reading model');
  
  Filename = tk_getfile('*.dat','.',title='Krigeage - Read model');
  if (Filename~='') then
    kModel = readKrig(Filename);
  end  
  set(m_text,'string','Krigeage');
endfunction

/////////////////////////////////////////////
// warning message                         //
//                                         //
// Display a dialog with a warning message //
/////////////////////////////////////////////

function warningMessage(wMessage)
  hf2 = scf();
  drawlater;
  hf2.position    = [10 10 250 100];
  hf2.figure_name = 'Warning';
  toolbar(hf2.figure_id,'off');
  delmenu(hf2.figure_id,gettext('&Fichier'));
  delmenu(hf2.figure_id,gettext('&Outils'));
  delmenu(hf2.figure_id,gettext('&Edition'));
  delmenu(hf2.figure_id,gettext('&?'));
  hf2_edit = uicontrol(hf2,'style','text','position',[20 60 230 40],'string',wMessage);
  MyString_OK = 'delete(gcf())';
  hf2_but_OK = uicontrol(hf2,'style','pushbutton', ...
                          'position',[100 20 50 20], ...
                          'string','OK', ...
                          'callback',MyString_OK);    
  drawnow;
endfunction

///////////////////////////////////////////////////////////////////////
// updateModel                                                       //
//                                                                   //
// update the name of the correlation model in the kriging structure //
///////////////////////////////////////////////////////////////////////

function updateModel(my_hf, my_edit_hf)
  global kModel;
  global listOfModels;
  
  kModel('model') = listOfModels(get(my_edit_hf,'value'));
  
  close(my_hf);
endfunction

///////////////////////////////////////////////////////////////////////
// SetModel                                                          //
//                                                                   //
// update the name of the correlation model in the kriging structure //
///////////////////////////////////////////////////////////////////////

function SetModel()
  global m_text;
  global listOfModels;
  
  set(m_text,'string','Updating model');

  hf2 = scf();
  drawlater;
  hf2.position    = [10 10 200 160];
  hf2.figure_name = 'Set model';
  toolbar(hf2.figure_id,'off');
  delmenu(hf2.figure_id,gettext('&Fichier'));
  delmenu(hf2.figure_id,gettext('&Outils'));
  delmenu(hf2.figure_id,gettext('&Edition'));
  delmenu(hf2.figure_id,gettext('&?'));
  
  hf2_edit = uicontrol(hf2,'style','listbox','position',[20 60 150 80],'string',listOfModels);
  
  MyString_OK     = 'updateModel(' + string(hf2) + ',' + string(hf2_edit) + ')';
  MyString_Cancel = 'delete(gcf())';
  
  hf2_but_OK = uicontrol(hf2,'style','pushbutton', ...
                          'position',[20 20 50 20], ...
                          'string','OK', ...
                          'callback',MyString_OK);
  hf2_but_Cancel = uicontrol(hf2,'style','pushbutton', ...
                             'position',[120, 20, 50, 20], ...
                             'string','Cancel', ...
                             'callback',MyString_Cancel);
                             
  set(m_text,'string','Krigeage');
  drawnow;
endfunction

////////////////////////////////////////////////////////////////////
// updateTendency                                                 //
//                                                                //
// update the expression of the tendency in the kriging structure //
////////////////////////////////////////////////////////////////////

function updateTendency(my_hf, my_edit_hf)
  global kModel;
  
  kModel('tendency') = get(my_edit_hf,'string');
  if isempty(kModel('tendency')) then
    // The tendency can't be empty. We must have, a minima, a constant value
    kModel('tendency') = '1';
  end
  close(my_hf);
endfunction

////////////////////////////////////////////////////////////////////
// SetTendency                                                    //
//                                                                //
// update the expression of the tendency in the kriging structure //
////////////////////////////////////////////////////////////////////

function SetTendency()
  global m_text;
  global kModel;
  
  set(m_text,'string','Updating tendency');

  hf2 = scf();
  drawlater;
  hf2.position    = [10 10 200 100];
  hf2.figure_name = 'Set tendency';
  toolbar(hf2.figure_id,'off');
  delmenu(hf2.figure_id,gettext('&Fichier'));
  delmenu(hf2.figure_id,gettext('&Outils'));
  delmenu(hf2.figure_id,gettext('&Edition'));
  delmenu(hf2.figure_id,gettext('&?'));
  hf2_edit = uicontrol(hf2,'style','edit','position',[20 60 150 20],'string',kModel('tendency'));
  
  MyString_OK     = 'updateTendency(' + string(hf2) + ',' + string(hf2_edit) + ')';
  MyString_Cancel = 'delete(gcf())';
  
  hf2_but_OK = uicontrol(hf2,'style','pushbutton', ...
                          'position',[20 20 50 20], ...
                          'string','OK', ...
                          'callback',MyString_OK);
  hf2_but_Cancel = uicontrol(hf2,'style','pushbutton', ...
                             'position',[120, 20, 50, 20], ...
                             'string','Cancel', ...
                             'callback',MyString_Cancel);

  set(m_text,'string','Krigeage');
  drawnow;
endfunction

////////////////////////////////////////////////////////////////////////
// updateListOfVar                                                    //
//                                                                    //
// update the expression of the list of vars in the kriging structure //
////////////////////////////////////////////////////////////////////////

function updateListOfVar(my_hf, my_edit_hf)
  global m_text;
  global kModel;
  global Dimension;
  global Dimension_Old;
  
  Dimension_Old = size(tokens(kModel('var'),' '),1);

  kModel('var') = get(my_edit_hf,'string');
  if isempty(kModel('var')) then
    // We must have, a minima, one variable
    kModel('var') = 'x1';
  end
  
  Dimension = size(tokens(kModel('var'),' '),1);
    
  close(my_hf);

  if (Dimension_Old~=Dimension) then
    clear DataMeas;
    clear DataValid;
    warningMessage('You must read again \n learn and validation files');
    Dimension_Old = Dimension;
  end
  
  set(m_text,'string','Krigeage');
endfunction

////////////////////////////////////////////////////////////////////////
// SetListOfVar                                                       //
//                                                                    //
// update the expression of the list of vars in the kriging structure //
////////////////////////////////////////////////////////////////////////

function SetListOfVar()
  global m_text;
  global kModel;
  
  set(m_text,'string','Updating list of var');

  hf2 = scf();
  drawlater;
  hf2.position    = [10 10 200 100];
  hf2.figure_name = 'Set list of var';
  toolbar(hf2.figure_id,'off');
  delmenu(hf2.figure_id,gettext('&Fichier'));
  delmenu(hf2.figure_id,gettext('&Outils'));
  delmenu(hf2.figure_id,gettext('&Edition'));
  delmenu(hf2.figure_id,gettext('&?'));
  
  hf2_edit = uicontrol(hf2,'style','edit','position',[20 60 150 20],'string',kModel('var'));
  
  MyString_OK     = 'updateListOfVar(' + string(hf2) + ',' + string(hf2_edit) + ')';
  MyString_Cancel = 'delete(gcf())';
  
  hf2_but_OK = uicontrol(hf2,'style','pushbutton', ...
                          'position',[20 20 50 20], ...
                          'string','OK', ...
                          'callback',MyString_OK);
  hf2_but_Cancel = uicontrol(hf2,'style','pushbutton', ...
                             'position',[120, 20, 50, 20], ...
                             'string','Cancel', ...
                             'callback',MyString_Cancel);

  set(m_text,'string','Krigeage');
  drawnow;
endfunction

////////////////////////////////////////////////////////////////////
// updateStepDistance                                             //
//                                                                //
// update the value of the step distance in the kriging structure //
////////////////////////////////////////////////////////////////////

function updateStepDistance(my_hf, my_edit_hf)
  global StepDistance;
  
  StepDistance = eval(get(my_edit_hf,'string'));
  
  close(my_hf);
endfunction

////////////////////////////////////////////////////////////////////
// SetStepDistance                                                //
//                                                                //
// update the value of the step distance in the kriging structure //
////////////////////////////////////////////////////////////////////

function SetStepDist()
  global m_text;
  global StepDistance;
  
  set(m_text,'string','Updating step distance');

  hf2 = scf();
  drawlater;
  hf2.position    = [10 10 200 100];
  hf2.figure_name = 'Set step distance';
  toolbar(hf2.figure_id,'off');
  delmenu(hf2.figure_id,gettext('&Fichier'));
  delmenu(hf2.figure_id,gettext('&Outils'));
  delmenu(hf2.figure_id,gettext('&Edition'));
  delmenu(hf2.figure_id,gettext('&?'));
  
  hf2_edit = uicontrol(hf2,'style','edit','position',[20 60 150 20],'string',string(StepDistance));
  
  MyString_OK     = 'updateStepDistance(' + string(hf2) + ',' + string(hf2_edit) + ')';
  MyString_Cancel = 'delete(gcf()))';
  
  hf2_but_OK = uicontrol(hf2,'style','pushbutton', ...
                          'position',[20 20 50 20], ...
                          'string','OK', ...
                          'callback',MyString_OK);
  hf2_but_Cancel = uicontrol(hf2,'style','pushbutton', ...
                             'position',[120, 20, 50, 20], ...
                             'string','Cancel', ...
                             'callback',MyString_Cancel);

  set(m_text,'string','Krigeage');
  drawnow;
endfunction

//////////////////////////////////////////////////////////////
// updateHorizon                                            //
//                                                          //
// update the value of the horizon in the kriging structure //
//////////////////////////////////////////////////////////////

function updateHorizon(my_hf, my_edit_hf)
  global Horizon;
  
  Horizon = eval(get(my_edit_hf,'string'));
  
  close(my_hf);
endfunction

//////////////////////////////////////////////////////////////
// SetHorizon                                               //
//                                                          //
// update the value of the horizon in the kriging structure //
//////////////////////////////////////////////////////////////

function SetHorizon()
  global m_text;
  global Horizon;
  
  set(m_text,'string','Update learn horizon length');

  hf2 = scf();
  drawlater;
  hf2.position    = [10 10 200 100];
  hf2.figure_name = 'Set Horizon';
  toolbar(hf2.figure_id,'off');
  delmenu(hf2.figure_id,gettext('&Fichier'));
  delmenu(hf2.figure_id,gettext('&Outils'));
  delmenu(hf2.figure_id,gettext('&Edition'));
  delmenu(hf2.figure_id,gettext('&?'));
  
  hf2_edit = uicontrol(hf2,'style','edit','position',[20 60 150 20],'string',string(Horizon));
  
  MyString_OK     = 'updateHorizon(' + string(hf2) + ',' + string(hf2_edit) + ')';
  MyString_Cancel = 'delete(gcf())';
  
  hf2_but_OK = uicontrol(hf2,'style','pushbutton', ...
                          'position',[20 20 50 20], ...
                          'string','OK', ...
                          'callback',MyString_OK);
  hf2_but_Cancel = uicontrol(hf2,'style','pushbutton', ...
                             'position',[120, 20, 50, 20], ...
                             'string','Cancel', ...
                             'callback',MyString_Cancel);

  set(m_text,'string','Krigeage');
  drawnow;
endfunction

////////////////////////////////////////////////////////////////////////////
// updateNbRestart                                                        //
//                                                                        //
// update the number of random restart for learning the correlation model //
////////////////////////////////////////////////////////////////////////////

function updateNbRestart(my_hf, my_edit_hf)
  global NbRestart;
  
  NbRestart = eval(get(my_edit_hf,'string'));
  
  close(my_hf);
endfunction

////////////////////////////////////////////////////////////////////////////
// SetNbRestart                                                           //
//                                                                        //
// update the number of random restart for learning the correlation model //
////////////////////////////////////////////////////////////////////////////

function SetNbRestart()
  global m_text;
  global NbRestart;
  
  set(m_text,'string','Update number of restart');

  hf2 = scf();
  drawlater;
  hf2.position    = [10 10 200 100];
  hf2.figure_name = 'Set Number of restart';
  toolbar(hf2.figure_id,'off');
  delmenu(hf2.figure_id,gettext('&Fichier'));
  delmenu(hf2.figure_id,gettext('&Outils'));
  delmenu(hf2.figure_id,gettext('&Edition'));
  delmenu(hf2.figure_id,gettext('&?'));
  
  hf2_edit = uicontrol(hf2,'style','edit','position',[20 60 150 20],'string',string(NbRestart));
  
  MyString_OK     = 'updateNbRestart(' + string(hf2) + ',' + string(hf2_edit) + ')';
  MyString_Cancel = 'delete(gcf())';
  
  hf2_but_OK = uicontrol(hf2,'style','pushbutton', ...
                          'position',[20 20 50 20], ...
                          'string','OK', ...
                          'callback',MyString_OK);
  hf2_but_Cancel = uicontrol(hf2,'style','pushbutton', ...
                             'position',[120, 20, 50, 20], ...
                             'string','Cancel', ...
                             'callback',MyString_Cancel);

  set(m_text,'string','Krigeage');
  drawnow;
endfunction

/////////////////////////////////////////////////
// ExportModelFile                             //
//                                             //
// Set the name of the kriging model to export //
/////////////////////////////////////////////////

function ExportModelFile(my_hf, my_edit_hf)
  global kModel;
  
  VarName = get(my_edit_hf,'string');

  saveKrig(VarName, kModel);
  
  close(my_hf);
endfunction

/////////////////////////////////////////////////
// SetExportFile                               //
//                                             //
// Set the name of the kriging model to export //
/////////////////////////////////////////////////

function SetExportFile()
  global m_text;

  set(m_text,'string','Export model to file');

  hf2 = scf();
  drawlater;
  hf2.position    = [10 10 200 100];
  hf2.figure_name = 'Export to file';
  toolbar(hf2.figure_id,'off');
  delmenu(hf2.figure_id,gettext('&Fichier'));
  delmenu(hf2.figure_id,gettext('&Outils'));
  delmenu(hf2.figure_id,gettext('&Edition'));
  delmenu(hf2.figure_id,gettext('&?'));
  
  hf2_edit = uicontrol(hf2,'style','edit','position',[20 60 150 20],'string','krig_name.dat');
  
  MyString_OK     = 'ExportModelFile(' + string(hf2) + ',' + string(hf2_edit) + ')';
  MyString_Cancel = 'delete(gcf())';
  
  hf2_but_OK = uicontrol(hf2,'style','pushbutton', ...
                          'position',[20 20 50 20], ...
                          'string','OK', ...
                          'callback',MyString_OK);
  hf2_but_Cancel = uicontrol(hf2,'style','pushbutton', ...
                             'position',[120, 20, 50, 20], ...
                             'string','Cancel', ...
                             'callback',MyString_Cancel);

  set(m_text,'string','Krigeage');
  drawnow;
endfunction

/////////////////////////////////////////////////
// ValidateKrig                                //
//                                             //
// Perform the validation of the kriging model //
/////////////////////////////////////////////////

function ValidateKrig()
  global m_text;
  global kModel;
  global DataValid;
  
  set(m_text,'string','Validation of the model');

  if ((size(DataValid,1)==0)&(size(kModel('data'),1)==0)) then
    warningMessage('You must read again \n validation file');
    return;
  elseif (size(DataValid,1)==0) then
    warningMessage('You must read again \n validation file');
    return;
  elseif (size(kModel('data'),1)==0) then
    warningMessage('You must read again \n learn file');
    return;
  end
   
  ieee_flag = ieee();
  ieee(2);

  WId = waitbar(0, 'Validation of the model');
  
  X = kModel('data')(:,1:$-1);
  Y = kModel('data')(:,$);
  NbPoints = size(kModel('data'),1);
  
  // We split the tendancy into a list of terms which will be used to build the matrixes P and p_x
  ListOfTerms = tokens(kModel('tendency'),' ');
  ListOfVar   = tokens(kModel('var'),' ');
  
  // Building matrix P
  // We replace all the le xj by X(i,j)
  for i=1:size(ListOfTerms,1)
    for j=1:size(ListOfVar,1)
      ListOfTerms(i) = strsubst(ListOfTerms(i),ListOfVar(j),'X(i,' + string(j) + ')');
    end
  end

  P = zeros(NbPoints,size(ListOfTerms,1));
  for i=1:NbPoints
    for j=1:size(ListOfTerms,1)
      P(i,j) = eval(ListOfTerms(j));
    end
  end

  Index_x = 0;

  k_x = zeros(NbPoints,1);
  
  for k=1:size(DataValid,1)
    waitbar(floor(1000*k/size(DataValid,1))/1000, WId);
  
    x = DataValid(k,1:$-1);
    
    Index_x = Index_x + 1;

    [T_mod(Index_x), var(Index_x)] = computeKrig(kModel, x);
  end

  winclose(WId);
  
  // If there is juste one input, we plot y=f(x). Otherwise, we plot y=f(n° of point)
  if (size(DataValid,2)==2) then
    X_mod  = DataValid(:,1);
    X_name = 'x';
  else
    X_mod  = 1:size(DataValid,1);
    X_name = 'rang';
  end
  
  hg = scf();
  drawlater;
  set(hg,'figure_name','Krigeage - Validation - Plot %d');

  plot(X_mod, T_mod, 'r');
  plot(X_mod, DataValid(:,size(DataValid,2)), 'g');
  plot(X_mod, T_mod + k_var*var, 'k-.');
  plot(X_mod, T_mod - k_var*var, 'k-.');
  xtitle('Estimation / Mesure - Krigeage',X_name,'y');
  strplus  = sprintf('Estimation + %d * var',k_var);
  strmoins = sprintf('Estimation - %d * var',k_var);
  legend(['Estimation','Mesure',strplus,strmoins]);
  drawnow;
  set(m_text,'string','Krigeage');

  ieee(ieee_flag);
  
  R2 = 1 - sum((DataValid(:,$) - T_mod).^2) / size(DataValid,1) / (stdev(DataValid(:,$))^2);
  printf('Validation R2 = %f\n', R2);
endfunction

///////////////////////////////////////////////
// LearnKrig                                 //
//                                           //
// Perform the learning of the kriging model //
///////////////////////////////////////////////

function LearnKrig()
  global m_text;

  global kModel;
  global Horizon;
  global StepDistance;
  global NbRestart;
  
  if (size(kModel('data'),1)==0) then
    warningMessage('You must read again \n learn file');
  end

  set(m_text,'string','Learning the model');

  kModel = learnKrig(kModel, %T, StepDistance, Horizon, NbRestart);

  set(m_text,'string','Krigeage');
endfunction

