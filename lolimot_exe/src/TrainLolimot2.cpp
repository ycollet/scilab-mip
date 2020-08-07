// Author: Y. Collette
// EMail: ycollet@freesurf.fr
// Date: 02/04/2007

/*
This file is part of Lolimot.

    Foobar is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <unistd.h>
#include <string.h>

#if !defined(WIN32)
#include <dlfcn.h>

#define dlhandle void *
#endif

#include <fstream>
#include <string>
#include <vector>
#include <ios>
#include <iostream>
#include <cstdlib>
#include <functional>
#include <limits>
#include <sstream>

#include <tokenize.hpp>

using namespace std;

#include <LL_Dimension.h>
#include <LL_Mesure.h>
#include <LL_Lolimot.h>
#include <LL_Partition.h>

#include <Learning.h>
#include <Filter.h>
#include <DefaultParam.h>
#include <DisplayData.h>
#include <Analysis.h>
#include <Post.h>
#include <ExportC.h>
#include <ExportCpp.h>
#include <ExportMatlab.h>
#include <TrainLolimotStruct.h>
#include <Validate.h>

#define NBMAXPARTITIONS     100
#define RESIDUGAPPERCENTAGE 0.0
#define LL_SIGMA            0.33
#define LL_NBCUT            2

#define BUFFERSIZE          8196

//#define DEBUG           1

#if !defined(WIN32)
typedef float (*lolimot_t)(vector<float> &);
#endif

// Structure is_comma
// Defines the delimitor as the comma during the analysis of the data files

struct is_comma : public std::unary_function<char,bool>
{
  bool operator() (char c) const
  { return (c==',');}
};

int main(int argc, char ** argv)
{
  vector<vector<float> >  InputData, AuxData;
  vector<vector<float> >  MeasureData;
  vector<float>           WeightData;
  vector<vector<float> >  TempoInput, TempoOutput;
  vector<vector<string> > List_TRAIN_FileName, List_VALID_FileName;
  vector<vector<unsigned int> > List_TRAIN_SkipLines, List_VALID_SkipLines;
  vector<string>          Aux_FileName;
  vector<unsigned int>    Aux_SkipLines, List_Retro_Input, Offset_For_Files;
  vector<Filter_Entry>    List_Filter;
  vector<Delta_Entry>     List_Delta;
  vector<PartToCut>       List_PartToCut;
  vector<Model_Coeff>     List_ModelCoeff;
  vector<Data_Skip>       List_SkipLearn, List_SkipValid;
  string                  Line, ModelName, ResultDir, ResName, AuxStr, Res_Filename, R2_Filename, Part_Filename;
  string                  Load_Model_Filename, Save_Model_Filename, LOO_Filename, SampledDataFile, Norm_Filename;
  string                  Current_Dir, Working_Dir, CustParamStr, CustRetBefStr, CustRetAfterStr, CrossValid_Filename;
  string                  DistribOfCutName;
  ifstream              * DataFile   = NULL;
  ifstream              * ConfigFile = NULL;
  ofstream              * ExportFile = NULL;
  ofstream              * ModelFile_Save = NULL;
  ifstream              * ModelFile_Load = NULL;
  vector<vector<string> > WList;
  vector<float>           DList;
  vector<float>           InputMin, MeasureMin;
  vector<float>           InputMax, MeasureMax;
  vector<float>           Learn_InputMin, Learn_MeasureMin;
  vector<float>           Learn_InputMax, Learn_MeasureMax;
  vector<float>           X_Test, Y_LL;
  vector<Model_Input>     List_Model_Input, List_Model_Output;
  Model_Input             One_Model_Input, One_Model_Output, Weight_Input;
  PartToCut               Aux_PartToCut;
  Model_Coeff             Aux_ModelCoeff;
  Data_Skip               Aux_DataSkip;
  stringstream            FileName;
  char                 ** Buffer    = NULL;
  char                  * CfgBuffer = NULL;
  LL_Lolimot            * Lolimot   = NULL;
  float                   StdResidu = 0.0;
  int                     NbMaxPart = -1, NbCut = -1, Measure = -1, Filter_Step = -1, MAverage = -1, Select_Step = -1;
  int                     CrossValid_NbFiles = -1;
  float                   Sigma = -1.0, ResiduGapPercentage  = -1.0;
  bool                    Result, Validation_File, BS_Light = false, Display = false;
  bool                    BS_ComputeR2 = false, BS_Classic = false, BS_NoMeth = false;
  bool                    UseDistribCutting = false, UseUniformCutting = true;
  bool                    UseOptimizeSigma = false, UsePartFilename = false, Learn_Recur = false, UseUpdate = false;
  bool                    SampledData_Exit = false, UseNormalize = false, UseDelta = false; 
  bool                    PrintPart = false, UseTransformLog = false, Use_SkipValid = false;
  bool                    Echant_Output_UseMin = true, Use_LOO = false, UseSampledData = true, PrintDistribOfCut = false;
  bool                    ExportC = false, ExportCpp = false, ExportMatlab = false, NbColSet = false, Use_SkipLearn = false;
  bool                    D_ExportC = false, D_ExportCpp = false, D_ExportMatlab = false, Echant_AddMinMax = false;
  bool                    DD_ExportC = false, DD_ExportCpp = false, DD_ExportMatlab = false;
  bool                    D_ExportC2 = false, D_ExportCpp2 = false, D_ExportMatlab2 = false, UseLasso = false;
  bool                    DD_ExportC2 = false, DD_ExportCpp2 = false, DD_ExportMatlab2 = false;
  bool                    Compute_C1 = false, Compute_C2 = false, Save_Model = false, Load_Model = false;
  bool                    CustomParam = false, CustomRet_Before = false, CustomRet_After = false, CrossValid = false;
  int                     Echant_Stat_NbBar = -1, Echant_Stat_NbPoints = -1, Echant_Dyn_NbBar = -1;
  int                     Echant_Delta_NbPoints = -1, Echant_Delta_NbTimes = -1, Echant_PartGreedDelta_NbPoints = -1;
  int                     Echant_Rand_NbPoints = -1, Echant_Greedy_NbPoints = -1, Echant_PartGreedDelta_NbTimes = -1;
  int                     Echant_Greedy_Fast_NbPointsAtTheEnd = -1, Echant_Greedy_Fast_NbPointsToSelect = -1;
  int                     Echant_Greedy_Add_Measure_NbPoints = -1;
  int                     Echant_Greedy_Fast_Add_Measure_NbPointsAtTheEnd = -1;
  int                     Echant_Greedy_Fast_Add_Measure_NbPointsToSelect = -1;
  int                     Echant_Output_NbPts = -1, Echant_Greedy_Partial_NbPtsToSelect = -1;
  int                     Nb_Iter = -1, Echant_Greedy_Partial_NbTimesToSelect = -1;
  int                     NbPartToExplore = -1, Echant_Dyn_NbPoints = -1;
  float                   BS_Proportion = -1, SigmaMin = 0.0, SigmaStep = 0.1, SigmaMax = 1.0;
  float                   CrossValid_StdResidual, CrossValid_MeanResidual, UseTransformLogAlpha = 1.0;
  float                   UseTransformLogEps = 0.0, MembershipThreshold = 0.0;
  int                     Rec_NbPart_Ph1 = -1, Rec_NbPart_Ph2 = -1, Valid_Start = -1, Valid_End = -1;
  int                     Learn_Start = -1, Learn_End = -1;
  unsigned int            NbMaxTempo = numeric_limits<unsigned int>::min(), RandSeed = 12345, ShiftInput = 0;
  unsigned int            ShiftOutput = 0;
  int                     ResidualType = -1, LineCount = 0, NbCol = 0, NbColAux = 0;
  unsigned int            i, j, k, n, m;
  struct Filter_Entry     Aux_Entry;

  Working_Dir = "./";
  Current_Dir = "./";

  CfgBuffer = new char[BUFFERSIZE];

  getcwd(CfgBuffer, BUFFERSIZE);
  Current_Dir = CfgBuffer;
  Current_Dir += "/";

  delete [] CfgBuffer;

  BS_Light = false;

  // We write the results in the current directory by default
  ResultDir = "./";

  // Some default parameters
  Res_Filename        = "Lolimot_Validation_Residual.dat";
  R2_Filename         = "Lolimot_R2.dat";
  Part_Filename       = "Lolimot_ListOfPart.dat";
  Load_Model_Filename = "";
  Save_Model_Filename = "Lolimot_Model.mod";

  MeasureData.resize(0);
  InputData.resize(0);
  WeightData.resize(0);

  List_Filter.resize(0);
  List_PartToCut.resize(0);
  List_ModelCoeff.resize(0);
  List_SkipLearn.resize(0);
  List_SkipValid.resize(0);

  // Opening the configuration file
  ConfigFile = new ifstream;

  CfgBuffer = new char[BUFFERSIZE];

  if (argc>1)
    {
      ConfigFile->open(argv[1], ios_base::in);
    } /* End If */
  else
    {
      cerr << "Error: you must give a config file - Usage: ./TrainLolimot2 ConfigFile" << endl;
      exit(1);
    } /* End Else */

  if (!ConfigFile->is_open())
    {
      cerr << "Error: can't open file " << argv[1] << "." << endl;
      exit(1);
    } /* End Else */

  List_Model_Output.resize(0);
  List_Model_Input.resize(0);
  
  List_TRAIN_FileName.resize(0);
  List_VALID_FileName.resize(0);

  List_TRAIN_SkipLines.resize(0);
  List_VALID_SkipLines.resize(0);

  List_Retro_Input.resize(0);

  WList.resize(1);

  while(ConfigFile->getline(CfgBuffer, BUFFERSIZE, '\n'))
    {
      if (strlen(CfgBuffer)==0) continue;
      
      if (Display) cout << "Buffer = " << CfgBuffer << endl;
      WList[0] = tokenize(CfgBuffer);

      if (WList[0][0][0]=='%')
	{
	  // C'est un commentaire, on ne fait rien
	} /* End If */
      else if (WList[0].size()==0)
	{
	  // C'est une ligne vide, on ne fait rien
	} /* End Else If */
      else if (WList[0][0]=="MODELNAME")
	{
	  if (WList[0].size()==1)
	    {
	      cerr << "ERROR: MODELNAME - Filename is missing" << endl;
	      exit(1);
	    } /* End If */
	  ModelName = WList[0][1];
	} /* End If */
      else if (WList[0][0]=="TRAIN")
	{
	  if (WList[0].size()==1)
	    {
	      cerr << "ERROR: TRAIN - Learning data filename is missing" << endl;
	      exit(1);
	    } /* End If */

	  Aux_FileName.clear();
	  Aux_SkipLines.clear();
	  for(i=1; i<WList[0].size(); i+=2)
	    {
	      Aux_FileName.push_back(WList[0][i]);
	      Aux_SkipLines.push_back((unsigned int)atoi(WList[0][i+1].c_str()));
	    } /* End For */
	  List_TRAIN_FileName.push_back(Aux_FileName);
	  List_TRAIN_SkipLines.push_back(Aux_SkipLines);
	} /* End Else If */
      else if (WList[0][0]=="VALID")
	{
	  if (WList[0].size()==1)
	    {
	      cerr << "ERROR: VALID - Validation data filename is missing" << endl;
	      exit(1);
	    } /* End If */

	  Aux_FileName.clear();
	  Aux_SkipLines.clear();
	  for(i=1; i<WList[0].size(); i+=2)
	    {
	      Aux_FileName.push_back(WList[0][i]);
	      Aux_SkipLines.push_back((unsigned int)atoi(WList[0][i+1].c_str()));
	    } /* End For */
	  List_VALID_FileName.push_back(Aux_FileName);
	  List_VALID_SkipLines.push_back(Aux_SkipLines);
	} /* End Else If */
      else if (WList[0][0]=="INPUT")
	{
	  if ((WList[0].size()!=5)&&
	      (WList[0].size()!=6)&&
	      (WList[0].size()!=7))
	    {
	      cerr << "ERROR: INPUT usage" << endl;
	      cerr << "       INPUT Name Pos Tempo FileNumber GroupNb DontCut" << endl;
	      cerr << "Name: name of the input" << endl;
	      cerr << "Pos:  column number of the input in the data file" << endl;
	      cerr << "Tempo: number of delay samples" << endl;
	      cerr << "FileNumber: file number of which input belongs" << endl;
	      cerr << "GroupNb: number of the group to which the input belongs (for R2 computation)" << endl;
	      cerr << "DontCut: 0 = false 1 = true. True if you dont wan't any cut with respect to this variable" << endl;
	      exit(1);
	    } /* End If */

	  One_Model_Input.Name       = WList[0][1];
	  One_Model_Input.Pos        = atoi(WList[0][2].c_str());
	  One_Model_Input.Tempo      = atoi(WList[0][3].c_str());
	  One_Model_Input.FileNumber = atoi(WList[0][4].c_str());
	  One_Model_Input.GroupNb    = 0;
	  One_Model_Input.DontCut    = false;
	  One_Model_Input.Type       = Input;

	  if (WList[0].size()==6)
	    {
	      One_Model_Input.GroupNb = atoi(WList[0][5].c_str());
	    } /* End If */

	  if (WList[0].size()==7)
	    {
	      One_Model_Input.GroupNb = atoi(WList[0][5].c_str());
	      One_Model_Input.DontCut = (atoi(WList[0][6].c_str())==1);
	    } /* End If */

	  List_Model_Input.push_back(One_Model_Input);
	  
	  if (NbMaxTempo<One_Model_Input.Tempo) NbMaxTempo = One_Model_Input.Tempo;
	} /* End Else If */
      else if (WList[0][0]=="OUTPUT")
	{
	  if (WList[0].size()!=5)
	    {
	      cerr << "ERROR: OUTPUT usage" << endl;
	      cerr << "       OUTPUT Name Pos Tempo FileNumber" << endl;
	      cerr << "Name: name of the input" << endl;
	      cerr << "Pos:  column number of the output in the data file" << endl;
	      cerr << "Tempo: number of delay samples" << endl;
	      cerr << "FileNumber: file number of which output belongs" << endl;
	      exit(1);
	    } /* End If */

	  One_Model_Output.Name       = WList[0][1];
	  One_Model_Output.Pos        = atoi(WList[0][2].c_str());
	  One_Model_Output.Tempo      = atoi(WList[0][3].c_str());
	  One_Model_Output.FileNumber = atoi(WList[0][4].c_str());
	  One_Model_Output.GroupNb    = 0;
	  One_Model_Output.DontCut    = false;
	  One_Model_Output.Type       = Output;

	  List_Model_Output.push_back(One_Model_Output);

	  if (NbMaxTempo<One_Model_Output.Tempo) NbMaxTempo = One_Model_Output.Tempo;
	} /* End Else If */
      else if (WList[0][0]=="RETROINPUT")
	{
	  if ((WList[0].size()!=5)&&
	      (WList[0].size()!=6)&&
	      (WList[0].size()!=7))
	    {
	      cerr << "ERROR: RETROINPUT usage" << endl;
	      cerr << "       RETROINPUT Name Pos Tempo FileNumber GroupNb DontCut" << endl;
	      cerr << "Name: name of the input" << endl;
	      cerr << "Pos:  column number of the input in the data file" << endl;
	      cerr << "Tempo: number of delay samples" << endl;
	      cerr << "FileNumber: file number of which input belongs" << endl;
	      cerr << "GroupNb: number of the group to which the input belongs (for R2 computation)" << endl;
	      cerr << "DontCut: 0 = false 1 = true. True if you dont wan't any cut with respect to this variable" << endl;
	      exit(1);
	    } /* End If */

	  One_Model_Input.Name       = WList[0][1];
	  One_Model_Input.Pos        = atoi(WList[0][2].c_str());
	  One_Model_Input.Tempo      = atoi(WList[0][3].c_str());
	  One_Model_Input.FileNumber = atoi(WList[0][4].c_str());
	  One_Model_Input.GroupNb    = 0;
	  One_Model_Input.DontCut    = false;
	  One_Model_Input.Type       = RetroInput;

	  if (One_Model_Input.Tempo==0)
	    {
	      cerr << "ERROR: RETROINPUT - Tempo must be greater or equal to 1" << endl;
	      exit(1);
	    } /* End If */

	  if (WList[0].size()==6)
	    {
	      One_Model_Input.GroupNb = atoi(WList[0][5].c_str());
	    } /* End If */

	  if (WList[0].size()==7)
	    {
	      One_Model_Input.GroupNb = atoi(WList[0][5].c_str());
	      One_Model_Input.DontCut = (atoi(WList[0][6].c_str())==1);
	    } /* End If */

	  List_Model_Input.push_back(One_Model_Input);

	  List_Retro_Input.push_back(List_Model_Input.size()-1);
	} /* End Else If */
      else if (WList[0][0]=="WEIGHT")
	{
	  if (WList[0].size()!=4)
	    {
	      cerr << "ERROR: WEIGHT usage" << endl;
	      cerr << "       WEIGHT Name Pos FileNumber" << endl;
	      cerr << "Name: name of the input" << endl;
	      cerr << "Pos:  column number of the input in the data file" << endl;
	      cerr << "FileNumber: file number of which input belongs" << endl;
	      exit(1);
	    } /* End If */

	  Weight_Input.Name       = WList[0][1];
	  Weight_Input.Pos        = atoi(WList[0][2].c_str());
	  Weight_Input.FileNumber = atoi(WList[0][3].c_str());
	  Weight_Input.Tempo      = 0;
	  Weight_Input.GroupNb    = 0;
	  Weight_Input.DontCut    = false;
	  Weight_Input.Type       = Weight;
	} /* End Else If */
      else if (WList[0][0]=="OPTIMIZESIGMA")
	{
	  UseOptimizeSigma = true;
	  if (WList[0].size()==1)
	    {
	      SigmaMin  = 0.1;
	      SigmaStep = 0.01;
	      SigmaMax  = 10.0*SigmaStep;
	    } /* End If */
	  else if (WList[0].size()==2)
	    {
	      SigmaMin  = atof(WList[0][1].c_str());
	      SigmaStep = 0.01;
	      SigmaMax  = SigmaMin + 10.0*SigmaStep;
	    } /* End If */
	  else if (WList[0].size()==3)
	    {
	      SigmaMin  = atof(WList[0][1].c_str());
	      SigmaStep = atof(WList[0][2].c_str());
	      SigmaMax  = SigmaMin + 10.0*SigmaStep;
	    } /* End If */
	  else if (WList[0].size()==4)
	    {
	      SigmaMin  = atof(WList[0][1].c_str());
	      SigmaStep = atof(WList[0][2].c_str());
	      SigmaMax  = atof(WList[0][3].c_str());
	    } /* End If */
	  else
	    {
	      cerr << "ERROR: OPTIMIZESIGMA usage" << endl;
	      cerr << "       OPTIMIZESIGMA SigmaMin SigmaStep SigmaMax" << endl;
	      exit(1);
	    } /* End Else */
	} /* End Else If */
      else if (WList[0][0]=="NBMAXPART")
	{
	  NbMaxPart = atoi(WList[0][1].c_str());
	} /* End Else If */
      else if (WList[0][0]=="CUT")
	{
	  if (WList[0].size()!=6)
	    {
	      cerr << "ERROR: CUT usage" << endl;
	      cerr << "       CUT NameOfPart DimensionNo CutValue LowPartName UpperPartName" << endl;
	      cerr << "NameOfPart : name of the partition to cut" << endl;
	      cerr << "DimensionNo: index of the dimension to cut (starting from 0)" << endl;
	      cerr << "CutValue   : where to put the cut" << endl;
	      cerr << "LowPartName: name of the partition name (partition with the lowest values)" << endl;
	      cerr << "UpPartName : name of the partition name (partition with the largest values)" << endl;
	      exit(1);
	    } /* End If */

	  Aux_PartToCut.PartName           = WList[0][1];
	  Aux_PartToCut.Dimension          = atoi(WList[0][2].c_str());
	  Aux_PartToCut.CutPosition        = atof(WList[0][3].c_str());
	  Aux_PartToCut.PartName_LowerPart = WList[0][4];
	  Aux_PartToCut.PartName_UpperPart = WList[0][5];
	  
	  List_PartToCut.push_back(Aux_PartToCut);
	} /* End Else If */
      else if (WList[0][0]=="SETMODELCOEFF")
	{
	  if (WList[0].size()!=4)
	    {
	      cerr << "ERROR: SETMODELCOEFF usage" << endl;
	      cerr << "       SETMODELCOEFF NameOfPart DimensionNo CoeffValue" << endl;
	      cerr << "NameOfPart : name of the partition to cut" << endl;
	      cerr << "DimensionNo: index of the dimension to cut (starting from 0). index -1 corresponds to the" << endl;
	      cerr << "             constant coefficient." << endl;
	      cerr << "CoeffValue : the value of the coefficient" << endl;
	      exit(1);
	    } /* End If */

	  Aux_ModelCoeff.PartName   = WList[0][1];
	  Aux_ModelCoeff.Dimension  = atoi(WList[0][2].c_str());
	  Aux_ModelCoeff.Coeff      = atof(WList[0][3].c_str());
	  
	  List_ModelCoeff.push_back(Aux_ModelCoeff);
	} /* End Else If */
      else if (WList[0][0]=="SAVESAMPLEDDATA")
	{
	  SampledDataFile = WList[0][1];
	  UseSampledData = true;
	  if (atoi(WList[0][2].c_str())==1) SampledData_Exit = true;
	} /* End Else If */
      else if (WList[0][0]=="SKIPLEARN")
	{
	  if (WList.size()!=3)
	    {
	      cerr << "ERROR: SKIPLEARN usage" << endl;
	      cerr << "       SKIPLEARN BeginLine EndLine" << endl;
	      cerr << "BeginLine: Index of the beginning of the sequence to remove" << endl;
	      cerr << "EndLine:   Index of the end of the sequence to remove" << endl;
	      exit(1);
	    } /* End If */
	  Use_SkipLearn = true;
	  Aux_DataSkip.Begin_Part = atoi(WList[0][1].c_str());
	  Aux_DataSkip.End_Part   = atoi(WList[0][2].c_str());
	  List_SkipLearn.push_back(Aux_DataSkip);
	} /* End Else If */
      else if (WList[0][0]=="SKIPVALID")
	{
	  if ((WList.size()<3)||(WList.size()>4))
	    {
	      cerr << "ERROR: SKIPVALID usage" << endl;
	      cerr << "       SKIPVALID NumOfValidFile BeginLine EndLine" << endl;
	      cerr << "NumOfValidFile (optional): Number of the valid file on which to remove a sequence of data" << endl;
	      cerr << "BeginLine: Index of the beginning of the sequence to remove" << endl;
	      cerr << "EndLine:   Index of the end of the sequence to remove" << endl;
	      exit(1);
	    } /* End If */
	  Use_SkipValid = true;
	  if (WList.size()==4)
	    {
	      Aux_DataSkip.NumOfFile  = atoi(WList[0][1].c_str());
	      Aux_DataSkip.Begin_Part = atoi(WList[0][2].c_str());
	      Aux_DataSkip.End_Part   = atoi(WList[0][3].c_str());
	    } /* End If */
	  if (WList.size()==3)
	    {
	      Aux_DataSkip.NumOfFile  = 0;
	      Aux_DataSkip.Begin_Part = atoi(WList[0][2].c_str());
	      Aux_DataSkip.End_Part   = atoi(WList[0][3].c_str());
	    } /* End If */
	  List_SkipValid.push_back(Aux_DataSkip);
	} /* End Else If */
      else if (WList[0][0]=="EXPORTC")
	{
	  ExportC = true;
	} /* End Else If */
      else if (WList[0][0]=="EXPORTCPP")
	{
	  ExportCpp = true;
	} /* End Else If */
      else if (WList[0][0]=="EXPORTMATLAB")
	{
	  ExportMatlab = true;
	} /* End Else If */
      else if (WList[0][0]=="DEXPORTC")
	{
	  D_ExportC = true;
	} /* End Else If */
      else if (WList[0][0]=="DEXPORTC2")
	{
	  D_ExportC2 = true;
	} /* End Else If */
      else if (WList[0][0]=="DEXPORTCPP")
	{
	  D_ExportCpp = true;
	} /* End Else If */
      else if (WList[0][0]=="DEXPORTCPP2")
	{
	  D_ExportCpp2 = true;
	} /* End Else If */
      else if (WList[0][0]=="DEXPORTMATLAB")
	{
	  D_ExportMatlab = true;
	} /* End Else If */
      else if (WList[0][0]=="DEXPORTMATLAB2")
	{
	  D_ExportMatlab2 = true;
	} /* End Else If */
      else if (WList[0][0]=="DDEXPORTC")
	{
	  DD_ExportC = true;
	} /* End Else If */
      else if (WList[0][0]=="DDEXPORTC2")
	{
	  DD_ExportC2 = true;
	} /* End Else If */
      else if (WList[0][0]=="DDEXPORTCPP")
	{
	  DD_ExportCpp = true;
	} /* End Else If */
      else if (WList[0][0]=="DDEXPORTCPP2")
	{
	  DD_ExportCpp2 = true;
	} /* End Else If */
      else if (WList[0][0]=="DDEXPORTMATLAB")
	{
	  DD_ExportMatlab = true;
	} /* End Else If */
      else if (WList[0][0]=="DDEXPORTMATLAB2")
	{
	  DD_ExportMatlab2 = true;
	} /* End Else If */
      else if (WList[0][0]=="SAVEMODEL")
	{
	  Save_Model_Filename = WList[0][1];
	  Save_Model     = true;
	} /* End Else If */
      else if (WList[0][0]=="LOADMODEL")
	{
	  Load_Model_Filename = WList[0][1];
	  Load_Model     = true;
	} /* End Else If */
      else if (WList[0][0]=="SETWD")
	{
	  Working_Dir = WList[0][1].c_str();
	} /* End Else If */
      else if (WList[0][0]=="CROSSVALID")
	{
	  if (WList[0].size()!=3)
	    {
	      cerr << "ERROR: CROSSVALID usage" << endl;
	      cerr << "       CUT NameOfFile NbSplit" << endl;
	      cerr << "NameOfFile : name of the file where we store the results" << endl;
	      cerr << "NbSplit    : Nb of division of the data file to perform the cross validation" << endl;
	      exit(1);
	    } /* End If */
	  CrossValid_Filename = WList[0][1].c_str();
	  CrossValid_NbFiles  = atoi(WList[0][2].c_str());
	  CrossValid          = true;
	} /* End Else If */
      else if (WList[0][0]=="RANDSEED")
	{
	  RandSeed = atoi(WList[0][1].c_str());
	} /* End Else If */
      else if (WList[0][0]=="USELASSO")
	{
	  UseLasso = true;
	} /* End Else If */
      else if (WList[0][0]=="RESIDUALTYPE")
	{
	  ResidualType = atoi(WList[0][1].c_str());
	} /* End Else If */
      else if (WList[0][0]=="RESIDUALFILENAME")
	{
	  Res_Filename = WList[0][1];
	} /* End Else If */
      else if (WList[0][0]=="R2FILENAME")
	{
	  R2_Filename = WList[0][1];
	} /* End Else If */
      else if (WList[0][0]=="PARTLIST")
	{
	  Part_Filename   = WList[0][1];
	  UsePartFilename = true;
	} /* End Else If */
      else if (WList[0][0]=="LEAVEONEOUT")
	{
	  LOO_Filename = WList[0][1];
	  Use_LOO      = true;
	} /* End Else If */
      else if (WList[0][0]=="MEMBERSHIPTHRESHOLD")
	{
	  MembershipThreshold = atof(WList[0][1].c_str());
	} /* End Else If */
      else if (WList[0][0]=="USEDELTA")
	{
	  UseDelta = true;
	} /* End Else If */
      else if (WList[0][0]=="LEARNMETH")
	{
	  BS_Classic   = false;
	  BS_Light     = false;
	  BS_NoMeth    = false;

	  if (WList[0][1]=="NO")
	    {
	      BS_NoMeth = true;
	    } /* End If */
	  if (WList[0][1]=="CLASSIC")
	    {
	      BS_Classic   = true;
	      BS_Light     = false;
	    } /* End If */
	  if (WList[0][1]=="BOOTSTRAPLIGHT")
	    {
	      BS_Classic = false;
	      BS_Light   = true;

	      if (WList[0].size()==2)
		{
		  BS_Proportion   = 0.2;
		  NbPartToExplore = -1;
		} /* End If */
	      else if (WList[0].size()==3)
		{
		  BS_Proportion   = atof(WList[0][2].c_str());
		  NbPartToExplore = -1;
		} /* End Else If */
	      else if (WList[0].size()==4)
		{
		  BS_Proportion   = atof(WList[0][2].c_str());
		  NbPartToExplore = atoi(WList[0][3].c_str());
		} /* End Else */
	      else
		{
		  cerr << "Error: BOOTSTRAPLIGHT usage:" << endl;
		  cerr << "       BOOTSTRAPLIGHT ValidProp NbPart" << endl;
		  cerr << "ValidProp: Quantity of validation points in percent of the total number of learning points" << endl;
		  cerr << "NbPart   : Number of partitions to generate" << endl;
		  exit(1);
		} /* End Else */
	    } /* End If */ 
	  if (WList[0][1]=="RECURSIVE")
	    {
	      Learn_Recur = true;
	      if (WList[0].size()==2)
		{
		  Rec_NbPart_Ph1 = 10;
		  Rec_NbPart_Ph2 = 10;
		} /* End If */
	      else if (WList[0].size()==3)
		{
		  Rec_NbPart_Ph1 = atoi(WList[0][2].c_str());
		  Rec_NbPart_Ph2 = Rec_NbPart_Ph1;
		} /* End Else If */
	      else if (WList[0].size()==4)
		{
		  Rec_NbPart_Ph1 = atoi(WList[0][2].c_str());
		  Rec_NbPart_Ph2 = atoi(WList[0][3].c_str());
		} /* End Else If */
	      else
		{
		  cerr << "Error: RECURSIVE usage:" << endl;
		  cerr << "       RECURSIVE [NbPart Phase1] [NbPart Phase 2]" << endl;
		  cerr << " Both arguments are optional. If no argument found then both parameters are set to 10." << endl;
		  cerr << " If one argument found, both parameters are set to";
		  cerr << " the same value as the first argument." << endl;
		  exit(-1);
		} /* End Else */
	    } /* End If */
	  if (WList[0][1]=="UPDATE")
	    {
	      if (WList[0].size()==3)
		{
		  UseUpdate = true;
		  Nb_Iter    = atoi(WList[0][2].c_str());
		} /* End If */
	      else
		{
		  cerr << "Error: UPDATE usage:" << endl;
		  cerr << "       UPDATE NbIter" << endl;
		  cerr << "NbIter is the number of partitions during which we perform updating";
		  cerr << "of the recursive model" << endl;
		  exit(-1);
		} /* End Else */
	    } /* End If */
	} /* End Else If */
      else if (WList[0][0]=="COMPUTER2")
	{
	  BS_ComputeR2 = true;
	} /* End If */
      else if (WList[0][0]=="COMPUTEC1")
	{
	  Compute_C1 = true;
	} /* End If */
      else if (WList[0][0]=="COMPUTEC2")
	{
	  Compute_C2 = true;
	} /* End If */
      else if (WList[0][0]=="VALIDSTART")
	{
	  Valid_Start = atoi(WList[0][1].c_str());
	} /* End If */
      else if (WList[0][0]=="VALIDEND")
	{
	  Valid_End = atoi(WList[0][1].c_str());
	} /* End If */
      else if (WList[0][0]=="LEARNSTART")
	{
	  Learn_Start = atoi(WList[0][1].c_str());
	} /* End If */
      else if (WList[0][0]=="LEARNEND")
	{
	  Learn_End = atoi(WList[0][1].c_str());
	} /* End If */
      else if (WList[0][0]=="SHIFTINPUTS")
	{
	  ShiftInput = atoi(WList[0][1].c_str());
	} /* End If */
      else if (WList[0][0]=="SHIFTOUTPUTS")
	{
	  ShiftOutput = atoi(WList[0][1].c_str());
	} /* End If */
      else if (WList[0][0]=="USETRANSFORMLOG")
	{
	  UseTransformLog      = true;
	  UseTransformLogAlpha = atof(WList[0][1].c_str());
	  if (WList[0].size()==3) UseTransformLogEps = atof(WList[0][2].c_str());
	  else                    UseTransformLogEps = 0.0;
	} /* End Else If */
      else if (WList[0][0]=="RESULTDIR")
	{
	  ResultDir = WList[0][1];
	} /* End Else If */
      else if (WList[0][0]=="CUSTOMPARAM")
	{
	  CustomParam = true;
	  // We concatenate the arguments. So, we will have a list of arguments separated by spaces
	  CustParamStr = "";
	  for(i=1; i<WList[0].size(); i++)
	    {
	      CustParamStr += WList[0][i];
	    } /* End For */
	} /* End Else If */
      else if (WList[0][0]=="CUSTOMRETURNBEFORE")
	{
	  CustomRet_Before = true;
	  // We concatenate the arguments. So, we will have a list of arguments separated by spaces
	  CustRetBefStr = "";
	  for(i=1; i<WList[0].size(); i++)
	    {
	      CustRetBefStr += WList[0][i];
	    } /* End For */
	} /* End Else If */
      else if (WList[0][0]=="CUSTOMRETURNAFTER")
	{
	  CustomRet_After = true;
	  // We concatenate the arguments. So, we will have a list of arguments separated by spaces
	  CustRetAfterStr = "";
	  for(i=1; i<WList[0].size(); i++)
	    {
	      CustRetAfterStr += WList[0][i];
	    } /* End For */
	} /* End Else If */
      else if (WList[0][0]=="MAVERAGE")
	{
	  MAverage = atoi(WList[0][1].c_str());
	} /* End Else If */
      else if (WList[0][0]=="NORMALIZEINOUT")
	{
	  UseNormalize = true;
	  Norm_Filename = WList[0][1];
	} /* End Else If */
      else if (WList[0][0]=="PRINTPARTITIONS")
	{
	  PrintPart = true;
	} /* End Else If */
      else if (WList[0][0]=="PRINTDISTRIBOFCUT")
	{
	  if (WList[0].size()==3)
	    {
	      cerr << "Error: PRINTDISTRIBOFCUT usage:" << endl;
	      cerr << "       PRINTDISTRIBOFCUT Filename" << endl;
	      cerr << "Filename is the name of the file where the cuts will be wrote" << endl;
	      exit(-1);
	    } /* End If */
	  PrintDistribOfCut = true;
	  DistribOfCutName  = WList[0][1];
	} /* End Else If */
      else if (WList[0][0]=="USEDISTRIBUTEDCUTTING")
	{
	  UseDistribCutting  = true;
	  UseUniformCutting  = false;
	} /* End Else If */
      else if (WList[0][0]=="USEUNIFORMCUTTING")
	{
	  UseDistribCutting  = false;
	  UseUniformCutting  = true;
	} /* End Else If */
      else if (WList[0][0]=="DISPLAY")
	{
	  Display = true;
	} /* End Else If */
      else if (WList[0][0]=="FILTERSTEP")
	{
	  Filter_Step = atoi(WList[0][1].c_str());
	} /* End Else If */
      else if (WList[0][0]=="SELECTSTEP")
	{
	  Select_Step = atoi(WList[0][1].c_str());
	} /* End Else If */
      else if (WList[0][0]=="RESGAPPERC")
	{
	  ResiduGapPercentage = atof(WList[0][1].c_str());
	} /* End Else If */
      else if (WList[0][0]=="ECHANTSTAT")
	{
	  Echant_Stat_NbBar    = atoi(WList[0][1].c_str());
	  Echant_Stat_NbPoints = atoi(WList[0][2].c_str());
	} /* End Else If */
      else if (WList[0][0]=="ECHANTDELTA")
	{
	  Echant_Delta_NbPoints = atoi(WList[0][1].c_str());
	  Echant_Delta_NbTimes  = atoi(WList[0][2].c_str());
	} /* End Else If */
      else if (WList[0][0]=="ECHANTADDMINMAX")
	{
	  Echant_AddMinMax = true;
	} /* End Else If */
      else if (WList[0][0]=="ECHANTDYN")
	{
	  Echant_Dyn_NbBar    = atoi(WList[0][1].c_str());
	  Echant_Dyn_NbPoints = atoi(WList[0][2].c_str());
	} /* End Else If */
      else if (WList[0][0]=="ECHANTRAND")
	{
	  Echant_Rand_NbPoints = atoi(WList[0][1].c_str());
	} /* End Else If */
      else if (WList[0][0]=="ECHANTGREEDY")
	{
	  Echant_Greedy_NbPoints = atoi(WList[0][1].c_str());
	} /* End Else If */
      else if (WList[0][0]=="ECHANTGREEDYFAST")
	{
	  Echant_Greedy_Fast_NbPointsAtTheEnd = atoi(WList[0][1].c_str());
	  Echant_Greedy_Fast_NbPointsToSelect = atoi(WList[0][2].c_str());
	} /* End Else If */
      else if (WList[0][0]=="ECHANTGREEDYADDMEASURE")
	{
	  Echant_Greedy_Add_Measure_NbPoints = atoi(WList[0][1].c_str());
	} /* End Else If */
      else if (WList[0][0]=="ECHANTGREEDYFASTADDMEASURE")
	{
	  Echant_Greedy_Fast_Add_Measure_NbPointsAtTheEnd = atoi(WList[0][1].c_str());
	  Echant_Greedy_Fast_Add_Measure_NbPointsToSelect = atoi(WList[0][2].c_str());
	} /* End Else If */
      else if (WList[0][0]=="ECHANTGREEDYPARTIALADDMEASURE")
	{
	  Echant_Greedy_Partial_NbPtsToSelect   = atoi(WList[0][1].c_str());
	  Echant_Greedy_Partial_NbTimesToSelect = atoi(WList[0][2].c_str());
	} /* End Else If */      
      else if (WList[0][0]=="ECHANTGREEDYPARTIALDELTA")
	{
	  Echant_PartGreedDelta_NbPoints = atoi(WList[0][1].c_str());
	  Echant_PartGreedDelta_NbTimes  = atoi(WList[0][2].c_str());
	} /* End Else If */      
      else if (WList[0][0]=="ECHANTOUTPUTMIN")
	{
	  Echant_Output_NbPts = atoi(WList[0][1].c_str());
	  Echant_Output_UseMin = true;
	} /* End Else If */
      else if (WList[0][0]=="ECHANTOUTPUTMAX")
	{
	  Echant_Output_NbPts = atoi(WList[0][1].c_str());
	  Echant_Output_UseMin = false;
	} /* End Else If */
      else if (WList[0][0]=="ECHANTINPUTADD")
	{
	  Aux_Entry.Min     = atof(WList[0][1].c_str());
	  Aux_Entry.Max     = atof(WList[0][2].c_str());
	  Aux_Entry.InputNb = atoi(WList[0][3].c_str());
	  Aux_Entry.Add     = true;

	  List_Filter.push_back(Aux_Entry);
	} /* End Else If */
      else if (WList[0][0]=="ECHANTINPUTREMOVE")
	{
	  Aux_Entry.Min     = atof(WList[0][1].c_str());
	  Aux_Entry.Max     = atof(WList[0][2].c_str());
	  Aux_Entry.InputNb = atoi(WList[0][3].c_str());
	  Aux_Entry.Add     = false;

	  List_Filter.push_back(Aux_Entry);
	} /* End Else If */
      else if (WList[0][0]=="SIGMA")
	{
	  Sigma = atof(WList[0][1].c_str());
	} /* End Else If */
      else if (WList[0][0]=="NBCUT")
	{
	  NbCut = atoi(WList[0][1].c_str());
	} /* End Else If */
      else if (WList[0][0]=="MEASURE")
	{
	  Measure = atoi(WList[0][1].c_str());
	} /* End Else If */
      else if (strlen(CfgBuffer)==0)
	{
	} /* End Else If */
      else
	{
	  cerr << "Lecture du fichier " << argv[1] << " : mauvais mot clef : Buffer = " << CfgBuffer << endl;
	  exit(1);
	} /* End Else */

      WList.resize(1);
      DList.resize(1);
    } /* End While */
  
  ConfigFile->close();
  
  delete ConfigFile;

  delete [] CfgBuffer;

  // Initialisation of the random generator
  srand(RandSeed);

  // Temporisation are measure starting from 1 in the configuration file
  if (NbMaxTempo!=0) NbMaxTempo--;

  cerr << "End of the reading of the configuration file" << endl;

  if (ModelName=="") ModelName = "Default";
  if (List_TRAIN_FileName.size()==0)
    {
      cerr << "No learning data files" << endl;
      exit(1);
    } /* End If */
  if (List_VALID_FileName.size()==0)
    {
      cerr << "No validation data files" << endl;
      Validation_File = false;
    } /* End If */
  else
    {
      Validation_File = true;
    } /* End Else */
  if (NbMaxPart==-1)
    {
      NbMaxPart = NBMAXPARTITIONS;
    } /* End If */
  if (ResiduGapPercentage==-1.0)
    {
      ResiduGapPercentage =  RESIDUGAPPERCENTAGE;
    } /* End If */
  if (Sigma==-1.0)
    {
      Sigma = LL_SIGMA;
    } /* End If */
  if (NbCut==-1)
    {
      NbCut = LL_NBCUT;
    } /* End If */
  if (Measure==-1.0)
    {
      Measure = 0;
    } /* End If */
  if (NbPartToExplore==-1)
    {
      NbPartToExplore=NbMaxPart;
    } /* End If */
  if (ResultDir.c_str()[ResultDir.length()-1]!='/') ResultDir = ResultDir + "/";

  if (Display)
    {
      cout << "Model : " << ModelName << endl;
      
      for(i=0; i<List_TRAIN_FileName.size(); i++)
	{
	  cout << "Learning data " << i;
	  for(j=0; j<List_TRAIN_FileName[i].size(); j++)
	    {
	      cout << " " << List_TRAIN_FileName[i][j];
	    } /* End For */
	  cout << endl;
	} /* End For */
      
      for(i=0; i<List_VALID_FileName.size(); i++)
	{
	  cout << "Validation data " << i;
	  for(j=0; j<List_VALID_FileName[i].size(); j++)
	    {
	      cout << " " << List_VALID_FileName[i][j];
	    } /* End For */
	  cout << endl;
	} /* End For */

      for(i=0; i<List_Model_Input.size(); i++)
	{
	  cout << "Input " << i;
	  cout << " : Name = "      << List_Model_Input[i].Name;
	  cout << " Pos = "         << List_Model_Input[i].Pos;
	  cout << " Tempo = "       << List_Model_Input[i].Tempo;
	  cout << " FileNumber = "  << List_Model_Input[i].FileNumber;
	  cout << endl;
	} /* End For */
      
      for(i=0; i<List_Model_Output.size(); i++)
	{
	  cout << "Output " << i;
	  cout << " : Name = "      << List_Model_Output[i].Name;
	  cout << " Pos = "         << List_Model_Output[i].Pos;
	  cout << " Tempo = "       << List_Model_Output[i].Tempo;
	  cout << " FileNumber = "  << List_Model_Output[i].FileNumber;
	  cout << endl;
	} /* End For */
      
      cout << "NbMaxPart           = " << NbMaxPart           << endl;
      cout << "ResiduGapPercentage = " << ResiduGapPercentage << endl;
      cout << "Sigma               = " << Sigma               << endl;
      cout << "NbCut               = " << NbCut               << endl;
      cout << "Measure             = " << Measure             << endl;
    } /* End If */

  // We change the working directory

  chdir(Working_Dir.c_str());

  // Opening the files

  cout << "Learning - Reading the data - Part 1 - We don't generate the offsets" << endl;

  Offset_For_Files.resize(List_TRAIN_FileName.size(), 0);

  for(n=0; n<List_TRAIN_FileName.size(); n++)
    {
      // We allocate as many file descriptors as there file name on a TRAIN line
      DataFile = new ifstream[List_TRAIN_FileName[n].size()];
      
      // We allocate a buffer by file descriptor
      typedef char * PtrChar;
      
      Buffer = new PtrChar[List_TRAIN_FileName[n].size()];
      
      // We open all the files a we initialize the buffers
      for(m=0; m<List_TRAIN_FileName[n].size(); m++)
	{
	  cout << "Opening the file " << List_TRAIN_FileName[n][m] << endl;
	  
	  DataFile[m].open(List_TRAIN_FileName[n][m].c_str(), ios_base::in | ios_base::binary);
	  
	  if (!DataFile[m].is_open())
	    {
	      cerr << "Learning - Can't open file " << List_TRAIN_FileName[n][m];
	      cerr << "." << endl;

	      // We change the working directory
	      chdir(Current_Dir.c_str());

	      exit(1);
	    } /* End If */
	  
	  Buffer[m] = new char[BUFFERSIZE];
	} /* End For */
      
      Result = true;
      WList.resize(List_TRAIN_FileName[n].size());
      
      // We skip the lines which contains the labels of the columns
      for(m=0; m<List_TRAIN_FileName[n].size(); m++)
	{
	  for(i=0; i<List_TRAIN_SkipLines[n][m]; i++)
	    {
	      DataFile[m].getline(Buffer[m], BUFFERSIZE, '\n');
	    } /* End For */
	} /* End For */
      
      LineCount = -1;
      
      while(Result)
	{
	  // We store a line of each files in each buffer
	  NbColAux = 0;
	  for(m=0; m<List_TRAIN_FileName[n].size(); m++)
	    {
	      Result   = Result & (bool)DataFile[m].getline(Buffer[m], BUFFERSIZE, '\n');
	      vector<string> tmpLine = tokenize(Buffer[m], is_space());
	      if (tmpLine.size()!=0) WList[m] = tmpLine; // Skip empty lines
	      
	      if (!NbColSet) NbCol    += WList[m].size();
	      else           NbColAux += WList[m].size();

	      if      ((n==0)&&(m==0)) Offset_For_Files[m] = WList[m].size();
	      else if ((n==0)&&(m!=0)) Offset_For_Files[m] = Offset_For_Files[m-1] + WList[m].size();
	    } /* End For */
	  
	  if (!NbColSet) 
	    {
	      NbColSet = true;
	      NbColAux = NbCol;
	    } /* End If */
	  
	  if ((Result)&&(NbCol!=NbColAux))
	    {
	      cerr << "Error while reading files " << List_TRAIN_FileName[n][0] << " ... " << endl;
	      cerr << "Error around line " << LineCount << endl;
	      cerr << "Number of column : " << NbCol << endl;
	      cerr << "Number of column read : " << NbColAux << endl;

	      // We change the working directory
	      chdir(Current_Dir.c_str());

	      exit(-1);
	    } /* End If */
	  
	  LineCount++;
	  
	  // If we arrive at the end of a file, we stop
	  // We Post-process InputData and MeasureData: we remove the NbMaxTemp first lines
	  // to avoid outlier data
	  
	  // We get all the data
	  DList.resize(0);
	  
	  for(j=0; j<WList.size(); j++)
	    {
	      for(k=0; k<WList[j].size(); k++)
		{
		  DList.push_back(atof(WList[j][k].c_str()));
		} /* End For */
	    } /* End For */
	  
	  InputData.push_back(DList);
	} /* End While */
      
      for(m=0; m<List_TRAIN_FileName[n].size(); m++)
	{
	  DataFile[m].close();
	  if (Buffer[m]) delete [] Buffer[m];
	} /* End For */
      
      if (DataFile) delete [] DataFile;
      if (Buffer)   delete [] Buffer;
    } /* End For n */

  // Modification of the position of the columns
  for(i=1; i<List_Model_Input.size(); i++)
    {
      if (List_Model_Input[i].FileNumber!=0)
	List_Model_Input[i].Pos += Offset_For_Files[List_Model_Input[i].FileNumber-1];
    } /* End For */

  for(i=0; i<List_Model_Input.size(); i++)
    {
      List_Model_Input[i].FileNumber = 0;
    } /* End For */

  for(i=1; i<List_Model_Output.size(); i++)
    {
      if (List_Model_Output[i].FileNumber!=0)
	List_Model_Output[i].Pos += Offset_For_Files[List_Model_Output[i].FileNumber-1];
    } /* End For */

  for(i=0; i<List_Model_Output.size(); i++)
    {
      List_Model_Output[i].FileNumber = 0;
    } /* End For */

  Weight_Input.Pos        += Offset_For_Files[List_Model_Input[0].FileNumber];
  Weight_Input.FileNumber  = 0;

  // Generation of a fake MeasureData set
  MeasureData.resize(InputData.size(), vector<float>(List_Model_Output.size(), 0.0));

  if (Filter_Step>0)
    {
      cout << "Learning - phase 1 - Filtering the data" << endl;

      Filter_Filter(InputData, MeasureData, Filter_Step);
    } /* End If */

  if (MAverage>0)
    {
      cout << "Learning - phase 1 - Moving average applied on the data" << endl;
      
      Filter_MovingAverage(InputData, MeasureData, MAverage);
    } /* End If */

  if (Select_Step>0)
    {
      cout << "Learning - phase 1 - Select of a data line out of " << Select_Step << endl;

      Filter_Select(InputData, MeasureData, Select_Step);
    } /* End If */

  cout << "Learning - Reading the data - Part 2 - generation of the offsets" << endl;

  // Opening the data files

  InputMin.resize(List_Model_Input.size());
  InputMax.resize(List_Model_Input.size());
  
  for(i=0; i<List_Model_Input.size(); i++)
    {
      InputMin[i]   = numeric_limits<float>::max();
      InputMax[i]   = - numeric_limits<float>::max();
    } /* End For */
  
  MeasureMin.resize(List_Model_Output.size());
  MeasureMax.resize(List_Model_Output.size());
  
  for(i=0; i<List_Model_Output.size(); i++)
    {
      MeasureMin[i] = numeric_limits<float>::max();
      MeasureMax[i] = - numeric_limits<float>::max();
    } /* End For */
  
  // Sizing of the input temporizations
  
  TempoInput.resize(List_Model_Input.size());
  
  for(i=0; i<List_Model_Input.size(); i++)
    {
      TempoInput[i].resize(List_Model_Input[i].Tempo + ShiftInput, 0.0);
    } /* End For */

  // Sizing of the output temporizations
  
  TempoOutput.resize(List_Model_Output.size());
  
  for(i=0; i<List_Model_Output.size(); i++)
    {
      TempoOutput[i].resize(List_Model_Output[i].Tempo + ShiftOutput, 0.0);
    } /* End For */

  AuxData.resize(InputData.size(), vector<float>(InputData[0].size(), 0.0));

  AuxData = InputData;

  WeightData.resize(0);
  InputData.resize(0);
  MeasureData.resize(0);

  for(i=0; i<AuxData.size(); i++)
    {
      if (i>=NbMaxTempo)
	{
	  // We get back the inputs
	  DList.resize(0);
	  
	  for(j=0; j<List_Model_Input.size(); j++)
	    {
	      DList.push_back(Temporize(TempoInput[j], AuxData[i][List_Model_Input[j].Pos]));
	    } /* End For */
	  
	  InputData.push_back(DList);
	  
	  for(j=0; j<DList.size(); j++)
	    {
	      if (InputMin[j]>DList[j]) InputMin[j] = DList[j];
	      if (InputMax[j]<DList[j]) InputMax[j] = DList[j];
	    } /* End For */
	  
	  // We get back the outputs
	  DList.resize(0);

	  for(j=0; j<List_Model_Output.size(); j++)
	    {
	      DList.push_back(Temporize(TempoOutput[j], AuxData[i][List_Model_Output[j].Pos]));
	    } /* End For */
	      
	  MeasureData.push_back(DList);
	  
	  for(j=0; j<DList.size(); j++)
	    {
	      if (MeasureMin[j]>DList[j]) MeasureMin[j] = DList[j];
	      if (MeasureMax[j]<DList[j]) MeasureMax[j] = DList[j];
	    } /* End For */

	  // We get back the weights
	  DList.resize(0);

	  for(j=0; j<WeightData.size(); j++)
	    {
	      WeightData.push_back(AuxData[i][Weight_Input.Pos]);
	    } /* End For */
	} /* End If */
      else if (i<NbMaxTempo)
	{
	  // We put the input data in the temporization line
	  for(j=0; j<List_Model_Input.size(); j++)
	    {
	      Temporize(TempoInput[j], AuxData[i][List_Model_Input[j].Pos]);
	    } /* End For */
	  
	  // We put the output data in the temporization line
	  for(j=0; j<List_Model_Output.size(); j++)
	    {
	      Temporize(TempoOutput[j], AuxData[i][List_Model_Output[j].Pos]);
	    } /* End For */
	} /* End Else */
    } /* End For */

  AuxData.clear();

  if (UseDelta)
    {
      // Processing inputs
      Delta_Entry AuxDelta;
      bool        InputInListDelta = false;
      int         AuxValue;

      List_Delta.resize(0);

      for(i=0; i<List_Model_Input.size()-1; i++)
	{
	  // Is the input i already in the List_Delta list ?

	  InputInListDelta = false;

	  for(j=0; j<List_Delta.size(); j++)
	    {
	      for(k=0; k<List_Delta[j].ListOfInputs.size(); k++)
		{
		  InputInListDelta = InputInListDelta && (i==List_Delta[j].ListOfInputs[k]);
		} /* End For */
	    } /* End For */

	  if (InputInListDelta) continue;

	  for(j=i+1; j<List_Model_Input.size(); j++)
	    {
	      if (List_Model_Input[i].Name==List_Model_Input[j].Name)
		{
		  AuxDelta.ListOfInputs.push_back(j);
		  AuxDelta.ListOfTempo.push_back(List_Model_Input[j].Tempo);
		} /* End If */
	    } /* End For */
	} /* End For */

      // We order each List_Delta entry with respect to the Tempo

      for(i=0; i<List_Delta.size(); i++)
	{
	  for(j=0; j<List_Delta[i].ListOfTempo.size()-1; j++)
	    {
	      for(k=j; k<List_Delta[i].ListOfTempo.size(); k++)
		{
		  if (List_Delta[i].ListOfTempo[j]>List_Delta[i].ListOfTempo[k])
		    {
		      AuxValue = List_Delta[i].ListOfTempo[j];
		      List_Delta[i].ListOfTempo[j] = List_Delta[i].ListOfTempo[k];
		      List_Delta[i].ListOfTempo[k] = AuxValue;

		      AuxValue = List_Delta[i].ListOfInputs[j];
		      List_Delta[i].ListOfInputs[j] = List_Delta[i].ListOfInputs[k];
		      List_Delta[i].ListOfInputs[k] = AuxValue;
		    } /* End If */
		} /* End For */
	    } /* End For */
	} /* End For */

      // InputData
      // We replace E(0) with 0.5*(E(0) + E(1))
      // We Replace E(i) with 0.5*(E(i) - E(0))

      for(i=0; i<InputData.size(); i++)
	{
	  for(j=0; j<List_Delta.size(); j++)
	    {
	      for(k=1; k<List_Delta[j].ListOfInputs.size(); k++)
		{
		  InputData[i][List_Delta[j].ListOfInputs[k]] = 0.5*(InputData[i][List_Delta[j].ListOfInputs[k]] - 
								     InputData[i][List_Delta[j].ListOfInputs[0]]);
		} /* End For */
	      InputData[i][List_Delta[j].ListOfInputs[0]] = 0.5*(InputData[i][List_Delta[j].ListOfInputs[1]] + 
								 InputData[i][List_Delta[j].ListOfInputs[0]]);
	    } /* End For */
	} /* End For */
    } /* End If */

  if ((unsigned int)Measure>=MeasureData[0].size())
    {
      cerr << "Measure must be included between 0 and " << MeasureData[0].size() << endl;

      // We change the current working directory
      chdir(Current_Dir.c_str());

      exit(1);
    } /* End If */

  Learn_InputMin   = InputMin;
  Learn_InputMax   = InputMax;
  Learn_MeasureMin = MeasureMin;
  Learn_MeasureMax = MeasureMax;

  if (UseNormalize)
    {
      FileName.str("");
      FileName << ResultDir << "_Train_" << Norm_Filename;
      
      Filter_Normalize(InputData, MeasureData, InputMin, InputMax,
		       MeasureMin, MeasureMax, FileName.str(), true);

      // Modification of the mins and maxs of the inputs / outputs
      for(i=0; i<InputMin.size(); i++)
	{
	  InputMin[i] = 0.0;
	  InputMax[i] = 1.0;
	} /* End For */
      
      for(i=0; i<MeasureMin.size(); i++)
	{
	  MeasureMin[i] = 0.0;
	  MeasureMax[i] = 1.0;
	} /* End For */
    } /* End If */

  if (Display)
    {
      Display_DisplayData(InputData, MeasureData);
    } /* End If */

  cout << "Learning - Sampling the data" << endl;
  
  if (MeasureData.size()!=InputData.size())
    {
      cerr << "Learning - Number of input lines is different than the number of output lines: ";
      cerr << InputData.size() << " - " << MeasureData.size() << endl;

      // We change the current working directory
      chdir(Current_Dir.c_str());

      exit(1);
    } /* End If */


  if (Use_SkipLearn)
    {
      cout << "Learning - Skip Sequences" << endl;

      Skip(InputData, MeasureData, List_SkipLearn, 0);
    } /* End If */

  if (List_Filter.size()>0)
    {
      for(i=0; i<List_Filter.size(); i++)
	{
	  cout << "Learning - Sampling with respect to the input " << List_Filter[i].InputNb << endl;
	  cout << "Learning - between " << List_Filter[i].Min << " and " << List_Filter[i].Max << endl;

	  if (List_Filter[i].Add)
	    {
	      Filter_Echant_Input_Add(InputData, MeasureData, InputMin, InputMax,
				      List_Filter[i].Min, List_Filter[i].Max, List_Filter[i].InputNb, Echant_AddMinMax);
	    } /* End If */
	  else
	    {
	      Filter_Echant_Input_Remove(InputData, MeasureData, InputMin, InputMax,
					 List_Filter[i].Min, List_Filter[i].Max, List_Filter[i].InputNb, Echant_AddMinMax);
	    } /* End Else */
	} /* End For */
    } /* End If */

  if ((Echant_Output_NbPts>0)&&(Echant_Output_UseMin))
    {
      cout << "Learning - Sampling with respect to the output : " << Echant_Output_NbPts;
      cout << " points of smallest values" << endl;
      Filter_Echant_Output(InputData, MeasureData, MeasureMin, MeasureMax, Measure, Echant_Output_NbPts, true,
			   Echant_AddMinMax);
    } /* End If */

  if ((Echant_Output_NbPts>0)&&(!Echant_Output_UseMin))
    {
      cout << "Learning - Sampling with respect to the output : " << Echant_Output_NbPts;
      cout << " points of biggest values" << endl;
      Filter_Echant_Output(InputData, MeasureData, MeasureMin, MeasureMax, Measure, Echant_Output_NbPts, false,
			   Echant_AddMinMax);
    } /* End If */

  if (Echant_Rand_NbPoints>0)
    {
      cout << "Learning - Random sampling of " << Echant_Rand_NbPoints << " points" << endl;
      Filter_EchantRand(InputData, MeasureData, Echant_Rand_NbPoints, Echant_AddMinMax);
    } /* End If */

  if (Echant_Stat_NbBar>0)
    {
      cout << "Learning - Static sampling : Histogram = " << Echant_Stat_NbBar << " to find " << Echant_Stat_NbPoints;
      cout << " points" << endl;

      Filter_EchantStat(InputData, MeasureData, Echant_Stat_NbBar, Echant_Stat_NbPoints, Echant_AddMinMax);
    } /* End If */

  if (Echant_Dyn_NbBar>0)
    {
      cout << "Learning - Dynamic sampling : Histogram = " << Echant_Dyn_NbBar << " to find " << Echant_Dyn_NbPoints;
      cout << " points" << endl;

      Filter_EchantDyn(InputData, MeasureData, Echant_Dyn_NbBar, Echant_Dyn_NbPoints, Echant_AddMinMax);
    } /* End If */

  if (Echant_Delta_NbPoints>0)
    {
      cout << "Learning - Delta sampling : " << Echant_Delta_NbPoints << " points 1 " << Echant_Delta_NbTimes << " times" << endl;

      Filter_EchantDelta(InputData, MeasureData, Echant_Delta_NbPoints, Echant_Delta_NbTimes, Echant_AddMinMax, Measure);
    } /* End If */

  if (Echant_Greedy_Fast_NbPointsAtTheEnd>0)
    {
      cout << "Learning - Fast greedy sampling of " << Echant_Greedy_Fast_NbPointsAtTheEnd;
      cout << " points in " << Echant_Greedy_Fast_NbPointsToSelect << " times" << endl;
      Filter_EchantGreedy_Fast(InputData, MeasureData, InputMin, InputMax, 
			       Echant_Greedy_Fast_NbPointsToSelect,
			       Echant_Greedy_Fast_NbPointsAtTheEnd,
			       Echant_AddMinMax);
    } /* End If */

  if (Echant_Greedy_Fast_Add_Measure_NbPointsAtTheEnd>0)
    {
      cout << "Learning - Fast greedy sampling of " << Echant_Greedy_Fast_Add_Measure_NbPointsAtTheEnd;
      cout << " points in " << Echant_Greedy_Fast_Add_Measure_NbPointsToSelect << " times (adding the measures)" << endl;
      Filter_EchantGreedy_Fast_AddMeasure(InputData, MeasureData, InputMin, InputMax, MeasureMin, MeasureMax,
					  Echant_Greedy_Fast_Add_Measure_NbPointsToSelect,
					  Echant_Greedy_Fast_Add_Measure_NbPointsAtTheEnd,
					  Echant_AddMinMax, Measure);
    } /* End If */

  if (Echant_Greedy_Partial_NbPtsToSelect>0)
    {
      cout << "Learning - Partiel fast greedy sampling of " << Echant_Greedy_Partial_NbPtsToSelect;
      cout << " points in " << Echant_Greedy_Partial_NbTimesToSelect << " times (adding the measures)" << endl;
      Filter_EchantPartialGreedy_AddMeasure(InputData, MeasureData, InputMin, InputMax, MeasureMin, MeasureMax,
					    Echant_Greedy_Partial_NbPtsToSelect,
					    Echant_Greedy_Partial_NbTimesToSelect,
					    Echant_AddMinMax, Measure);
    } /* End If */

  if (Echant_PartGreedDelta_NbPoints>0)
    {
      cout << "Learning - Delta partiel greedy sampling of " << Echant_PartGreedDelta_NbPoints;
      cout << " points in " << Echant_PartGreedDelta_NbTimes << " times" << endl;
      Filter_EchantPartialGreedyDelta(InputData, MeasureData, InputMin, InputMax, MeasureMin, MeasureMax,
				      Echant_PartGreedDelta_NbPoints,
				      Echant_PartGreedDelta_NbTimes,
				      Echant_AddMinMax, Measure);
    } /* End If */

  if (Echant_Greedy_NbPoints>0)
    {
      cout << "Learning - Greedy sampling of " << Echant_Greedy_NbPoints << " points" << endl;
      Filter_EchantGreedy(InputData, MeasureData, InputMin, InputMax, Echant_Greedy_NbPoints, Echant_AddMinMax);
    } /* End If */

  if (Echant_Greedy_Add_Measure_NbPoints>0)
    {
      cout << "Learning - Greedy sampling of " << Echant_Greedy_Add_Measure_NbPoints << " points ";
      cout << "(adding the measures)" << endl;
      Filter_EchantGreedy_AddMeasure(InputData, MeasureData,
				     InputMin, InputMax, MeasureMin, MeasureMax, 
				     Echant_Greedy_Add_Measure_NbPoints,
				     Echant_AddMinMax, Measure);
    } /* End If */

  if (Display)
    {
      Display_DisplayData(InputData, MeasureData);
    } /* End If */

  if (UseSampledData)
    {
      FileName.str("");
      FileName << ResultDir << SampledDataFile;

      ExportFile = new ofstream;
      ExportFile->open(FileName.str().c_str());

      for(i=0; i<InputData.size(); i++)
	{
	  for(j=0; j<InputData[i].size(); j++)
	    {
	      (*ExportFile) << InputData[i][j] << " ";
	    } /* End For */
	  (*ExportFile) << MeasureData[i][Measure] << endl;
	} /* End For */

      ExportFile->close();

      delete ExportFile;
    } /* End If */

  if (SampledData_Exit)
    {
      cout << "TrainLolimot2:: Exiting after data sampling" << endl;
      exit(1);
    } /* End If */

  cout << "Learning - Initialisation of the LOLIMOT network" << endl;
  
  Lolimot = new LL_Lolimot;
  
  Lolimot->cleanMesures();
  Lolimot->cleanPartitions();
  Lolimot->cleanDimensions();
  
  Lolimot->setUseTransformLog(UseTransformLog);
  Lolimot->setTransformLogAlpha(UseTransformLogAlpha);
  Lolimot->setTransformLogEps(UseTransformLogEps);
  Lolimot->setTransformLogMin(MeasureMin[Measure]);
  Lolimot->setTransformLogMax(MeasureMax[Measure]);
  Lolimot->setMembershipThreshold(MembershipThreshold);
  
  if (Load_Model)
    {
      unsigned int Load_NbCoeff, Load_NbDimension, Load_NbPartition;
      string       Load_Name;
      
      FileName.str("");
      FileName << ResultDir << Load_Model_Filename;
      
      cout << "Loading the Lolimot network : " << FileName.str() << endl;
      
      ModelFile_Load = new ifstream;
      ModelFile_Load->open(FileName.str().c_str());
      
      if (!ModelFile_Load->is_open())
	{
	  cerr << "Load Model Error: can't open file " << FileName.str() << endl;
	  
	  // We change the current directory
	  chdir(Current_Dir.c_str());
	  
	  exit(1);
	} /* End If */
      
      (*ModelFile_Load) >> Load_NbPartition;
      (*ModelFile_Load) >> Load_NbDimension;
      (*ModelFile_Load) >> Load_NbCoeff;
      
      Lolimot->cleanPartitions();
      Lolimot->cleanDimensions();
      Lolimot->cleanMesures();
      
      for(i=0; i<Load_NbPartition; i++)
	{
	  Lolimot->addPartition();
	} /* End For */
      
      for(i=0; i<Load_NbPartition; i++)
	{
	  Lolimot->getPartition(i)->setNbDimension(Load_NbDimension);
	  
	  // We load the boundaries of the partitions
	  for(j=0; j<Load_NbDimension; j++)
	    {
	      (*ModelFile_Load) >> StdResidu;
	      Lolimot->getPartition(i)->setDimensionMin(j, StdResidu);
	      (*ModelFile_Load) >> StdResidu;
	      Lolimot->getPartition(i)->setDimensionMax(j, StdResidu);
	      (*ModelFile_Load) >> StdResidu;
	      Lolimot->getPartition(i)->setInhibe(j, (bool)StdResidu);
	    } /* End For */
	  
	  // We store the model corresponding to the partition
	  (*ModelFile_Load) >> StdResidu;
	  Lolimot->getPartition(i)->setCoeff0(StdResidu);
	  (*ModelFile_Load) >> StdResidu;
	  Lolimot->getPartition(i)->setFreezeCoeff0((bool)StdResidu);
	  
	  for(k=0; k<Load_NbCoeff; k++)
	    {
	      (*ModelFile_Load) >> StdResidu;
	      Lolimot->getPartition(i)->setCoeff(k, StdResidu);
	      (*ModelFile_Load) >> StdResidu;
	      Lolimot->getPartition(i)->setFreezeCoeff(k, (bool)StdResidu);
	    } /* End For */
	  
	  // We store some complementary informations related to the partitions
	  (*ModelFile_Load) >> StdResidu;
	  Lolimot->getPartition(i)->setResidu(StdResidu);
	  (*ModelFile_Load) >> StdResidu;
	  Lolimot->getPartition(i)->SetDoNotCut((bool)StdResidu);
	  (*ModelFile_Load) >> Load_Name;
	  Lolimot->getPartition(i)->setName(Load_Name);
	} /* End For */
      
	  // We store some informations related to the dimensions
      
      for(i=0; i<Load_NbDimension; i++)
	{
	  Lolimot->addDimension("", 0.0, 0.0, 0.0, 0);
	} /* End For */
      
      for(i=0; i<Load_NbDimension; i++)
	{
	  (*ModelFile_Load) >> StdResidu;
	  Lolimot->getDimension(i)->setMin(StdResidu);
	  (*ModelFile_Load) >> StdResidu;
	  Lolimot->getDimension(i)->setMax(StdResidu);
	  (*ModelFile_Load) >> StdResidu;
	  Lolimot->getDimension(i)->setEcartType(StdResidu);
	  (*ModelFile_Load) >> StdResidu;
	  Lolimot->getDimension(i)->setNbDiscretisations((int)StdResidu);
	  (*ModelFile_Load) >> Load_Name;
	  Lolimot->getDimension(i)->setName(Load_Name);
	  (*ModelFile_Load) >> StdResidu;
	  Lolimot->setDontCutVar(i, (bool)StdResidu);
	} /* End For */
      
      // We store some informations related to the lolimot network
      (*ModelFile_Load) >> StdResidu;
      Lolimot->setNbMaxPartitions((int)StdResidu);
      (*ModelFile_Load) >> StdResidu;
      Lolimot->setResiduGapPercentage(StdResidu);
      (*ModelFile_Load) >> StdResidu;
      Lolimot->setSigma(StdResidu);
      
      (*ModelFile_Load) >> Load_Name;
      if (Load_Name=="none")
	Lolimot->setCustomArguments("");
      else
	Lolimot->setCustomArguments(Load_Name);
      
      (*ModelFile_Load) >> Load_Name;
      if (Load_Name=="none")
	Lolimot->setCustomReturn_Before("");
      else
	Lolimot->setCustomReturn_Before(Load_Name);
      
      (*ModelFile_Load) >> Load_Name;
      if (Load_Name=="none")
	Lolimot->setCustomReturn_After("");
      else
	Lolimot->setCustomReturn_After(Load_Name);
      
      (*ModelFile_Load) >> StdResidu;
      Lolimot->setUseTransformLog((bool)StdResidu);
      (*ModelFile_Load) >> StdResidu;
      Lolimot->setTransformLogAlpha((double)StdResidu);
      (*ModelFile_Load) >> StdResidu;
      Lolimot->setMembershipThreshold((double)StdResidu);
      

      (*ModelFile_Load) >> StdResidu;
      Lolimot->setTransformLogMin(StdResidu);
      (*ModelFile_Load) >> StdResidu;
      Lolimot->setTransformLogMax(StdResidu);
      (*ModelFile_Load) >> StdResidu;
      Lolimot->setTransformLogEps(StdResidu);

      ModelFile_Load->close();
      
      delete ModelFile_Load;
    } /* End If */
    
  Lolimot->updatePartitions();
  
  Lolimot->setNbMaxPartitions(NbMaxPart);
  Lolimot->setResiduGapPercentage(ResiduGapPercentage);
  Lolimot->setUseLasso(UseLasso);

  Lolimot->setNbParam(0);
  
  switch(ResidualType)
    {
    case 1:
      Lolimot->addFctResidu(DefFctResidu_L1);
      break;
    case 2:
      Lolimot->addFctResidu(DefFctResidu_L1_Rel);
      break;
    case 3:
      Lolimot->addFctResidu(DefFctResidu_L2);
      break;
    case 4:
      Lolimot->addFctResidu(DefFctResidu_L2_Rel);
      break;
    case 5:
      Lolimot->addFctResidu(DefFctResidu_NoNeg);
      break;
    default:
      Lolimot->addFctResidu(DefFctResidu);
      break;
    } /* End Switch */
    
  Lolimot->setUseDistributedCutting(UseDistribCutting);
  Lolimot->setUseUniformCutting(UseUniformCutting);
  
  if (CustomParam)      Lolimot->setCustomArguments(CustParamStr);
  if (CustomRet_Before) Lolimot->setCustomReturn_Before(CustRetBefStr);
  if (CustomRet_After)  Lolimot->setCustomReturn_After(CustRetAfterStr);
  
  // Manual cutting
  if (List_PartToCut.size()!=0)
    {
      cout << "Learning : Performing a manual cutting" << endl;
      
      cout << "Manual cutting: verification of the coherency of the cuttings" << endl;
      
      for(i=0; i<List_PartToCut.size(); i++)
	{
	  if (InputMin[List_PartToCut[i].Dimension]>List_PartToCut[i].CutPosition)
	    {
	      cerr << "ERROR: the min value of the dimension " << List_PartToCut[i].Dimension;
	      cerr << " is greater than the position of the cutting." << endl;
	      cerr << "Error on the cutting number " << i << endl;
	      
	      // We change the current directory
	      chdir(Current_Dir.c_str());
	      
	      exit(1);
	    } /* End If */
	} /* End For */
      
      for(i=0; i<List_Model_Input.size(); i++)
	{
	  Lolimot->addDimension(List_Model_Input[i].Name.c_str(), InputMin[i], InputMax[i], Sigma, NbCut);
	  if (List_Model_Input[i].DontCut) Lolimot->setDontCutVar(i, true);
	} /* End For */
      
      Load_Model = true;
      
      for(i=0; i<List_PartToCut.size(); i++)
	{
	  if (Lolimot->getPartitionSet().size()==0)
	    {
	      // We create a second partition for the initialization
	      LL_Partition * newPartition_Low = Lolimot->addPartition();
	      LL_Partition * newPartition_Up  = Lolimot->addPartition();
	      LL_Dimension * curDimension     = NULL;
	      
	      // We initialize this partition with the bounds of the dimensions
	      for(j = 0; j < Lolimot->getDimensionSet().size(); j++)
		{
		  curDimension = Lolimot->getDimension(j);
		  
		  newPartition_Low->setDimensionMin(j, curDimension->getMin());
		  newPartition_Low->setDimensionMax(j, curDimension->getMax());
		  
		  newPartition_Up->setDimensionMin(j, curDimension->getMin());
		  newPartition_Up->setDimensionMax(j, curDimension->getMax());
		} /* End For */
	      
	      newPartition_Low->setDimensionMax(List_PartToCut[i].Dimension,
						List_PartToCut[i].CutPosition);
	      newPartition_Up->setDimensionMin(List_PartToCut[i].Dimension,
					       List_PartToCut[i].CutPosition);
	      newPartition_Low->setName(List_PartToCut[i].PartName_LowerPart);
	      newPartition_Up->setName(List_PartToCut[i].PartName_UpperPart);
	      
	      newPartition_Low->updateEcartType(Lolimot);
	      newPartition_Up->updateEcartType(Lolimot);
	    } /* End If */
	  else
	    {
	      for(j=0; j<Lolimot->getPartitionSet().size(); j++)
		{
		  if (List_PartToCut[i].PartName==Lolimot->getPartition(j)->getName())
		    {
		      // We add a new partition
		      LL_Partition * AuxPart = Lolimot->addPartition();
		      AuxPart->setDimensions(Lolimot->getPartition(j), Lolimot);
		      
		      Lolimot->getPartition(j)->setDimensionMax(List_PartToCut[i].Dimension,
								List_PartToCut[i].CutPosition);
		      AuxPart->setDimensionMin(List_PartToCut[i].Dimension,
					       List_PartToCut[i].CutPosition);
		      
		      Lolimot->getPartition(j)->setName(List_PartToCut[i].PartName_LowerPart);
		      AuxPart->setName(List_PartToCut[i].PartName_UpperPart);
		      
		      Lolimot->getPartition(j)->updateEcartType(Lolimot);
		      AuxPart->updateEcartType(Lolimot);
		      
		      j = Lolimot->getPartitionSet().size(); // We quit the loop
		    } /* End If */
		} /* End For */
	    } /* End Else */
	} /* End For */
      
      // Updating the models
      Lolimot->updatePartitions();
      
      LL_Partition * AuxPart = NULL;
      
      std::cerr << "List of manually created partitions" << std::endl;
      
      for(i=0; i<Lolimot->getPartitionSet().size(); i++)
	{
	  AuxPart = Lolimot->getPartition(i);
	  
	  std::cerr << "n " << i+1 << " " << AuxPart->printDimensions() << std::endl;
	} /* End For */
    } /* End If */
  
  for(i=0; i<List_ModelCoeff.size(); i++)
    {
      if (Lolimot->getPartitionSet().size()==0)
	{
	  cerr << "ERROR: There is no partitions !" << endl;
	  exit(1);
	} /* End If */
      else
	{
	  for(j=0; j<Lolimot->getPartitionSet().size(); j++)
	    {
	      if (List_ModelCoeff[i].PartName==Lolimot->getPartition(j)->getName())
		{
		  // We modify the coefficients of the partition
		  if (List_ModelCoeff[i].Dimension==-1)
		    {
		      Lolimot->getPartition(j)->setFreezeCoeff0(true);
		      Lolimot->getPartition(j)->setCoeff0(List_ModelCoeff[i].Coeff);
		    } /* End If */
		  else
		    {
		      Lolimot->getPartition(j)->setFreezeCoeff(List_ModelCoeff[i].Dimension, true);
		      Lolimot->getPartition(j)->setCoeff(List_ModelCoeff[i].Dimension, List_ModelCoeff[i].Coeff);
		    } /* End Else */
		} /* End If */
	    } /* End For */
	} /* End Else */
    } /* End For */
  
  if (PrintPart)
    {
      cout << "List of partitions" << endl;
      for(i=0; i<Lolimot->getPartitionSet().size(); i++)
	{
	  Lolimot->getPartition(i)->printDimensions();
	} /* End For */
    } /* End If */
  
  cout << "Learning the LOLIMOT network" << endl;
  
  if (BS_Light)
    {
      cout << "Learning - Bootstrap Light method selected" << endl;
      
      Learning_Bootstrap_Light(Lolimot, InputData, MeasureData, WeightData, List_Model_Input, MeasureMin, MeasureMax,
			       InputMin, InputMax, Sigma, NbCut,
			       BS_Proportion, NbPartToExplore, Measure, ResultDir, Res_Filename,
			       Load_Model);
      
    } /* End If */
  
  if (Learn_Recur)
    {
      cout << "Learning - Recursive method selected" << endl;
      
      Learning_Recursive(Lolimot, InputData, MeasureData, WeightData, List_Model_Input,
			 MeasureMin, MeasureMax, InputMin, InputMax,
			 Sigma, NbCut, Rec_NbPart_Ph1, Rec_NbPart_Ph2, Measure, 
			 Load_Model);
    } /* End If */
  
  if (BS_NoMeth)
    {
      cout << "Learning - No method selected" << endl;
      
      Learning_NoMeth(Lolimot, InputData, MeasureData, WeightData, List_Model_Input, InputMin, InputMax,
		      MeasureMin, MeasureMax, Sigma, NbCut, Measure, Load_Model);
    } /* End If */

  if (BS_Classic)
    {
      cout << "Learning - Classical method selected" << endl;
      
      Learning_Classic(Lolimot, InputData, MeasureData, WeightData, List_Model_Input, InputMin, InputMax,
		       MeasureMin, MeasureMax, Sigma, NbCut, Measure, Load_Model);
    } /* End If */
  
  if (UseUpdate)
    {
      cout << "Learning - Update method selected" << endl;
      
      Learning_Update(Lolimot, InputData, MeasureData, WeightData, List_Model_Input, InputMin, InputMax,
		      MeasureMin, MeasureMax, Sigma, NbCut, Measure, Nb_Iter, Load_Model);
    } /* End If */

  if (BS_ComputeR2)
    {
      cout << "Learning - Compute R2 method selected" << endl;
      
      Learning_Compute_R2(Lolimot, InputData, MeasureData, WeightData, List_Model_Input, MeasureMin, MeasureMax,
			  InputMin, InputMax, Sigma, NbCut, Measure, ResultDir, R2_Filename);
    } /* End If */
  
  // We search for the best sigma value for this validation file
  if (UseOptimizeSigma)
    {
      cout << "Learning - Optimisation of the Sigma parameter" << endl;
      
      Optimize_Sigma(Lolimot, InputData, MeasureData, WeightData, List_Model_Input, Measure, SigmaMin, SigmaMax, SigmaStep);
    } /* End If */
  
  if (Use_LOO)
    {
      cout << "Learning - Leave One Out" << endl;
      
      FileName.str("");
      FileName << ResultDir << LOO_Filename;
      
      Learning_LeavOneOut(Lolimot, InputData, MeasureData, WeightData, List_Model_Input,
			  InputMin, InputMax, MeasureMin, MeasureMax, Measure, FileName.str());
    } /* End If */
  
  if (UsePartFilename)
    {
      cout << "Learning - Exportation of all the partitions" << endl;
      
      Lolimot->exportAllPartitions(ResultDir, Part_Filename);
    } /* End If */
  
  if (PrintDistribOfCut)
    {
      string AuxName = ResultDir + DistribOfCutName;
      
      cout << "Learning - Exportation of the distribution of cuttings in the file " << AuxName << endl;
      
      Analysis_List_DistribOfCut(Lolimot, AuxName);
    } /* End If */
  
  if (CrossValid)
    {
      Learning_CrossValidation(Lolimot, InputData, MeasureData, WeightData, List_Model_Input,
			       MeasureMin, MeasureMax, InputMin, InputMax, CrossValid_NbFiles, Measure, ResultDir,
			       CrossValid_Filename, CrossValid_MeanResidual, CrossValid_StdResidual);
    } /* End If */
  
  if (Save_Model)
    {
      FileName.str("");
      FileName << ResultDir << Save_Model_Filename;
      
      cout << "Saving the LOLIMOT network : " << FileName.str() << endl;
      
      ModelFile_Save = new ofstream;
      ModelFile_Save->open(FileName.str().c_str());

      if (!ModelFile_Save->is_open())
	{
	  cerr << "Save Model Error: can't open file " << FileName.str() << endl;

	  // We change the current directory
	  chdir(Current_Dir.c_str());

	  exit(1);
	} /* End If */

      (*ModelFile_Save) << Lolimot->getPartitionSet().size() << endl;
      (*ModelFile_Save) << Lolimot->getDimensionSet().size() << endl;
      (*ModelFile_Save) << Lolimot->getPartition(0)->getCoeffSet().size() << endl;
      
      for(i=0; i<Lolimot->getPartitionSet().size(); i++)
	{
	  // We store the bounds of the partitions
	  for(j=0; j<Lolimot->getDimensionSet().size(); j++)
	    {
	      (*ModelFile_Save) << Lolimot->getPartition(i)->getDimensionMin(j) << endl;
	      (*ModelFile_Save) << Lolimot->getPartition(i)->getDimensionMax(j) << endl;
	      (*ModelFile_Save) << Lolimot->getPartition(i)->getInhibe(j) << endl;
	    } /* End For */

	  // We store the model corresponding to the partition
	  (*ModelFile_Save) << Lolimot->getPartition(i)->getCoeff0() << endl;
	  (*ModelFile_Save) << Lolimot->getPartition(i)->getFreezeCoeff0() << endl;

	  for(k=0; k<Lolimot->getPartition(i)->getCoeffSet().size(); k++)
	    {
	      (*ModelFile_Save) << Lolimot->getPartition(i)->getCoeff(k) << endl;
	      (*ModelFile_Save) << Lolimot->getPartition(i)->getFreezeCoeff(k) << endl;
	    } /* End For */
	  
	  // We store some complementary informations related to the partitions
	  (*ModelFile_Save) << Lolimot->getPartition(i)->getResidu() << endl;
	  (*ModelFile_Save) << Lolimot->getPartition(i)->DoNotCut() << endl;
	  (*ModelFile_Save) << Lolimot->getPartition(i)->getName() << endl;
	} /* End For */

      // We store some informations related to the dimensions
      for(i=0; i<Lolimot->getDimensionSet().size(); i++)
	{
	  (*ModelFile_Save) << Lolimot->getDimension(i)->getMin() << endl;
	  (*ModelFile_Save) << Lolimot->getDimension(i)->getMax() << endl;
	  (*ModelFile_Save) << Lolimot->getDimension(i)->getEcartType() << endl;
	  (*ModelFile_Save) << Lolimot->getDimension(i)->getNbDiscretisations() << endl;
	  (*ModelFile_Save) << Lolimot->getDimension(i)->getName() << endl;
	  (*ModelFile_Save) << Lolimot->getDontCutVar(i) << endl;

	} /* End For */

      // We store some informations related to the LOLIMOT network
      (*ModelFile_Save) << Lolimot->getNbMaxPartitions() << endl;
      (*ModelFile_Save) << Lolimot->getResiduGapPercentage() << endl;
      (*ModelFile_Save) << Lolimot->getSigma() << endl;

      // We store the customisation
      if (Lolimot->getCustomArguments()!="")
	(*ModelFile_Save) << Lolimot->getCustomArguments() << endl;
      else
	(*ModelFile_Save) << "none" << endl;

      if (Lolimot->getCustomReturn_Before()!="")
	(*ModelFile_Save) << Lolimot->getCustomReturn_Before() << endl;
      else
	(*ModelFile_Save) << "none" << endl;

      if (Lolimot->getCustomReturn_After()!="")
	(*ModelFile_Save) << Lolimot->getCustomReturn_After() << endl;
      else
	(*ModelFile_Save) << "none" << endl;

      (*ModelFile_Save) << Lolimot->getUseTransformLog() << endl;
      (*ModelFile_Save) << Lolimot->getTransformLogAlpha() << endl;
      (*ModelFile_Save) << Lolimot->getMembershipThreshold() << endl;
      (*ModelFile_Save) << Lolimot->getTransformLogMin() << endl;
      (*ModelFile_Save) << Lolimot->getTransformLogMax() << endl;
      (*ModelFile_Save) << Lolimot->getTransformLogEps() << endl;

      ModelFile_Save->close();

      delete ModelFile_Save;
    } /* End If */

  // Validation of the lolimot net

  Validate(Lolimot,
	   InputData, MeasureData, MeasureMin, MeasureMax, InputMin, InputMax, Measure, 
	   ResultDir, ModelName, TempoInput, List_Model_Output, List_Model_Input,
	   List_Retro_Input, Learn_Start, Learn_End, Display,
	   CrossValid, CrossValid_NbFiles, CrossValid_MeanResidual, CrossValid_StdResidual,
	   "Learning", UseTransformLogEps, UseTransformLogAlpha, false, 0,
	   Compute_C1, Compute_C2);
    
  if (ExportMatlab)
    {
      FileName.str("");
      FileName << ModelName << List_Model_Output[Measure].Name;
      cout << "Export LOLIMOT model to matlab : " << FileName.str() << endl;
      exportFunctionInMatlab(Lolimot, ResultDir, FileName.str(), List_Model_Input);
    } /* End If */
  
  if (ExportCpp)
    {
      FileName.str("");
      FileName << ModelName << List_Model_Output[Measure].Name;
      cout << "Export LOLIMOT model to c++ : " << FileName.str() << endl;
      exportFunctionInCpp(Lolimot, ResultDir, FileName.str(), List_Model_Input);
      exportHeaderInCpp(Lolimot, ResultDir, FileName.str(), List_Model_Input);
    } /* End If */
  
  if (ExportC)
    {
      FileName.str("");
      FileName << ModelName << List_Model_Output[Measure].Name;
      cout << "Export LOLIMOT model to c : " << FileName.str() << endl;
      exportFunctionInC(Lolimot, ResultDir, FileName.str(), List_Model_Input);
      exportHeaderInC(Lolimot, ResultDir, FileName.str(), List_Model_Input);
    } /* End If */
  
  if (D_ExportMatlab2)
    {
      FileName.str("");
      FileName << "D2_" << ModelName << List_Model_Output[Measure].Name;
      cout << "Export type 2 partial derivatives LOLIMOT model to matlab : " << FileName.str() << endl;
      exportDerivativeFunctionInMatlab2(Lolimot, ResultDir, FileName.str(), List_Model_Input);
    } /* End If */
  
  if (D_ExportCpp2)
    {
      FileName.str("");
      FileName << "D2_" << ModelName << List_Model_Output[Measure].Name;
      cout << "Export type 2 partial derivatives LOLIMOT model to c++ : " << FileName.str() << endl;
      exportDerivativeFunctionInCpp2(Lolimot, ResultDir, FileName.str(), List_Model_Input);
      exportHeaderInCpp2(Lolimot, ResultDir, FileName.str(), List_Model_Input);
    } /* End If */
  
  if (D_ExportC2)
    {
      FileName.str("");
      FileName << "D2_" << ModelName << List_Model_Output[Measure].Name;
      cout << "Export type 2 partial derivatives LOLIMOT model to c : " << FileName.str() << endl;
      exportDerivativeFunctionInC2(Lolimot, ResultDir, FileName.str(), List_Model_Input);
      exportHeaderInC2(Lolimot, ResultDir, FileName.str(), List_Model_Input);
    } /* End If */
  
  if (DD_ExportMatlab2)
    {
      FileName.str("");
      FileName << "DD2_" << ModelName << List_Model_Output[Measure].Name;
      cout << "Export type 2 second partial derivatives LOLIMOT model to matlab : " << FileName.str() << endl;
      exportSecondDerivativeFunctionInMatlab2(Lolimot, ResultDir, FileName.str(), List_Model_Input);
    } /* End If */
  
  if (DD_ExportCpp2)
    {
      FileName.str("");
      FileName << "DD2_" << ModelName << List_Model_Output[Measure].Name;
      cout << "Export type 2 second partial derivatives LOLIMOT model to c++ : " << FileName.str() << endl;
      exportSecondDerivativeFunctionInCpp2(Lolimot, ResultDir, FileName.str(), List_Model_Input);
      exportHeaderInCpp2(Lolimot, ResultDir, FileName.str(), List_Model_Input);
    } /* End If */
  
  if (DD_ExportC2)
    {
      FileName.str("");
      FileName << "DD2_" << ModelName << List_Model_Output[Measure].Name;
      cout << "Export type 2 second partial derivatives LOLIMOT model to c : " << FileName.str() << endl;
      exportSecondDerivativeFunctionInC2(Lolimot, ResultDir, FileName.str(), List_Model_Input);
      exportHeaderInC2(Lolimot, ResultDir, FileName.str(), List_Model_Input);
    } /* End If */
  
  for(i=0; i<Lolimot->getDimensionSet().size(); i++)
    {
      if (D_ExportMatlab)
	{
	  FileName.str("");
	  FileName << "D_" << i << "_" << ModelName << List_Model_Output[Measure].Name;
	  cout << "Export partial derivatives LOLIMOT model to matlab : " << FileName.str() << endl;
	  exportDerivativeFunctionInMatlab(Lolimot, ResultDir, FileName.str(), i+1, List_Model_Input);
	} /* End If */
      
      if (D_ExportCpp)
	{
	  FileName.str("");
	  FileName << "D_" << i << "_" << ModelName << List_Model_Output[Measure].Name;
	  cout << "Export partial derivatives LOLIMOT model to c++ : " << FileName.str() << endl;
	  exportDerivativeFunctionInCpp(Lolimot, ResultDir, FileName.str(), i, List_Model_Input);
	  exportHeaderInCpp(Lolimot, ResultDir, FileName.str(), List_Model_Input);
	} /* End If */
      
      if (D_ExportC)
	{
	  FileName.str("");
	  FileName << "D_" << i << "_" << ModelName << List_Model_Output[Measure].Name;
	  cout << "Export partial derivatives LOLIMOT model to c : " << FileName.str() << endl;
	  exportDerivativeFunctionInC(Lolimot, ResultDir, FileName.str(), i, List_Model_Input);
	  exportHeaderInC(Lolimot, ResultDir, FileName.str(), List_Model_Input);
	} /* End If */
    } /* End For */
  
  
  for(i=0; i<Lolimot->getDimensionSet().size(); i++)
    {
      for(j=0; j<Lolimot->getDimensionSet().size(); j++)
	{
	  if (DD_ExportMatlab)
	    {
	      FileName.str("");
	      FileName << "DD_" << i << "_" << j << "_" << ModelName << List_Model_Output[Measure].Name;
	      cout << "Export second partial derivatives LOLIMOT model to matlab : " << FileName.str() << endl;
	      exportSecondDerivativeFunctionInMatlab(Lolimot, ResultDir, FileName.str(), i+1, j+1, List_Model_Input);
	    } /* End If */
	  
	  if (DD_ExportCpp)
	    {
	      FileName.str("");
	      FileName << "DD_" << i << "_" << j << "_" << ModelName << List_Model_Output[Measure].Name;
	      cout << "Export second partial derivatives LOLIMOT model to c++ : " << FileName.str() << endl;
	      exportSecondDerivativeFunctionInCpp(Lolimot, ResultDir, FileName.str(), i, j, List_Model_Input);
	      exportHeaderInCpp(Lolimot, ResultDir, FileName.str(), List_Model_Input);
	    } /* End If */
	  
	  if (DD_ExportC)
	    {
	      FileName.str("");
	      FileName << "DD_" << i << "_" << j << "_" << ModelName << List_Model_Output[Measure].Name;
	      cout << "Export second partial derivatives LOLIMOT model to c : " << FileName.str() << endl;
	      exportSecondDerivativeFunctionInC(Lolimot, ResultDir, FileName.str(), i, j, List_Model_Input);
	      exportHeaderInC(Lolimot, ResultDir, FileName.str(), List_Model_Input);
	    } /* End If */
	} /* End For */
    } /* End For */

  //
  // Cleaning the lists containing some data
  //
	  
  MeasureData.clear();
  InputData.clear();

  if (Validation_File)
    {
      //
      // Validation of the LOLIMOT network
      //
      // Opening the data files
	  
      for(n=0; n<List_VALID_FileName.size(); n++)
	{
	  MeasureData.resize(0);
	  InputData.resize(0);

	  // We allocate as many file descriptors as there are file name on a VALID line
	  DataFile = new ifstream[List_VALID_FileName[n].size()];
	      
	  // We allocate a buffer per file descriptor
	  typedef char * PtrChar;
	      
	  Buffer = new PtrChar[List_VALID_FileName[n].size()];
	      
	  // We open all the files and we initialize the buffers
	  for(m=0; m<List_VALID_FileName[n].size(); m++)
	    {
	      cout << "Opening file " << List_VALID_FileName[n][m] << endl;
		  
	      DataFile[m].open(List_VALID_FileName[n][m].c_str(), ios_base::in | ios_base::binary);
		  
	      if (!DataFile[m].is_open())
		{
		  cerr << "Validation - can't open file " << List_VALID_FileName[n][m];
		  cerr << "." << endl;

		  // We change te current working directory
		  chdir(Current_Dir.c_str());

		  exit(1);
		} /* End If */
		  
	      Buffer[m] = new char[BUFFERSIZE];
	    } /* End For */
	      
	  Result = true;
	  WList.resize(List_VALID_FileName[n].size());
	      
	  // We skip the lines which contain the labels of the columns
	  for(m=0; m<List_VALID_FileName[n].size(); m++)
	    {
	      for(i=0; i<List_VALID_SkipLines[n][m]; i++)
		{
		  DataFile[m].getline(Buffer[m], BUFFERSIZE, '\n');
		} /* End For */
	    } /* End For */

	  LineCount = -1;

	  while(Result)
	    {
	      NbColAux = 0;
	      // We store a line of each file in each buffer
	      for(m=0; m<List_VALID_FileName[n].size(); m++)
		{
		  Result = Result & (bool)DataFile[m].getline(Buffer[m], BUFFERSIZE, '\n');
		  vector<string> tmpLine = tokenize(Buffer[m], is_space());
		  if (tmpLine.size()!=0) WList[m] = tmpLine; // Skip empty lines
		  NbColAux += WList[m].size();
		} /* End For */
		
	      if (WList[m].size()==0)
		{
		  cout << "TrainLolimot2: Empty line " << LineCount + 1 << " ... skipping" << endl;
		  continue;
		} /* End If */

	      if ((Result)&&(NbCol!=NbColAux))
		{
		  cerr << "Error while reading files " << List_VALID_FileName[n][0] << " ... " << endl;
		  cerr << "Error around line " << LineCount << endl;
		  cerr << "Number of column : " << NbCol << endl;
		  cerr << "Number of column read : " << NbColAux << endl;

		  // We change the current working directory
		  chdir(Current_Dir.c_str());

		  exit(-1);
		} /* End If */

	      LineCount++;

	      // We get everything
	      DList.resize(0);

	      for(j=0; j<WList.size(); j++)
		{
		  for(k=0; k<WList[j].size(); k++)
		    {
		      DList.push_back(atof(WList[j][k].c_str()));
		    } /* End For */
		} /* End For */

	      InputData.push_back(DList);
	    } /* End While */

	  for(m=0; m<List_VALID_FileName[n].size(); m++)
	    {
	      DataFile[m].close();
	      if (Buffer[m]) delete [] Buffer[m];
	    } /* End For */
	      
	  if (DataFile) delete [] DataFile;
	  if (Buffer)   delete [] Buffer;

	  cout << "Validation - Reading data - Part 2 - time steps generation" << endl;
  
	  InputMin.resize(List_Model_Input.size());
	  InputMax.resize(List_Model_Input.size());
	  
	  for(i=0; i<List_Model_Input.size(); i++)
	    {
	      InputMin[i] = numeric_limits<float>::max();
	      InputMax[i] = - numeric_limits<float>::max();
	    } /* End For */
	  
	  MeasureMin.resize(List_Model_Output.size());
	  MeasureMax.resize(List_Model_Output.size());
	  
	  for(i=0; i<List_Model_Output.size(); i++)
	    {
	      MeasureMin[i] = numeric_limits<float>::max();
	      MeasureMax[i] = - numeric_limits<float>::max();
	    } /* End For */
	  
	  // Sizing the input temporizations
	  
	  TempoInput.resize(List_Model_Input.size());
	  
	  for(j=0; j<List_Model_Input.size(); j++)
	    {
	      TempoInput[j].resize(List_Model_Input[j].Tempo + ShiftInput, 0.0);
	    } /* End For */
	  
	  // Sizing the output temporizations
	  
	  TempoOutput.resize(List_Model_Output.size());
	  
	  for(j=0; j<List_Model_Output.size(); j++)
	    {
	      TempoOutput[j].resize(List_Model_Output[j].Tempo + ShiftOutput, 0.0);
	    } /* End For */

	  AuxData.resize(InputData.size(), vector<float>(InputData[0].size(), 0.0));
	  
	  AuxData = InputData;
	  
	  InputData.clear();
	  MeasureData.clear();
	  InputData.resize(0);
	  MeasureData.resize(0);

	  for(i=0; i<AuxData.size(); i++)
	    {
	      if (i>=NbMaxTempo)
		{
		  // We get all the inputs
		  DList.resize(0);
		  
		  for(j=0; j<List_Model_Input.size(); j++)
		    {
		      DList.push_back(Temporize(TempoInput[j], AuxData[i][List_Model_Input[j].Pos]));
		    } /* End For */
		  
		  InputData.push_back(DList);
		  
		  for(j=0; j<DList.size(); j++)
		    {
		      if (InputMin[j]>DList[j]) InputMin[j] = DList[j];
		      if (InputMax[j]<DList[j]) InputMax[j] = DList[j];
		    } /* End For */
		  
		  // We get all the outputs
		  DList.resize(0);
		  
		  for(j=0; j<List_Model_Output.size(); j++)
		    {
		      DList.push_back(Temporize(TempoOutput[j], AuxData[i][List_Model_Output[j].Pos]));
		    } /* End For */
		  
		  MeasureData.push_back(DList);
		  
		  for(j=0; j<DList.size(); j++)
		    {
		      if (MeasureMin[j]>DList[j]) MeasureMin[j] = DList[j];
		      if (MeasureMax[j]<DList[j]) MeasureMax[j] = DList[j];
		    } /* End For */
		} /* End If */
	      else if (i<NbMaxTempo)
		{
		  // We store the inputs in the temporization line
		  for(j=0; j<List_Model_Input.size(); j++)
		    {
		      Temporize(TempoInput[j], AuxData[i][List_Model_Input[j].Pos]);
		    } /* End For */
		  
		  // We store the outputs in the temporization line
		  for(j=0; j<List_Model_Output.size(); j++)
		    {
		      Temporize(TempoOutput[j], AuxData[i][List_Model_Output[j].Pos]);
		    } /* End For */
		} /* End Else */
	    } /* End For */

	  AuxData.clear();

	  if (UseDelta)
	    {
	      // Processing inputs
	      Delta_Entry AuxDelta;
	      bool        InputInListDelta = false;
	      int         AuxValue;
	      
	      List_Delta.resize(0);
	      
	      for(i=0; i<List_Model_Input.size()-1; i++)
		{
		  // Is the input i already in the List_Delta list ?
		  
		  InputInListDelta = false;
		  
		  for(j=0; j<List_Delta.size(); j++)
		    {
		      for(k=0; k<List_Delta[j].ListOfInputs.size(); k++)
			{
			  InputInListDelta = InputInListDelta && (i==List_Delta[j].ListOfInputs[k]);
			} /* End For */
		    } /* End For */
		  
		  if (InputInListDelta) continue;
		  
		  for(j=i+1; j<List_Model_Input.size(); j++)
		    {
		      if (List_Model_Input[i].Name==List_Model_Input[j].Name)
			{
			  AuxDelta.ListOfInputs.push_back(j);
			  AuxDelta.ListOfTempo.push_back(List_Model_Input[j].Tempo);
			} /* End If */
		    } /* End For */
		} /* End For */
	      
	      // We order each List_Delta entry with respect to the Tempo
	      
	      for(i=0; i<List_Delta.size(); i++)
		{
		  for(j=0; j<List_Delta[i].ListOfTempo.size()-1; j++)
		    {
		      for(k=j; k<List_Delta[i].ListOfTempo.size(); k++)
			{
			  if (List_Delta[i].ListOfTempo[j]>List_Delta[i].ListOfTempo[k])
			    {
			      AuxValue = List_Delta[i].ListOfTempo[j];
			      List_Delta[i].ListOfTempo[j] = List_Delta[i].ListOfTempo[k];
			      List_Delta[i].ListOfTempo[k] = AuxValue;
			      
			      AuxValue = List_Delta[i].ListOfInputs[j];
			      List_Delta[i].ListOfInputs[j] = List_Delta[i].ListOfInputs[k];
			      List_Delta[i].ListOfInputs[k] = AuxValue;
			    } /* End If */
			} /* End For */
		    } /* End For */
		} /* End For */
	    } /* End If */
	  
	  if (UseNormalize)
	    {
	      FileName.str("");
	      FileName << ResultDir << "_Valid_" << n << "_" << Norm_Filename;
	      
	      InputMin   = Learn_InputMin;
	      InputMax   = Learn_InputMax;
	      MeasureMin = Learn_MeasureMin;
	      MeasureMax = Learn_MeasureMax;

	      Filter_Normalize(InputData, MeasureData, InputMin, InputMax, 
			       MeasureMin, MeasureMax, FileName.str(), false);
	    } /* End If */

	  if (Display)
	    {
	      Display_DisplayData(InputData, MeasureData);
	    } /* End If */
	  
	  cout << "Validation - Sampling data - File " << n << endl;
	  
	  if (MeasureData.size()!=InputData.size())
	    {
	      cerr << "Validation - Number of inputs line is different than the number of output lines : " << InputData.size();
	      cerr << " - " << MeasureData.size() << endl;

	      // We change the working directory
	      chdir(Current_Dir.c_str());

	      exit(1);
	    } /* End If */
	      
	  if (Use_SkipValid)
	    {
	      cout << "Validation - Skip Sequences" << endl;
	      
	      Skip(InputData, MeasureData, List_SkipValid, n);
	    } /* End If */

	  if (Filter_Step>0)
	    {
	      cout << "Validation - Filtering the data" << endl;
	      
	      Filter_Filter(InputData, MeasureData, Filter_Step);
	    } /* End If */
	  
	  if (MAverage>0)
	    {
	      cout << "Validation - Moving average applied on the data" << endl;
	      
	      Filter_MovingAverage(InputData, MeasureData, MAverage);
	    } /* End If */
	  
	  if (Select_Step>0)
	    {
	      cout << "Validation - Selection of a data out of " << Select_Step << endl;
	      
	      Filter_Select(InputData, MeasureData, Select_Step);
	    } /* End If */
	  
	  cout << "Validation - Test of the LOLIMOT network - File : " << List_VALID_FileName[n][0] << endl;
	  
	  Validate(Lolimot,
		   InputData, MeasureData, MeasureMin, MeasureMax, Learn_InputMin, Learn_InputMax, Measure, 
		   ResultDir, ModelName, TempoInput, List_Model_Output, List_Model_Input,
		   List_Retro_Input, Valid_Start, Learn_End, Display,
		   CrossValid, CrossValid_NbFiles, CrossValid_MeanResidual, CrossValid_StdResidual,
		   "Validation", UseTransformLogEps, UseTransformLogAlpha, true, n,
		   Compute_C1, Compute_C2);
	} /* End For n */
    } /* End If */
  
  //
  // We clean some storage area
  //
    
  InputData.clear();
  MeasureData.clear();
      
  Lolimot->cleanDimensions();
  Lolimot->cleanMesures();
  Lolimot->cleanPartitions();
  Lolimot->cleanWeight();

  if (Lolimot) delete Lolimot;

  // We change the working directory
  chdir(Current_Dir.c_str());

  return 0;
}
