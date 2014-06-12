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

#include <fstream>
#include <string>
#include <vector>
#include <ios>
#include <iostream>
#include <cstdlib>
#include <functional>
#include <limits>
#include <sstream>

#include <unistd.h>
#include <string.h>

using namespace std;

#include <LL_Dimension.h>
#include <LL_Mesure.h>
#include <LL_Lolimot.h>

#include <Learning.h>
#include <DefaultParam.h>
#include <DisplayData.h>
#include <Analysis.h>
#include <Post.h>
#include <ExportC.h>
#include <ExportCpp.h>
#include <ExportMatlab.h>
#include <TrainLolimotStruct.h>
#include <Validate.h>

#include <stack-c.h>
#include <api_common.h>
#include <api_double.h>
#include <api_list.h>
#include <api_string.h>

#define NBMAXPARTITIONS     100
#define RESIDUGAPPERCENTAGE 0.0
#define LL_SIGMA            0.33
#define LL_NBCUT            2

#define TRAIN_DATA_IN 1
#define IOLIST_INT    2
#define CUTLIST_INT   3
#define COEFFLIST_INT 4
#define PARAM_IN      5
#define LAST_PARAM    PARAM_IN

#define DEBUG  1

#if !defined(WIN32)
//typedef float (*lolimot_t)(float *);
typedef float (*lolimot_t)(vector<float> &);
#endif

#define CHECK_ERROR_API_SCILAB if(_SciErr.iErr)			    \
                             {					    \
			       printError(&_SciErr, 0);		    \
			       return 0;			    \
			     }					    

// Structure is_comma
// Defines the delimitor as the comma during the analysis of the data files

struct is_comma : public std::unary_function<char,bool>
{
  bool operator() (char c) const
  { return (c==',');}
};

extern "C" int sci_learn_lolimot(char * fname)
{
  vector<vector<float> >  InputData, AuxData;
  vector<vector<float> >  MeasureData;
  vector<float>           WeightData;
  vector<vector<float> >  TempoInput, TempoOutput;
  vector<unsigned int>    List_Retro_Input, Offset_For_Files;
  vector<PartToCut>       List_PartToCut;
  vector<Model_Coeff>     List_ModelCoeff;
  string                  Line, ModelName, ResName, AuxStr, Res_Filename, R2_Filename, Part_Filename;
  string                  Load_Model_Filename, Save_Model_Filename, LOO_Filename, Norm_Filename;
  string                  CustParamStr, CustRetBefStr, CustRetAfterStr, CrossValid_Filename;
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
  stringstream            FileName;
  char                 ** Buffer    = NULL;
  char                  * CfgBuffer = NULL;
  LL_Lolimot            * Lolimot   = NULL;
  float                   StdResidu = 0.0;
  int                     NbMaxPart = -1, NbCut = -1, Measure = -1;
  float                   Sigma = -1.0, ResiduGapPercentage  = -1.0;
  bool                    Result, Display = false;
  bool                    UseDistribCutting = false, UseUniformCutting = true;
  bool                    UseOptimizeSigma = false, UsePartFilename = false, Learn_Recur = false, UseUpdate = false;
  bool                    WeightInputSet = false, PrintPart = false, UseTransformLog = false;
  bool                    CustomParam = false, CustomRet_Before = false, CustomRet_After = false, CrossValid = false;
  int                     Nb_Iter = -1;
  int                     NbPartToExplore = -1;
  float                   SigmaMin = 0.0, SigmaStep = 0.1, SigmaMax = 1.0;
  float                   CrossValid_StdResidual, CrossValid_MeanResidual, UseTransformLogAlpha = 1.0;
  float                   UseTransformLogEps = 0.0, MembershipThreshold = 0.0;
  int                     Rec_NbPart_Ph1 = -1, Rec_NbPart_Ph2 = -1;
  int                     Learn_Start = -1, Learn_End = -1;
  unsigned int            NbMaxTempo = numeric_limits<unsigned int>::min(), RandSeed = 12345, ShiftInput = 0;
  unsigned int            ShiftOutput = 0;
  int                     ResidualType = -1, LineCount = 0, NbCol = 0, NbColAux = 0;
  unsigned int            i, j, k, n, m;
  ScilabStream scicout(std::cout);
  ScilabStream scicerr(std::cerr);
  SciErr _SciErr;
  StrCtx _StrCtx;

  MeasureData.resize(0);
  WeightData.resize(0);

  List_PartToCut.resize(0);
  List_ModelCoeff.resize(0);

  // Opening the configuration file

  List_Model_Output.resize(0);
  List_Model_Input.resize(0);
  
  List_Retro_Input.resize(0);

  // Get the train data set
  int p_train_nb_rows, p_train_nb_cols;
  double * p_train_matrix = NULL;
  int * p_train_in_address = NULL;

  getVarAddressFromPosition(TRAIN_DATA_IN,&p_train_in_address);
  res = getMatrixOfDouble(p_train_in_address, &p_train_nb_rows, &p_train_nb_cols, &p_train_matrix);

  if ((p_learn_nb_rows!=0)&&(p_learn_nb_cols!=0))
    {
      InputData.resize(p_train_nb_rows,vector<double>(p_train_nb_cols, 0.0));
      
      for(i=0;i<p_train_nb_rows; i++)
	{
	  for(j=0;j<p_train_nb_cols; j++)
	    {
	      InputData[i][j] = p_train_matrix[i + j*p_train_nb_rows];
	    }
	}
    }

  // Get the IOLIst
  int * p_iolist_address;
  int * p_iolist_item_address;
  int p_iolist_item_nb_cols, p_iolist_item_nb_rows, p_iolist_nb_items;

  getVarAddressFromPosition(IOLIST_IN,&p_iolist_address);

  res = getListItemNumber(p_iolist_address, &p_iolist_nb_items);

  List_Model_Input.resize(0);

  for(i=1;i<=p_iolist_nb_items;i++)
    {
      int pirow, picol, *pilen, *pilistaddress, size_pistring = 0;
      char ** pstStrings;

      res = getListItemAddress(p_iolist_address, i, &p_iolist_item_address);
      res = getVarDimension(p_iolist_item_address, &p_iolist_item_nb_rows, &p_iolist_item_nb_cols);

      _SciErr = getMatrixOfStringInList(pvApiCtx, sci_x, 2, &pirow, &picol, NULL, NULL);  CHECK_ERROR_API_SCILAB;
      pilen = (int *)MALLOC(pirow*picol*sizeof(int));
      _SciErr= getMatrixOfStringInList(pvApiCtx,sci_x,1,&pirow,&picol,pilen,NULL); CHECK_ERROR_API_SCILAB;
      pstStrings = (char **)MALLOC(pirow*picol*sizeof(char *));
      for(i=0;i<pirow*picol;i++)
	{
	  pstStrings[i] = (char *)MALLOC((pilen[i]+1)*sizeof(char));
	}
      _SciErr = getMatrixOfStringInList(pvApiCtx,sci_x,1,&pirow,&picol,pilen,pstStrings); CHECK_ERROR_API_SCILAB;
      One_Model_Input.Name.str(pstStrings[0]);

      _SciErr = getMatrixOfDoubleInList(pvApiCtx, sci_x, 3, &pirow, &picol, &pdblDataID); CHECK_ERROR_API_SCILAB;
      One_Model_Input.Pos = int(pdblDataID[0]);

      _SciErr = getMatrixOfDoubleInList(pvApiCtx, sci_x, 4, &pirow, &picol, &pdblDataID); CHECK_ERROR_API_SCILAB;
      One_Model_Input.Tempo = int(pdblDataID[0]);

      _SciErr = getMatrixOfDoubleInList(pvApiCtx, sci_x, 5, &pirow, &picol, &pdblDataID); CHECK_ERROR_API_SCILAB;
      One_Model_Input.GroupNb = int(pdblDataID[0]);

      _SciErr = getMatrixOfDoubleInList(pvApiCtx, sci_x, 6, &pirow, &picol, &pdblDataID); CHECK_ERROR_API_SCILAB;
      One_Model_Input.DontCut = bool(pdblDataID[0]);

      _SciErr = getMatrixOfDoubleInList(pvApiCtx, sci_x, 7, &pirow, &picol, &pdblDataID); CHECK_ERROR_API_SCILAB;
      switch((int)*pdblDataID)
	{
	case 0:
	  One_Model_Input.Type = Input; // Input Output RetroInput Weight
	  break;
	case 1:
	  One_Model_Input.Type = Output; // Input Output RetroInput Weight
	  break;
	case 2:
	  One_Model_Input.Type = RetroInput; // Input Output RetroInput Weight
	  break;
	case 3:
	  One_Model_Input.Type = Weight; // Input Output RetroInput Weight
	  break;
	default:
	  Scierror(999,"%s: Wrong type for input or output type\n", fname);
	  return 0;
	}
      List_Model_Input.push_back(One_Model_Input);
    }

  // Get the cutlist
  int * p_cutlist_address;
  int * p_cutlist_item_address;
  int p_cutlist_item_nb_cols, p_cutlist_item_nb_rows, p_cutlist_nb_items;

  getVarAddressFromPosition(IOLIST_IN,&p_cutlist_address);

  res = getListItemNumber(p_cutlist_address, &p_cutlist_nb_items);

  List_PartToCut.resize(0);

  for(i=1;i<=p_cutlist_nb_items;i++)
    {
      res = getListItemAddress(p_cutlist_address, i, &p_cutlist_item_address);
      res = getVarDimension(p_cutlist_item_address, &p_cutlist_item_nb_rows, &p_cutlist_item_nb_cols);

      if (p_cutlist_item_nb_rows*p_cutlist_item_nb_cols!=5)
	{
	  sciprint("%s: ERROR: CUT usage - [NameOfPart DimensionNo CutValue LowPartName UpperPartName]\n", fname);
	  sciprint("%s: NameOfPart : name of the partition to cut\n", fname);
	  sciprint("%s: DimensionNo: index of the dimension to cut (starting from 0)\n", fname);
	  sciprint("%s: CutValue   : where to put the cut\n", fname);
	  sciprint("%s: LowPartName: name of the partition name (partition with the lowest values)\n", fname);
	  sciprint("%s: UpPartName : name of the partition name (partition with the largest values)\n", fname);
          Scierror(999,"%s: bad option usage\n", fname);

	  return 0;
	}

      _SciErr = getMatrixOfStringInList(pvApiCtx, sci_x, 2, &pirow, &picol, NULL, NULL);  CHECK_ERROR_API_SCILAB;
      pilen = (int *)MALLOC(pirow*picol*sizeof(int));
      _SciErr= getMatrixOfStringInList(pvApiCtx,sci_x,2,&pirow,&picol,pilen,NULL); CHECK_ERROR_API_SCILAB;
      pstStrings = (char **)MALLOC(pirow*picol*sizeof(char *));
      for(i=0;i<pirow*picol;i++)
	{
	  pstStrings[i] = (char *)MALLOC((pilen[i]+1)*sizeof(char));
	}
      _SciErr = getMatrixOfStringInList(pvApiCtx,sci_x,2,&pirow,&picol,pilen,pstStrings); CHECK_ERROR_API_SCILAB;
      Aux_PartToCut.PartName.str(pstStrings[0]);
      // YC: ajouter les free

      _SciErr = getMatrixOfDoubleInList(pvApiCtx, sci_x, 3, &pirow, &picol, &pdblDataID); CHECK_ERROR_API_SCILAB;
      Aux_PartToCut.Dimension = int(pdblDataID[0]);

      _SciErr = getMatrixOfDoubleInList(pvApiCtx, sci_x, 4, &pirow, &picol, &pdblDataID); CHECK_ERROR_API_SCILAB;
      Aux_PartToCut.CutPosition = int(pdblDataID[0]);

      _SciErr = getMatrixOfStringInList(pvApiCtx, sci_x, 5, &pirow, &picol, NULL, NULL);  CHECK_ERROR_API_SCILAB;
      pilen = (int *)MALLOC(pirow*picol*sizeof(int));
      _SciErr= getMatrixOfStringInList(pvApiCtx,sci_x,5,&pirow,&picol,pilen,NULL); CHECK_ERROR_API_SCILAB;
      pstStrings = (char **)MALLOC(pirow*picol*sizeof(char *));
      for(i=0;i<pirow*picol;i++)
	{
	  pstStrings[i] = (char *)MALLOC((pilen[i]+1)*sizeof(char));
	}
      _SciErr = getMatrixOfStringInList(pvApiCtx,sci_x,5,&pirow,&picol,pilen,pstStrings); CHECK_ERROR_API_SCILAB;
      Aux_PartToCut.PartName_LowerPart.str(pstStrings[0]);

      _SciErr = getMatrixOfStringInList(pvApiCtx, sci_x, 6, &pirow, &picol, NULL, NULL);  CHECK_ERROR_API_SCILAB;
      pilen = (int *)MALLOC(pirow*picol*sizeof(int));
      _SciErr= getMatrixOfStringInList(pvApiCtx,sci_x,6,&pirow,&picol,pilen,NULL); CHECK_ERROR_API_SCILAB;
      pstStrings = (char **)MALLOC(pirow*picol*sizeof(char *));
      for(i=0;i<pirow*picol;i++)
	{
	  pstStrings[i] = (char *)MALLOC((pilen[i]+1)*sizeof(char));
	}
      _SciErr = getMatrixOfStringInList(pvApiCtx,sci_x,6,&pirow,&picol,pilen,pstStrings); CHECK_ERROR_API_SCILAB;
      Aux_PartToCut.PartName_UpperPart.str(pstStrings[0]);

      List_PartToCut.push_back(Aux_PartToCut);
    }

  // Get the coefflist
  int * p_coefflist_address;
  int * p_coefflist_item_address;
  int p_coefflist_item_nb_cols, p_coefflist_item_nb_rows, p_coefflist_nb_items;

  getVarAddressFromPosition(COEFFLIST_IN,&p_coefflist_address);

  res = getListItemNumber(p_coefflist_address, &p_coefflist_nb_items);

  for(i=1;i<=p_coefflist_nb_items;i++)
    {
      res = getListItemAddress(p_coefflist_address, i, &p_coefflist_item_address);
      res = getVarDimension(p_coefflist_item_address, &p_coefflist_item_nb_rows, &p_coefflist_item_nb_cols);

      if ((p_coefflist_item_nb_rows*p_coefflist_item_nb_cols!=3)
	{
	  sciprint("ERROR: Coeff usage - [PartName Dimension Coeff]\n", fname);
	  sciprint("PartName  : name of the partition of which we change the coeff", fname);
	  sciprint("Dimension : index of the dimension to change the coeff (starting from 0)", fname);
	  sciprint("Coeff     : value of the coeffcient", fname);
          Scierror(999,"%s: bad option usage\n", fname);
	  return 0;
	}

      _SciErr = getMatrixOfStringInList(pvApiCtx, sci_x, 2, &pirow, &picol, NULL, NULL);  CHECK_ERROR_API_SCILAB;
      pilen = (int *)MALLOC(pirow*picol*sizeof(int));
      _SciErr= getMatrixOfStringInList(pvApiCtx,sci_x,2,&pirow,&picol,pilen,NULL); CHECK_ERROR_API_SCILAB;
      pstStrings = (char **)MALLOC(pirow*picol*sizeof(char *));
      for(i=0;i<pirow*picol;i++)
	{
	  pstStrings[i] = (char *)MALLOC((pilen[i]+1)*sizeof(char));
	}
      _SciErr = getMatrixOfStringInList(pvApiCtx,sci_x,2,&pirow,&picol,pilen,pstStrings); CHECK_ERROR_API_SCILAB;
      Aux_ModelCoeff.PartName.str(pstStrings[0]);

      _SciErr = getMatrixOfDoubleInList(pvApiCtx, sci_x, 3, &pirow, &picol, &pdblDataID); CHECK_ERROR_API_SCILAB;
      Aux_ModelCoeff.Dimension = int(pdblDataID[0]);

      _SciErr = getMatrixOfDoubleInList(pvApiCtx, sci_x, 4, &pirow, &picol, &pdblDataID); CHECK_ERROR_API_SCILAB;
      Aux_ModelCoeff.Coeff = pdblDataID[0];
    }

  // YC; on stocke les paramètres dans une plist
  int     tmp_int, tmp_res, * tmp_vec_int, tmp_size;
  double  tmp_double, * tmp_vec_double;
  char *  tmp_char;

  INIT_PARAM(PARAM_IN);

  // modelname option
  GET_PARAM_STRING("modelname", tmp_char, "lolimot", tmp_res);
  ModelName.str(tmp_char);

  // optimize sigma option
  GET_PARAM_VEC_DOUBLE("optimizesigma", tmp_vec_double, tmp_size, tmp_res);
  if (tmpsize==1) 
    {
      SigmaMin  = tmp_vec_double[0]; // Default value: 0.1
      SigmaStep = 0.01;
      SigmaMax  = SigmaMin + 10*SigmaStep;
    }
  else if (tmpsize==2)
    {
      SigmaMin  = tmp_vec_double[0]; // Default value: 0.1
      SigmaStep = tmp_vec_double[1]; // Default value: 0.01
      SigmaMax  = SigmaMin + 10*SigmaStep;
    }
  else if (tmpsize>=3)
    {
      SigmaMin  = tmp_vec_double[0]; // Default value: 0.1
      SigmaStep = tmp_vec_double[1]; // Default value: 0.01
      SigmaMax  = tmp_vec_double[2]; // Default value: SigmaMin + 10*SigmaStep
    }
  else
    {
      Scierror(999,"%s: Wrong number of parameters for the option %s\n",fname,"optimizesigma");
      return 0;
    }

  // nbmaxpart option
  GET_PARAM_DOUBLE("nbmaxpart", tmp_double, NBMAXPARTITIONS, tmp_res);
  NbMaxPart = tmp_double;

  // randseed option
  GET_PARAM_DOUBLE("randseed", tmp_double, 12345, tmp_res);
  RandSeed = tmp_double;

  // uselasso option
  GET_PARAM_DOUBLE("uselasso", tmp_double, 0, tmp_res);
  if (tmp_res!=-1)
    {
      if (tmp_double!=0)
	{
	  UseLasso = true;
	}
    }
  
  // residualtype option
  GET_PARAM_DOUBLE("residualtype", tmp_double, 0, tmp_res);
  if (tmp_res!=-1)
    {
      ResidualType = tmp_double;
    }

  // residualfilename option
  GET_PARAM_STRING("residualfilename", tmp_char, "data_residual.dat", tmp_res);
  Res_Filename.str(tmp_char);

  // r2filename option
  GET_PARAM_STRING("r2filename", tmp_char, "data_r2.dat", tmp_res);
  R2_Filename.str(tmp_char);

  // partlist option
  GET_PARAM_STRING("partlist", tmp_char, "data_partlist.dat", tmp_res);
  Part_Filename.str(tmp_char);
  if (tmp_res!=-1)
    {
      UsePartFilename = true;
    }

  // leaveoneout option
  GET_PARAM_STRING("leaveoneout", tmp_char, "data_leaveoneout.dat", tmp_res);
  LOO_Filename.str(tmp_char);
  if (tmp_res!=-1)
    {
      Use_LOO = true;
    }

  // membershipthreshold option
  GET_PARAM_DOUBLE("membershipthreshold", tmp_double, 0, tmp_res);
  if (tmp_res!=-1)
    {
      MembershipThreshold = tmp_double;
    }

  // learnstart option
  GET_PARAM_DOUBLE("learnstart", tmp_double, 0, tmp_res);
  Learn_Start = tmp_double;

  // validend option
  GET_PARAM_DOUBLE("learnend", tmp_double, 0, tmp_res);
  Learn_End = tmp_double;

  // shiftinputs option
  GET_PARAM_DOUBLE("shiftinputs", tmp_double, 0, tmp_res);
  ShiftInputs = tmp_double;
  
  // shiftoutputs option
  GET_PARAM_DOUBLE("shiftoutputs", tmp_double, 0, tmp_res);
  ShiftOutputs = tmp_double;

  // usetransformlog option
  GET_PARAM_VEC_DOUBLE("usetransformlog", tmp_vec_double, tmp_size, tmp_res);

  if (tmp_res!=-1)
    {
      UseTransformLog      = true;
      if (tmp_size==1)
	{
	  UseTransformLogAlpha = tmp_vec_double[0];
	  UseTransformLogEps   = 0.0;
	}
      else if (tmp_size==2)
	{
	  UseTransformLogAlpha = tmp_vec_double[0];
	  UseTransformLogEps   = tmp_vec_double[1];
	}
      else
	{
	  Scierror(999,"%s: Wrong number of parameters for option %s\n", fname, "usetransformlog");
	  return 0;
	}
    }

  // printpartitions option
  GET_PARAM_DOUBLE("printpartitions", tmp_double, 0, tmp_res);
  if (tmp_res!=-1) PrintPart = true;

  // usedistributedcutting option
  GET_PARAM_DOUBLE("usedistributedcutting", tmp_double, 0, tmp_res);
  if (tmp_res!=-1)
    {
      UseDistribCutting  = true;
      UseUniformCutting  = false;
    }

  // useuniformcutting option
  GET_PARAM_DOUBLE("useuniformcutting", tmp_double, 0, tmp_res);
  if (tmp_res!=-1)
    {
      UseDistribCutting  = false;
      UseUniformCutting  = true;
    }

  // display option
  GET_PARAM_DOUBLE("display", tmp_double, 0, tmp_res);
  if (tmp_res!=-1)
    {
      Display = true;
    }
  
  // resgapperc option
  GET_PARAM_DOUBLE("resgapperc", tmp_double, RESIDUGAPPERCENTAGE, tmp_res);
  ResiduGapPercentage = tmp_double;
     
  // sigma option
  GET_PARAM_DOUBLE("sigma", tmp_double, LL_SIGMA, tmp_res);
  Sigma = tmp_double;

  // nbcut option
  GET_PARAM_DOUBLE("nbcut", tmp_double, LL_NBCUT, tmp_res);
  NbCut = tmp_double;

  // measure option
  GET_PARAM_DOUBLE("measure", tmp_double, 0, tmp_res);
  if (tmp_res!=-1)
    {
      Measure = tmp_double;
    }

  // Initialisation of the random generator
  srand(RandSeed);

  // Temporisation are measure starting from 1 in the configuration file
  if (NbMaxTempo!=0) NbMaxTempo--;

  if (ModelName=="") ModelName = "Default";

  if (NbPartToExplore==-1)
    {
      NbPartToExplore=NbMaxPart;
    }

  if (Display)
    {
      sciprint("Model : %s\n",ModelName.str().c_str());
      
      sciprint("Learning data set - size: %d x %d\n", p_train_nb_rows, p_train_nb_cols);
      sciprint("Validation data set - size: %d x %d\n", p_valid_nb_rows, p_valid_nb_cols);
      
      for(i=0; i<List_Model_Input.size(); i++)
	{
	  sciprint("Input %d", i);
	  sciprint(" : Name = %s", List_Model_Input[i].Name.str().c_str());
	  sciprint(" Pos = %d", List_Model_Input[i].Pos);
	  sciprint(" Tempo = %d\n", List_Model_Input[i].Tempo);
	}
      
      for(i=0; i<List_Model_Output.size(); i++)
	{
	  sciprint("Output %d", i);
	  sciprint(" : Name = %s", List_Model_Output[i].Name.str().c_str());
	  sciprint(" Pos = %d", List_Model_Output[i].Pos);
	  sciprint(" Tempo = %d\n", List_Model_Output[i].Tempo);
	}
      
      sciprint("NbMaxPart           = %d\n", NbMaxPart);
      sciprint("ResiduGapPercentage = %d\n", ResiduGapPercentage);
      sciprint("Sigma               = %d\n", Sigma);
      sciprint("NbCut               = %d\n", NbCut);
      sciprint("Measure             = %d\n", Measure;
    }

  // Modification of the position of the columns
  for(i=1; i<List_Model_Input.size(); i++)
    {
      if (List_Model_Input[i].FileNumber!=0)
	List_Model_Input[i].Pos += Offset_For_Files[List_Model_Input[i].FileNumber-1];
    }

  for(i=0; i<List_Model_Input.size(); i++)
    {
      List_Model_Input[i].FileNumber = 0;
    }

  for(i=1; i<List_Model_Output.size(); i++)
    {
      if (List_Model_Output[i].FileNumber!=0)
	List_Model_Output[i].Pos += Offset_For_Files[List_Model_Output[i].FileNumber-1];
    }

  for(i=0; i<List_Model_Output.size(); i++)
    {
      List_Model_Output[i].FileNumber = 0;
    }

  if (List_Model_Input[0].FileNumber!=0)
    Weight_Input.Pos += Offset_For_Files[List_Model_Input[0].FileNumber-1];

  Weight_Input.FileNumber  = 0;

  // Generation of a fake MeasureData set
  MeasureData.resize(InputData.size(), vector<float>(List_Model_Output.size(), 0.0));

  if (Display) sciprint("Learning - Reading the data - Part 2 - generation of the offsets\n");

  // Opening the data files

  InputMin.resize(List_Model_Input.size());
  InputMax.resize(List_Model_Input.size());
  
  for(i=0; i<List_Model_Input.size(); i++)
    {
      InputMin[i] = numeric_limits<float>::max();
      InputMax[i] = - numeric_limits<float>::max();
    }
  
  MeasureMin.resize(List_Model_Output.size());
  MeasureMax.resize(List_Model_Output.size());
  
  for(i=0; i<List_Model_Output.size(); i++)
    {
      MeasureMin[i] = numeric_limits<float>::max();
      MeasureMax[i] = - numeric_limits<float>::max();
    }
  
  // Sizing of the input temporizations
  
  TempoInput.resize(List_Model_Input.size());
  
  for(i=0; i<List_Model_Input.size(); i++)
    {
      TempoInput[i].resize(List_Model_Input[i].Tempo + ShiftInput, 0.0);
    }

  // Sizing of the output temporizations
  
  TempoOutput.resize(List_Model_Output.size());
  
  for(i=0; i<List_Model_Output.size(); i++)
    {
      TempoOutput[i].resize(List_Model_Output[i].Tempo + ShiftOutput, 0.0);
    }

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
	    }
	  
	  InputData.push_back(DList);
	  
	  for(j=0; j<DList.size(); j++)
	    {
	      if (InputMin[j]>DList[j]) InputMin[j] = DList[j];
	      if (InputMax[j]<DList[j]) InputMax[j] = DList[j];
	    }
	  
	  // We get back the outputs
	  DList.resize(0);

	  for(j=0; j<List_Model_Output.size(); j++)
	    {
	      DList.push_back(Temporize(TempoOutput[j], AuxData[i][List_Model_Output[j].Pos]));
	    }
	      
	  MeasureData.push_back(DList);
	  
	  for(j=0; j<DList.size(); j++)
	    {
	      if (MeasureMin[j]>DList[j]) MeasureMin[j] = DList[j];
	      if (MeasureMax[j]<DList[j]) MeasureMax[j] = DList[j];
	    }

	  // We get back the weights
	  DList.resize(0);

	  for(j=0; j<WeightData.size(); j++)
	    {
	      WeightData.push_back(AuxData[i][Weight_Input.Pos]);
	    }
	}
      else if (i<NbMaxTempo)
	{
	  // We put the input data in the temporization line
	  for(j=0; j<List_Model_Input.size(); j++)
	    {
	      Temporize(TempoInput[j], AuxData[i][List_Model_Input[j].Pos]);
	    }
	  
	  // We put the output data in the temporization line
	  for(j=0; j<List_Model_Output.size(); j++)
	    {
	      Temporize(TempoOutput[j], AuxData[i][List_Model_Output[j].Pos]);
	    }
	}
    }

  AuxData.clear();

  if ((unsigned int)Measure>=MeasureData[0].size())
    {
      Scierror(999,"Measure must be included between 0 and %d\n", MeasureData[0].size());

      return 0;
    }

  Learn_InputMin   = InputMin;
  Learn_InputMax   = InputMax;
  Learn_MeasureMin = MeasureMin;
  Learn_MeasureMax = MeasureMax;

  if (Display) Display_DisplayData(InputData, MeasureData);
  
  if (Display) sciprint("Learning - Sampling the data\n");
  
  if (MeasureData.size()!=InputData.size())
    {
      Scierror(999,"Learning - Number of input lines is different than the number of output lines: %d - %d\n", InputData.size(), MeasureData.size());

      return 0;
    }

  if (Display)
    {
      Display_DisplayData(InputData, MeasureData);
    }

  if (Display) sciprint("Learning - Initialisation of the LOLIMOT network\n");
  
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
    }
    
  Lolimot->setUseDistributedCutting(UseDistribCutting);
  Lolimot->setUseUniformCutting(UseUniformCutting);
  
  if (CustomParam)      Lolimot->setCustomArguments(CustParamStr);
  if (CustomRet_Before) Lolimot->setCustomReturn_Before(CustRetBefStr);
  if (CustomRet_After)  Lolimot->setCustomReturn_After(CustRetAfterStr);
  
  // Manual cutting
  if (List_PartToCut.size()!=0)
    {
      if (Display)
	{
	  sciprint("Learning : Performing a manual cutting\n");
	  sciprint("Manual cutting: verification of the coherency of the cuttings\n");
	}
      
      for(i=0; i<List_PartToCut.size(); i++)
	{
	  if (InputMin[List_PartToCut[i].Dimension]>List_PartToCut[i].CutPosition)
	    {
	      sciprint("ERROR: the min value of the dimension %d", List_PartToCut[i].Dimension);
	      sciprint(" is greater than the position of the cutting.\n");
	      sciprint("Error on the cutting number %d\n", i);
	      
	      return 0;
	    }
	}
      
      for(i=0; i<List_Model_Input.size(); i++)
	{
	  Lolimot->addDimension(List_Model_Input[i].Name.c_str(), InputMin[i], InputMax[i], Sigma, NbCut);
	  if (List_Model_Input[i].DontCut) Lolimot->setDontCutVar(i, true);
	}
            
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
		}
	      
	      newPartition_Low->setDimensionMax(List_PartToCut[i].Dimension,
						List_PartToCut[i].CutPosition);
	      newPartition_Up->setDimensionMin(List_PartToCut[i].Dimension,
					       List_PartToCut[i].CutPosition);
	      newPartition_Low->setName(List_PartToCut[i].PartName_LowerPart);
	      newPartition_Up->setName(List_PartToCut[i].PartName_UpperPart);
	      
	      newPartition_Low->updateEcartType(Lolimot);
	      newPartition_Up->updateEcartType(Lolimot);
	    }
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
		    }
		}
	    }
	}
      
      // Updating the models
      Lolimot->updatePartitions();
      
      LL_Partition * AuxPart = NULL;
      
      if (Display)
	{
	  sciprint("List of manually created partitions\n");

	  for(i=0; i<Lolimot->getPartitionSet().size(); i++)
	    {
	      AuxPart = Lolimot->getPartition(i);
	      
	      sciprint("n %d %s\n", i+1, AuxPart->printDimensions());
	    }
	}
    }
  
  for(i=0; i<List_ModelCoeff.size(); i++)
    {
      if (Lolimot->getPartitionSet().size()==0)
	{
	  Scierror(999,"ERROR: There is no partitions !\n");

	  return 0;
	}
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
		    }
		  else
		    {
		      Lolimot->getPartition(j)->setFreezeCoeff(List_ModelCoeff[i].Dimension, true);
		      Lolimot->getPartition(j)->setCoeff(List_ModelCoeff[i].Dimension, List_ModelCoeff[i].Coeff);
		    }
		}
	    }
	}
    }
  
  if (PrintPart && Display)
    {
      sciprint("List of partitions\n");
      for(i=0; i<Lolimot->getPartitionSet().size(); i++)
	{
	  Lolimot->getPartition(i)->printDimensions();
	}
    }
  
  if (Display)
    {
      sciprint("Learning the LOLIMOT network\n");
    }

  if (Display)
    {
      sciprint("Learning - Classical method selected\n");
    }
  
  Learning_Classic(Lolimot, InputData, MeasureData, WeightData, List_Model_Input, InputMin, InputMax,
                   MeasureMin, MeasureMax, Sigma, NbCut, Measure, Load_Model);
  
  if (UsePartFilename)
    {
      if (Display)
	{
	  sciprint("Learning - Exportation of all the partitions\n");
	}
      
      Lolimot->exportAllPartitions(Part_Filename);
    }
  
  //
  // Cleaning the lists containing some data
  //
	  
  MeasureData.clear();
  InputData.clear();

  Lolimot->cleanDimensions();
  Lolimot->cleanMesures();
  Lolimot->cleanPartitions();
  Lolimot->cleanWeight();

  if (Lolimot) delete Lolimot;

  return 0;
}
