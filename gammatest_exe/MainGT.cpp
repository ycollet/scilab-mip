// MainGT.cpp
// Author: Y. Collette
// Date: 19/11/2005
// Renault - DR - 64240

#include <fstream>
#include <string>
#include <vector>
#include <ios>
#include <iostream>
#include <cstdlib>
#include <functional>
#include <limits>
#include <sstream>
#include <utility>

#include <math.h>

#ifdef WIN32
#include <kdtree_static.h>
#endif

#include <Arguments.h>
#include <GammaTest.h>
#include <GammaTest_Estimator.h>
#include <tokenize.hpp>
#include <Filter.h>

#define DEBUG

using namespace std;

enum Mask {MASK_0, MASK_1, MASK_X};

#define BUFFERSIZE 8048

int main(int argc, const char ** argv)
{
  vector<vector<double> >  InputData, Selected_InputData;
  vector<vector<double> >  MeasureData, Selected_MeasureData;
  vector<double>           Y_Output;
  vector<double>           List_Sensitivity;
  vector<double>           ListDelta;
  vector<double>           ListGamma;
  vector<double>           ListAllDelta;
  vector<double>           ListAllGamma, MTest_a, MTest_b;
  vector<double>           ListA, ListB, ListVR, ListR2;
  vector<double>           InputMin, InputMax, MeasureMin, MeasureMax;
  vector<pair<int, Mask> > InputMask;
  vector<string>           ListMask;
  string                   List_FileName, StringMask;
  string                   DataAllDG_FN, DataDG_FN, DataRes_FN, DataMTest_FN;
  unsigned int             List_SkipLines = 0;
  ifstream               * DataFile = NULL;
  ofstream               * OutFile = NULL;
  vector<string>           WList;
  vector<double>           DList;
  stringstream             SString;
  char                   * Buffer = NULL;
  double                   a = 0.0, b = 0.0, Proportion = 0.0;
  unsigned int             p = 20, CountMask = 0;
  int                      InputRemoved = 0;
  unsigned int             i, j, k, l, m, n;
  bool                     S_Select = false, S_Rand = false, S_Greedy = false, S_Greedy_Fast = false, UseUP = false;
  bool                     S_Distrib = false, S_Output = false, S_Input_Select_Exhaust = false, UseVersion2 = false;
  bool                     S_Min = true,  S_Input_Select_Greedy_Forward = false, S_Mask = false, S_NoGammaTest = false;
  bool                     S_MTest = false, S_Norm = false, S_Estimator = false, S_Input_Select_Greedy_Backward = false;
  bool                     S_BucketSize = false, S_Sensitivity = false;
  int                      S_Histo = 100, BucketSize = 5;
  int                      S_Nb_Select = 0, S_Nb_Rand = 0, S_Nb_Greedy = 0, S_Nb_Greedy_Fast = 0;
  int                      S_Nb_Distrib = 0, S_Nb_Output = 0;
  unsigned int             S_Start_Point = 0, S_End_Point = 0;
  unsigned int             Retenue = 0, Index = 0, IndexUP = 0, RandSeed = 0;
  int                      IndexMin, Starting_n = 0, MTest_Min = 0, MTest_Max = 0;
  double                   Min_a = numeric_limits<double>::max(), Min_b = numeric_limits<double>::max();
  double                   Min_VR = numeric_limits<double>::max();
  double                   Last_Min_a = numeric_limits<double>::max(), Last_Min_b = numeric_limits<double>::max();
  double                   Last_Min_VR = numeric_limits<double>::max();
  string                   Last_Min_Mask;
  double                   Variance = 0.0, Mean = 0.0, Max = - numeric_limits<double>::max();
  bool                     OneInputSelected = true;

  const char*      help_message =
    "Gamma Test Statistic \n"
    "Usage: MainGT [OPTIONS] -i [Input file name]\n"
    "Available Options:\n\n"
    "--help                        show this help.\n"
    "-i    [string]                Input file name (required).\n"
    "-oall [string]                Output file name where are stored all the computed points for the Gamma Test \n"
    "                              statistic (optional. Default value DataAllDG.dat).\n"
    "-oreg [string]                Output file name where are stored all the points which are used for the computation \n"
    "                              of the regression line (optional. Default value DataDG.dat).\n"
    "-ores [string]                Ouput file name where are stored the equation of the regression line (optional.\n"
    "                              Default value DataRes.dat).\n"
    "-omtest [string]              Output file name where are stored the results of the m-test (optional. Default\n"
    "                              value DataMTest.dat).\n"
    "-n [int value]                Number of neighborhood to be considered for the computation of the Gamma Test\n"
    "                              (default value 20).\n"
    "-p [float value]              Proportion of points to be saved into the file specified through the -oall option\n"
    "                              (default value 1.0). This value must reside in the 0 - 1 interval.\n"
    "-s_output_nb   [int value]    Select nb points w.r.t the output if -s_output_min has a 1 parameter then it selects\n"
    "                              the nb minima values w.r.t the output.\n"
    "-s_output_min  [bool value]   Defines if we select the minima or maxima value w.r.t the output (default value 1.\n"
    "                              Minima values are selected\n"
    "-s_greedy      [int value]    Select nb points using a greedy method w.r.t the input parameters.\n"
    "-s_greedy_fast [int value]    Select nb points using a fast greedy method w.r.t the input parameters.\n"
    "-s_norm                       Normalize input and output data.\n"
    "-s_rand        [int value]    Select nb points at random without duplicated points.\n"
    "-s_randseed    [int value]    Set the rand seed value.\n"
    "-s_distrib_nb  [int value]    Select nb points at random using the inputs parameter values as a probability distribution.\n"
    "-s_distrib_histo [int value]  Use a histogram with nb bar to compute the probability distribution (default value 100).\n"
    "-s_select      [int value]    Select one point each nb steps.\n"
    "-s_input_select_exhaust       Perform the selection of inputs via the Gamma Test (may be long to compute).\n"
    "-s_input_select_greedy_forward  Perform the selection of inputs via the Gamma Test using a greedy method.\n"
    "                                The input selection process starts with no inputs and add one input\n"
    "                                after the other.\n"
    "-s_input_select_greedy_backward Perform the selection of inputs via the Gamma Test using a greedy method.\n"
    "                                The input selection process starts with all inputs and remove one input\n"
    "                                after the other.\n"
    "-s_mtest_min [int value]      \n"
    "-s_mtest_max [int value]      Perform the M-test between min and max. If just min or max set then min=max.\n"
    "-no_gamma_test                Skip the computation of the gamma test (if you just want to perform input selection.\n"
    "-use_version_2                Use the moving average for the computation of the Gamma Test.\n"
    "-use_unique_point [int value] Use the neighborhood of a single point to compute the Gamma Test.\n"
    "-starting_n [int value]       Start computing the Gamma Test after the nth neighbor.\n"
    "-input_mask    [string]       InputMask[i] = 0 - the input is always not selected.\n"
    "                              InputMask[i] = 1 - the input is always selected.\n"
    "                              InputMask[i] = X - the input is tested by the method.\n"
    "-s_estimator                  Compute the estimated value of y for each y using gamma test.\n"
    "-s_bucket_size [int]          Set the bucket size of the KD-tree.\n"
    "-s_sensitivity                Perform sensitivity test on the data set.\n"
    "-s_start_point [int]          First point to read in the data file.\n"
    "-s_end_point   [int]          Last point to read in the data file.\n\n";
  
  Arguments::setArgumentsWithSpaces("-i|-oall|-ores|-oreg|-omtest|-p|-n|-s_output_nb|-s_output_min|-s_greedy|-s_greedy_fast|-s_rand|-s_distrib_nb|-s_distrib_histo|-s_select|-s_input_select|-use_unique_point|-s_mtest_min|-s_mtest_max|-s_randseed|-s_bucket_size|-s_start_point|-s_end_point|-input_mask");
  
  const Arguments Args(argc, argv);
    
  if (Args.size() == 0)
    {
      cerr << "Error: No arguments received.\n"
	   << "See 'MainGT --help' for options." << endl;    
      exit(1);
    } /* End If */

  if (Args.has("-h|--help"))
    {
      cout << help_message << endl;      
      exit(0);
    } /* End If */

  //
  // Ici, il faut que InputFilename soit du type string
  //

  if (!Args.get("-i", List_FileName, " "))	
    {
      cerr << "Error: -i argument is missing or invalid.  "
	   << "(-i requires a file name.)" << endl;
      exit(1);
    } /* End If */

  Args.get("-oall", DataAllDG_FN, "DataAllDG.dat");
  Args.get("-oreg", DataDG_FN, "DataDG.dat");
  Args.get("-ores", DataRes_FN, "DataRes.dat");
  Args.get("-omtest", DataMTest_FN, "DataMTest.dat");
  Args.get("-n", p, 20);
  Args.get("-p", Proportion, 1.0);

  S_Select = Args.get("-s_select", S_Nb_Select, 1);
  S_Output = Args.get("-s_output_nb", S_Nb_Output, 1);
  S_Min    = Args.get("-s_output_min", S_Min, true);
  S_Greedy = Args.get("-s_greedy", S_Nb_Greedy, 1);
  S_Greedy_Fast = Args.get("-s_greedy_fast", S_Nb_Greedy_Fast, 1);
  S_Rand    = Args.get("-s_rand", S_Nb_Rand, 1);
  S_Distrib = Args.get("-s_distrib_nb", S_Nb_Distrib, 1);
  Args.get("-s_distrib_histo", S_Histo, 100);
  S_Input_Select_Exhaust = Args.has("-s_input_select_exhaust");
  S_Input_Select_Greedy_Forward  = Args.has("-s_input_select_greedy_forward");
  S_Input_Select_Greedy_Backward = Args.has("-s_input_select_greedy_backward");
  S_Mask        = Args.get("-input_mask", StringMask, "");
  S_NoGammaTest = Args.has("-no_gamma_test");
  UseVersion2   = Args.has("-use_version_2");
  UseUP         = Args.get("-use_unique_point", IndexUP, 0);
  Args.get("-starting_n", Starting_n, 1);
  S_MTest       = Args.get("-s_mtest_min", MTest_Min, -1);
  S_Norm        = Args.has("-s_norm");
  Args.get("-s_randseed", RandSeed, 0);
  S_Estimator   = Args.has("-s_estimator");
  S_BucketSize  = Args.get("-s_bucket_size", BucketSize, 5);
  S_Sensitivity = Args.has("-s_sensitivity");
  Args.get("-s_start_point", S_Start_Point, 0);
  Args.get("-s_end_point", S_End_Point, 0);
#ifdef DEBUG
  if (S_Mask) cout << "S_Mask OK" << endl;
  else        cout << "S_Mask not set" << endl;
  cout << "StringMask = " << StringMask << endl;
#endif
  srand(RandSeed);

  if (S_MTest)
    {
      Args.get("-s_mtest_max", MTest_Max, -1);
    } /* End If */

  if ((MTest_Max==-1)&&(MTest_Min==-1)) 
    {
      MTest_Max = p;
      MTest_Min = p;
    } /* End If */
  else if ((MTest_Max==-1)&&(MTest_Min!=-1))
    {
      MTest_Max = MTest_Min;
    } /* End Else If */
  else if ((MTest_Max!=-1)&&(MTest_Min==-1))
    {
      MTest_Min = MTest_Max;
    } /* End Else If */

  MeasureData.resize(0);
  InputData.resize(0);

  ListDelta.resize(0);
  ListGamma.resize(0);
  ListAllDelta.resize(0);
  ListAllGamma.resize(0);

  WList.resize(0);

  InputMin.resize(0);
  InputMax.resize(0);
  MeasureMin.resize(0);
  MeasureMax.resize(0);

  InputMask.resize(0);
  ListA.resize(0);
  ListB.resize(0);
  ListVR.resize(0);
  ListR2.resize(0);
  ListMask.resize(0);

  MTest_a.resize(0);
  MTest_b.resize(0);

  // Ouvertures des fichiers

  cout << "Ouverture du fichier " << List_FileName << "." << endl;
      
  DataFile = new ifstream;

  DataFile->open(List_FileName.c_str());
      
  if (!DataFile->is_open())
    {
      cerr << "Problème lors de l'ouverture du fichier " << List_FileName << "." << endl;
      exit(1);
    } /* End If */
	  
  Buffer = new char[BUFFERSIZE];
      
  // On passe les lignes contenant les intitulés des colonnes
  for(i=0; i<List_SkipLines; i++)
    {
      DataFile->getline(Buffer, BUFFERSIZE, '\n');
    } /* End For */

  while(DataFile->getline(Buffer, BUFFERSIZE, '\n'))
    {
      WList = tokenize(Buffer, is_space());

      // Si on est arrivé à la fin d'un fichier, on arrète
      // On post-traite InputData et MeasureData: on supprime les NbMaxTempo premières lignes
      // pour éviter de créer des données abhérantes.
      
      // On récupère les entrées
      DList.resize(0);
      
      for(j=0; j<WList.size()-1; j++)
	{
	  DList.push_back(atof(WList[j].c_str()));
	} /* End For */
      
      InputData.push_back(DList);
      
      if (InputMin.size()==0)
	{
	  InputMin.resize(InputData[0].size(), numeric_limits<double>::max());
	  InputMax.resize(InputData[0].size(), numeric_limits<double>::min());
	  MeasureMin.resize(1, numeric_limits<double>::max());
	  MeasureMax.resize(1, numeric_limits<double>::min());
	} /* End If */

      for(i=0; i<InputMin.size(); i++)
	{
	  if (InputMin[i]>DList[i]) InputMin[i] = DList[i];
	  if (InputMax[i]<DList[i]) InputMax[i] = DList[i];
	} /* End For */

      DList.resize(0);

      DList.push_back(atof(WList[WList.size()-1].c_str()));
      MeasureData.push_back(DList);

      if (MeasureMin[0]>DList[0]) MeasureMin[0] = DList[0];
      if (MeasureMax[0]<DList[0]) MeasureMax[0] = DList[0];
    } /* End While */

  DataFile->close();

  if (Buffer)   delete [] Buffer;
  if (DataFile) delete DataFile;

  // Retrieving the S_Start_Point - S_End_Point interval
  
  if ((S_Start_Point==0)&&(S_End_Point==0))
    {
      S_Start_Point = 0;
      S_End_Point   = InputData.size();
    } /* End If */
  else if (S_End_Point==0)
    {
      S_End_Point = InputData.size();
    } /* End Else If */
  else if (S_End_Point<S_Start_Point)
    {
      cout << "Error: S_Start_Point must be smaller than S_End_Point" << endl;
      exit(1);
    }  /* End Else If */

  Selected_InputData.resize(0);
  Selected_MeasureData.resize(0);

  for(i=S_Start_Point; i<S_End_Point; i++)
    {
      Selected_InputData.push_back(InputData[i]);
      Selected_MeasureData.push_back(MeasureData[i]);
    } /* End For */

  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;

  Selected_InputData.resize(0);
  Selected_MeasureData.resize(0);

  // Normalisation

  if (S_Norm)
    {
      cout << "Normalisation of the data set" << endl;

      for(i=0; i<InputData.size(); i++)
	{
	  for(j=0; j<InputData[i].size(); j++)
	    {
	      InputData[i][j] = (InputData[i][j] - InputMin[j])/(double)(InputMax[j] - InputMin[j]);
	    } /* End For */
	  MeasureData[i][0] = (MeasureData[i][0] - MeasureMin[0])/(double)(MeasureMax[0] - MeasureMin[0]);
	} /* End For */

      for(i=0; i<InputMin.size(); i++)
	{
	  InputMin[i] = 0;
	  InputMax[i] = 1;
	} /* End For */

      MeasureMin[0] = 0;
      MeasureMax[0] = 1;
    } /* End If */

  // Traitement du masque des entrées
  if (S_Mask)
    {
      if (StringMask.length()!=InputData[0].size())
	{
	  cerr << "The size of the mask if different from the number of inputs" << endl;
	  cerr << "Mask size = " << StringMask.length() << endl;
	  cerr << "Number of inputs = " << InputData[0].size() << endl;
	  exit(1);
	} /* End If */

      cout << "Input Mask Value : ";
      InputMask.resize(StringMask.length());
      for(i=0; i<InputMask.size(); i++)
	{
	  cout << StringMask[i];
	  if (StringMask[i]=='0')
	    {
	      InputMask[i].first = 0;
	      InputMask[i].second = MASK_0;
	    } /* End For */
	  else if (StringMask[i]=='1')
	    {
	      InputMask[i].first = 1;
	      InputMask[i].second = MASK_1;
	    } /* End If */
	  else
	    {
	      InputMask[i].first = 0;
	      InputMask[i].second = MASK_X;
	    } /* End Else */
	} /* End For */
      cout << endl;
    } /* End If */
  else
    {
      InputMask.resize(InputData[0].size());
      for(i=0; i<InputMask.size(); i++)
	{
	  InputMask[i].first  = 0;
	  InputMask[i].second = MASK_X;
	} /* End For */
    } /* End Else */

  if (S_Select)      Filter_Select(InputData, MeasureData, S_Nb_Select);
  if (S_Output)      Filter_Echant_Output(InputData, MeasureData, MeasureMin, MeasureMax, 0, S_Nb_Output, S_Min);
  if (S_Greedy)      Filter_EchantGreedy(InputData, MeasureData, InputMin, InputMax, S_Nb_Greedy);
  if (S_Greedy_Fast) Filter_EchantGreedy_Fast(InputData, MeasureData, InputMin, InputMax,
					      (InputData.size() - S_Nb_Greedy_Fast)/100, S_Nb_Greedy_Fast);
  if (S_Rand)        Filter_EchantRand(InputData, MeasureData, S_Nb_Rand);
  if (S_Distrib)     Filter_EchantStat(InputData, MeasureData, S_Histo, S_Nb_Distrib);

  // Computation of the Variance Ratio

  for(i=0; i<MeasureData.size(); i++)
    {
      Mean += MeasureData[i][0];
    } /* End For */

  Mean /= (double)MeasureData.size();

  for(i=0; i<MeasureData.size(); i++)
    {
      Variance += (MeasureData[i][0] - Mean) * (MeasureData[i][0] - Mean);
    } /* End For */

  Variance /= (double)MeasureData.size();

  if (!S_NoGammaTest)
    {
      GammaTest(InputData, MeasureData, ListDelta, ListGamma, ListAllDelta, ListAllGamma,
		p, a, b, Proportion, 0, UseVersion2, UseUP, IndexUP, Starting_n, S_MTest, MTest_Min, MTest_Max,
		MTest_a, MTest_b, BucketSize);
      
      if (S_MTest)
	{
	  OutFile = new ofstream;

	  OutFile->open(DataMTest_FN.c_str());

	  (*OutFile) << "% p - Delta - Gamma " << endl;
	  for(i=0, j=MTest_Min; i<MTest_a.size(); i++, j++)
	    {
	      (*OutFile) << i << " " << MTest_a[i] << " " << MTest_b[i] << endl;
	    } /* End For */

	  OutFile->close();

	  delete OutFile;
	} /* End If */

      OutFile = new ofstream;
      
      OutFile->open(DataDG_FN.c_str());
      
      (*OutFile) << "% Delta - Gamma " << endl;
      for(i=0; i<ListDelta.size(); i++)
	{
	  (*OutFile) << ListDelta[i] << " " << ListGamma[i] << endl;
	} /* End For */
      
      OutFile->close();
      
      if (OutFile) delete OutFile;
      
      OutFile = new ofstream;
      
      OutFile->open(DataAllDG_FN.c_str());
      
      (*OutFile) << "% Delta - Gamma" << endl;
      for(i=0; i<ListAllDelta.size(); i++)
	{
	  (*OutFile) << ListAllDelta[i] << " " << ListAllGamma[i] << endl;
	} /* End For */
      
      OutFile->close();
      
      if (OutFile) delete OutFile;
      
      OutFile = new ofstream;
      
      OutFile->open(DataRes_FN.c_str());
      
      (*OutFile) << "Gamma = G * Delta + A" << endl;
      (*OutFile) << "G = " << b << " A = " << a << endl;
      (*OutFile) << "R2 = 1 - abs(A)/Var(Output) = " << 1 - fabs(a)/Variance << endl;
      (*OutFile) << "VR = abs(A)/Var(Output) = " << fabs(a)/Variance << endl;

      OutFile->close();
      
      if (OutFile) delete OutFile;
    } /* End If */
  
  if (S_Input_Select_Exhaust)
    {
      ListA.resize(0);
      ListB.resize(0);
      ListVR.resize(0);
      ListR2.resize(0);
      ListMask.resize(0);

      // Initialisation of the mask

      InputMask.resize(InputData[0].size());

      if (!S_Mask)
	{
	  for(i=0; i<InputMask.size(); i++)
	    {
	      InputMask[i].first = 0;
	      InputMask[i].second = MASK_X;
	    } /* End For */
	} /* End If */

      // Selection of inputs
      
      cout << "Selection of Inputs" << endl;
      
      Retenue = 0;
      Index   = 0;

      while(!Retenue)
	{
	  while((InputMask[Index].second!=MASK_X)&&(Index<InputMask.size())) Index++;

	  InputMask[Index].first += 1;
	  if (InputMask[Index].first>1)
	    {
	      InputMask[Index].first = 0;
	      Retenue = 1;
	    } /* End If */
	  else
	    {
	      Retenue = 0;
	    } /* End Else */
	  
	  for(i=1; i<InputMask.size(); i++)
	    {
	      if (InputMask[i].second!=MASK_X) continue;

	      InputMask[i].first += Retenue;
	      if (InputMask[i].first>1)
		{
		  InputMask[i].first = 0;
		  Retenue = 1;
		} /* End If */
	      else
		{
		  Retenue = 0;
		} /* End Else */
	    } /* End For i */
	
	  if (Retenue) break;

	  cout << "Mask value : ";

	  SString.str("");

	  for(i=0; i<InputMask.size(); i++)
	    {
	      cout << InputMask[i].first;
	      SString << InputMask[i].first;
	    } /* End For */
	  cout << endl;

	  for(j=0; j<InputMask.size(); j++)
	    {
	      switch(InputMask[j].second)
		{
		case MASK_0:
		  cout << 0;
		  break;
		case MASK_1:
		  cout << 1;
		  break;
		case MASK_X:
		  cout << "X";
		  break;
		} /* End Switch */
	    } /* End For */
	  cout << endl;
	  
	  // Si a ce moment, la retenue vaut 1 -> on met le
	  // drapeau de retenue a 1 de la variable correspondante
	  
	  cout << "Preparing the data with selected inputs" << endl;

	  Selected_InputData.resize(0);
	  
	  for(i=0; i<InputData.size(); i++)
	    {
	      DList.resize(0);
	      for(j=0; j<InputMask.size(); j++)
		{
		  if (InputMask[j].first==1)
		    {
		      DList.push_back(InputData[i][j]);
		    } /* End If */
		} /* End For */
	      Selected_InputData.push_back(DList);
	    } /* End For */
	  
	  cout << "Performing the Gamma test" << endl;
	  
	  ListDelta.resize(0);
	  ListGamma.resize(0);
	  ListAllDelta.resize(0);
	  ListAllGamma.resize(0);

	  GammaTest(Selected_InputData, MeasureData, ListDelta, ListGamma, ListAllDelta, ListAllGamma,
		    p, a, b, Proportion, 0, UseVersion2, UseUP, IndexUP, Starting_n, false, MTest_Min, MTest_Max,
		    MTest_a, MTest_b, BucketSize);
	  
	  cout << "Value of A (Variance of the modelisation noise) : " << a << endl;
	  cout << "Value of B                                      : " << b << endl;
	  cout << "Value of the Variance ratio                     : " << fabs(a)/Variance     << endl;
	  cout << "Value of the R2                                 : " << 1 - fabs(a)/Variance << endl;

	  ListA.push_back(a);
	  ListB.push_back(b);
	  ListVR.push_back(fabs(a)/Variance);
	  ListR2.push_back(1 - fabs(a)/Variance);
	  ListMask.push_back(SString.str());
	} /* End While */
      
      // Computation of the sensitivity of the inputs

      cout << "Computing the summed-up information related to the selected inputs" << endl;
      
      Selected_InputData.resize(0);

      Max = - numeric_limits<double>::max();

      for(m=0, n=0; m<InputMask.size(); m++, n++)
	{
	  if (InputMask[m].first==0) continue;

	  for(k=0; k<InputData.size(); k++)
	    {
	      DList.resize(0);
	      for(l=0; l<InputMask.size(); l++)
		{
		  if ((InputMask[l].first==1)&&(l!=m))
		    {
		      DList.push_back(InputData[k][l]);
		    } /* End If */
		} /* End For */
	      Selected_InputData.push_back(DList);
	    } /* End For */
	  
	  cout << "Performing the Gamma test without input " << n << endl;
	  
	  ListDelta.resize(0);
	  ListGamma.resize(0);
	  ListAllDelta.resize(0);
	  ListAllGamma.resize(0);
	  
	  GammaTest(Selected_InputData, MeasureData, ListDelta, ListGamma, ListAllDelta, ListAllGamma,
		    p, a, b, Proportion, 0, UseVersion2, UseUP, IndexUP, Starting_n, false, MTest_Min, MTest_Max,
		    MTest_a, MTest_b, BucketSize);
	  
	  List_Sensitivity.push_back(1-fabs(a)/Variance);

	  if (Max<(1-fabs(a)/Variance)) Max = 1-fabs(a)/Variance;
	} /* End For */

      OutFile = new ofstream;
      
      OutFile->open("DataSelectInputsAB.dat");
      
      (*OutFile) << "% Y = B*X+A" << endl;
      (*OutFile) << "% Mask -- B -- A -- R -- R2" << endl;

      for(i=0; i<ListA.size(); i++)
	{
	  (*OutFile) << ListMask[i] << " -- " << ListB[i] << " -- " << ListA[i] << " -- ";
	  (*OutFile) << ListVR[i]   << " -- " << ListR2[i] << endl;
	} /* End For */
      
      if (List_Sensitivity.size()!=0)
	{
	  cout << "Sensitivity of the output with respect to the inputs: " << endl;

	  for(m=0; m<List_Sensitivity.size(); m++)
	    {
	      cout << "Input " << m << " removed : decrease of R2 = " << (List_Sensitivity[m] - Max)/Max * 100.0 << endl;
	    } /* End For */
	} /* End If */

      OutFile->close();
      
      delete OutFile;
    } /* End If */

  if (S_Input_Select_Greedy_Forward)
    {
      OneInputSelected = false;
      Last_Min_a  = numeric_limits<double>::max();
      Last_Min_b  = numeric_limits<double>::max();
      Last_Min_VR = numeric_limits<double>::max();
      Last_Min_Mask = "";

      ListA.resize(0);
      ListB.resize(0);
      ListVR.resize(0);
      ListR2.resize(0);
      ListMask.resize(0);

      // Initialisation of the mask

      CountMask = 0;

      for(i=0; i<InputMask.size(); i++)
	{
	  if (InputMask[i].second==MASK_X)
	    {
	      InputMask[i].first = 0;
	      CountMask++;
	    } /* End If */
	} /* End For */

      // Selection of inputs
      
      cout << "Selection of Inputs using a forward greedy method" << endl;
      
      for(i=0; i<InputMask.size(); i++)
	{
	  IndexMin = -1;
	  Min_a  = numeric_limits<double>::max();
	  Min_b  = numeric_limits<double>::max();
	  Min_VR = numeric_limits<double>::max();

	  for(j=0; j<InputMask.size(); j++)
	    {
	      if (InputMask[j].second!=MASK_X) continue;

 	      InputMask[j].first = 1;

	      SString.str("");
	      
	      cout << "Mask value : ";
	      for(k=0; k<InputMask.size(); k++)
		{
		  cout << InputMask[k].first;
		  SString << InputMask[k].first;
		} /* End For */
	      cout << endl;
	      
	      cout << "Mask value : ";
	      for(k=0; k<InputMask.size(); k++)
		{
		  switch(InputMask[k].second)
		    {
		    case MASK_0:
		      cout << 0;
		      break;
		    case MASK_1:
		      cout << 1;
		      break;
		    case MASK_X:
		      cout << "X";
		      break;
		    } /* End Switch */
		} /* End For */
	      cout << endl;
	      
	      // Si a ce moment, la retenue vaut 1 -> on met le
	      // drapeau de retenue a 1 de la variable correspondante
	      
	      cout << "Preparing the data with selected inputs" << endl;

	      Selected_InputData.resize(0);
	      
	      for(k=0; k<InputData.size(); k++)
		{
		  DList.resize(0);
		  for(l=0; l<InputMask.size(); l++)
		    {
		      if (InputMask[l].first==1)
			{
			  DList.push_back(InputData[k][l]);
			} /* End If */
		    } /* End For */
		  Selected_InputData.push_back(DList);
		} /* End For */

	      cout << "Performing the Gamma test" << endl;

	      ListDelta.resize(0);
	      ListGamma.resize(0);
	      ListAllDelta.resize(0);
	      ListAllGamma.resize(0);

	      GammaTest(Selected_InputData, MeasureData, ListDelta, ListGamma, ListAllDelta, ListAllGamma,
			p, a, b, Proportion, 0, UseVersion2, UseUP, IndexUP, Starting_n, false, MTest_Min, MTest_Max,
			MTest_a, MTest_b, BucketSize);
	      
	      cout << "Value of A (Variance of the modelisation noise) : " << a << endl;
	      cout << "Value of B                                      : " << b << endl;
	      cout << "Value of the variance ratio                     : " << fabs(a)/Variance << endl;

	      ListA.push_back(a);
	      ListB.push_back(b);
	      ListVR.push_back(fabs(a)/Variance);
	      ListR2.push_back(1 - fabs(a)/Variance);
	      ListMask.push_back(SString.str());

	      if (fabs(a)/Variance<Min_VR)
		{
		  Min_a    = a;
		  Min_b    = b;
		  Min_VR   = fabs(a)/Variance;
		  IndexMin = j;
		  OneInputSelected = true;
		} /* End If */

	      if (Last_Min_VR>Min_VR)
		{
		  Last_Min_VR = Min_VR;
		  Last_Min_a  = Min_a;
		  Last_Min_b  = Min_b;
		  Last_Min_Mask = SString.str();
		} /* End If */

	      // After having removed one  input we set all MASK_X inputs to 0

	      for(k=0; k<InputMask.size(); k++)
		{
		  if (InputMask[k].second==MASK_X) InputMask[k].first = 0;
		} /* End For */
	    } /* End For */
	  if (IndexMin!=-1)
	    {
	      InputMask[IndexMin].first  = 1;
	      InputMask[IndexMin].second = MASK_1;
	    } /* End If */
	  else 
	    {
	      // No inputs have been selected: we stop the selection process
	      break;
	    } /* End Else */
	} /* End For */
      
      // Computation of the sensitivity of the inputs
      
      List_Sensitivity.resize(0);

      Max = - numeric_limits<double>::max();

      vector<unsigned int> ListOfSelectedInputs(0);

      for(m=0; m<InputMask.size(); m++)
	{
	  if (Last_Min_Mask[m]=='1') ListOfSelectedInputs.push_back(m);
	} /* End For */

      if (ListOfSelectedInputs.size()>1)
	{
	  cout << "Computing the summed-up information related to the selected inputs" << endl;

	  for(m=0; m<ListOfSelectedInputs.size(); m++)
	    {
	      Selected_InputData.resize(0);
	      
	      for(k=0; k<InputData.size(); k++)
		{
		  DList.resize(0);
		  for(l=0; l<ListOfSelectedInputs.size(); l++)
		    {
		      if (l!=m)
			{
			  DList.push_back(InputData[k][ListOfSelectedInputs[l]]);
			} /* End If */
		    } /* End For */
		  Selected_InputData.push_back(DList);
		} /* End For */
	      
	      cout << "Performing the Gamma test without the selected input n° " << m << endl;
	      
	      ListDelta.resize(0);
	      ListGamma.resize(0);
	      ListAllDelta.resize(0);
	      ListAllGamma.resize(0);
	      
	      GammaTest(Selected_InputData, MeasureData, ListDelta, ListGamma, ListAllDelta, ListAllGamma,
			p, a, b, Proportion, 0, UseVersion2, UseUP, IndexUP, Starting_n, false, MTest_Min, MTest_Max,
			MTest_a, MTest_b, BucketSize);
	      
	      List_Sensitivity.push_back(1-fabs(a)/Variance);
	      
	      if (Max<(1-fabs(a)/Variance)) Max = 1-fabs(a)/Variance;
	    } /* End For */
	} /* End If */
      
      OutFile = new ofstream;
      
      OutFile->open("DataSelectInputsAB.dat");

      if (OneInputSelected)
	{
	  (*OutFile) << "Best selection :" << endl;
	  (*OutFile) << "Best Mask : " << Last_Min_Mask << endl;
	  (*OutFile) << "Y=B*X+A" << endl;
	  (*OutFile) << "B -- A -- VR -- R2" << endl;
	  (*OutFile) << Last_Min_b << " -- " << Last_Min_a << " -- " << Last_Min_VR << " -- " << 1 - Last_Min_VR << endl;
	} /* End If */
      else
	{
	  (*OutFile) << "Best selection :" << endl;
	  (*OutFile) << "No input selected !! " << endl;
	} /* End Else */

      (*OutFile) << "Mask value : ";
      for(j=0; j<InputMask.size(); j++)
	{
	  (*OutFile) << InputMask[j].first;
	} /* End For */
      (*OutFile) << endl;
      
      for(j=0; j<InputMask.size(); j++)
	{
	  switch(InputMask[j].second)
	    {
	    case MASK_0:
	      (*OutFile) << 0;
	      break;
	    case MASK_1:
	      (*OutFile) << 1;
	      break;
	    case MASK_X:
	      (*OutFile) << "X";
	      break;
	    } /* End Switch */
	} /* End For */
      (*OutFile) << endl;

      (*OutFile) << "% Y = B*X+A" << endl;
      (*OutFile) << "% Mask -- B -- A -- VR -- R2" << endl;

      for(i=0; i<ListA.size(); i++)
	{
	  (*OutFile) << ListMask[i] << " -- " << ListB[i] << " -- " << ListA[i] << " -- ";
	  (*OutFile) << ListVR[i]   << " -- " << ListR2[i] << endl;
	} /* End For */

      if (ListOfSelectedInputs.size()>1)
	{
	  if (List_Sensitivity.size()!=0)
	    {
	      cout << "Sensitivity of the output with respect to the inputs: " << endl;
	      
	      for(m=0; m<List_Sensitivity.size(); m++)
		{
		  cout << "Input " << m << " removed : decrease of R2 = " << (List_Sensitivity[m] - Max)/Max * 100.0 << endl;
		} /* End For */
	    } /* End If */
	} /* End If */

      OutFile->close();

      delete OutFile;
    } /* End If */

  if (S_Input_Select_Greedy_Backward)
    {
      OneInputSelected = false;
      Last_Min_a  = numeric_limits<double>::max();
      Last_Min_b  = numeric_limits<double>::max();
      Last_Min_VR = numeric_limits<double>::max();
      Last_Min_Mask = "";

      ListA.resize(0);
      ListB.resize(0);
      ListVR.resize(0);
      ListR2.resize(0);
      ListMask.resize(0);

      // Initialisation of the mask

      CountMask = 0;

      for(i=0; i<InputMask.size(); i++)
	{
	  if (InputMask[i].second==MASK_X)
	    {
	      InputMask[i].first = 1;
	      CountMask++;
	    } /* End If */
	} /* End For */

      // Selection of inputs
      
      cout << "Selection of Inputs using a backward greedy method" << endl;
      
      // Computation of the Gamma Test values for the whole input data set

      SString.str("");
      
      cout << "Mask value : ";
      for(k=0; k<InputMask.size(); k++)
	{
	  cout << InputMask[k].first;
	  SString << InputMask[k].first;
	} /* End For */
      cout << endl;
      
      cout << "Mask value : ";
      for(k=0; k<InputMask.size(); k++)
	{
	  switch(InputMask[k].second)
	    {
	    case MASK_0:
	      cout << 0;
	      break;
	    case MASK_1:
	      cout << 1;
	      break;
	    case MASK_X:
	      cout << "X";
	      break;
	    } /* End Switch */
	} /* End For */
      cout << endl;

      cout << "Preparing the data with selected inputs" << endl;
      
      Selected_InputData.resize(0);
      
      for(k=0; k<InputData.size(); k++)
	{
	  DList.resize(0);
	  for(l=0; l<InputMask.size(); l++)
	    {
	      if (InputMask[l].first==1)
		{
		  DList.push_back(InputData[k][l]);
		} /* End If */
	    } /* End For */
	  Selected_InputData.push_back(DList);
	} /* End For */
      
      cout << "Performing the Gamma test" << endl;
      
      ListDelta.resize(0);
      ListGamma.resize(0);
      ListAllDelta.resize(0);
      ListAllGamma.resize(0);
      
      GammaTest(Selected_InputData, MeasureData, ListDelta, ListGamma, ListAllDelta, ListAllGamma,
		p, a, b, Proportion, 0, UseVersion2, UseUP, IndexUP, Starting_n, false, MTest_Min, MTest_Max,
		MTest_a, MTest_b, BucketSize);
      
      cout << "Value of A (Variance of the modelisation noise) : " << a << endl;
      cout << "Value of B                                      : " << b << endl;
      cout << "Value of the variance ratio                     : " << fabs(a)/Variance << endl;
      
      ListA.push_back(a);
      ListB.push_back(b);
      ListVR.push_back(fabs(a)/Variance);
      ListR2.push_back(1 - fabs(a)/Variance);
      ListMask.push_back(SString.str());
      
      Min_a    = a;
      Min_b    = b;
      Min_VR   = fabs(a)/Variance;
      
      InputRemoved = CountMask;

      for(i=0; i<InputMask.size(); i++)
	{
	  for(j=0; j<InputMask.size(); j++)
	    {
	      if (InputMask[j].second!=MASK_X) continue;
	      
	      if ((InputMask.size()!=1)&&(InputRemoved!=1)) InputMask[j].first = 0;

	      SString.str("");
	      
	      cout << "Mask value : ";
	      for(k=0; k<InputMask.size(); k++)
		{
		  cout << InputMask[k].first;
		  SString << InputMask[k].first;
		} /* End For */
	      cout << endl;
	      
	      cout << "Mask value : ";
	      for(k=0; k<InputMask.size(); k++)
		{
		  switch(InputMask[k].second)
		    {
		    case MASK_0:
		      cout << 0;
		      break;
		    case MASK_1:
		      cout << 1;
		      break;
		    case MASK_X:
		      cout << "X";
		      break;
		    } /* End Switch */
		} /* End For */
	      cout << endl;
	      
	      // Si a ce moment, la retenue vaut 1 -> on met le
	      // drapeau de retenue a 1 de la variable correspondante
	      
	      cout << "Preparing the data with selected inputs" << endl;

	      Selected_InputData.resize(0);
	      
	      for(k=0; k<InputData.size(); k++)
		{
		  DList.resize(0);
		  for(l=0; l<InputMask.size(); l++)
		    {
		      if (InputMask[l].first==1)
			{
			  DList.push_back(InputData[k][l]);
			} /* End If */
		    } /* End For */
		  Selected_InputData.push_back(DList);
		} /* End For */

	      cout << "Performing the Gamma test" << endl;

	      ListDelta.resize(0);
	      ListGamma.resize(0);
	      ListAllDelta.resize(0);
	      ListAllGamma.resize(0);

	      GammaTest(Selected_InputData, MeasureData, ListDelta, ListGamma, ListAllDelta, ListAllGamma,
			p, a, b, Proportion, 0, UseVersion2, UseUP, IndexUP, Starting_n, false, MTest_Min, MTest_Max,
			MTest_a, MTest_b, BucketSize);
	      
	      cout << "Value of A (Variance of the modelisation noise) : " << a << endl;
	      cout << "Value of B                                      : " << b << endl;
	      cout << "Value of the variance ratio                     : " << fabs(a)/Variance << endl;

	      ListA.push_back(a);
	      ListB.push_back(b);
	      ListVR.push_back(fabs(a)/Variance);
	      ListR2.push_back(1 - fabs(a)/Variance);
	      ListMask.push_back(SString.str());

	      if (fabs(a)/Variance<Min_VR)
		{
		  Min_a    = a;
		  Min_b    = b;
		  Min_VR   = fabs(a)/Variance;
		  IndexMin = j;
		  OneInputSelected = true;
		} /* End If */

	      if (Last_Min_VR>Min_VR)
		{
		  Last_Min_VR = Min_VR;
		  Last_Min_a  = Min_a;
		  Last_Min_b  = Min_b;
		  Last_Min_Mask = SString.str();
		} /* End If */

	      // After having removed one  input we set all MASK_X inputs to 1

	      for(k=0; k<InputMask.size(); k++)
		{
		  if (InputMask[k].second==MASK_X) InputMask[k].first = 1;
		} /* End For */
	    } /* End For */
	  if (IndexMin!=-1)
	    {
	      InputMask[IndexMin].first  = 0;
	      InputMask[IndexMin].second = MASK_0;
	      InputRemoved--;
	    } /* End If */

	  IndexMin = -1;
	  Min_a  = numeric_limits<double>::max();
	  Min_b  = numeric_limits<double>::max();
	  Min_VR = numeric_limits<double>::max();
	} /* End For */
      
      // Computation of the sensitivity of the inputs

      List_Sensitivity.resize(0);
      
      Max = - numeric_limits<double>::max();

      vector<unsigned int> ListOfSelectedInputs(0);

      for(m=0; m<InputMask.size(); m++)
	{
	  if (Last_Min_Mask[m]=='1') ListOfSelectedInputs.push_back(m);
	} /* End For */

      if (ListOfSelectedInputs.size()>1)
	{
	  cout << "Computing the summed-up information related to the selected inputs" << endl;

	  for(m=0; m<ListOfSelectedInputs.size(); m++)
	    {
	      Selected_InputData.resize(0);
	      
	      for(k=0; k<InputData.size(); k++)
		{
		  DList.resize(0);
		  for(l=0; l<ListOfSelectedInputs.size(); l++)
		    {
		      if (l!=m)
			{
			  DList.push_back(InputData[k][ListOfSelectedInputs[l]]);
			} /* End If */
		    } /* End For */
		  Selected_InputData.push_back(DList);
		} /* End For */
	      
	      cout << "Performing the Gamma test without the selected input n° " << m << endl;
	      
	      ListDelta.resize(0);
	      ListGamma.resize(0);
	      ListAllDelta.resize(0);
	      ListAllGamma.resize(0);
	      
	      GammaTest(Selected_InputData, MeasureData, ListDelta, ListGamma, ListAllDelta, ListAllGamma,
			p, a, b, Proportion, 0, UseVersion2, UseUP, IndexUP, Starting_n, false, MTest_Min, MTest_Max,
			MTest_a, MTest_b, BucketSize);
	      
	      List_Sensitivity.push_back(1-fabs(a)/Variance);
	      
	      if (Max<(1-fabs(a)/Variance)) Max = 1-fabs(a)/Variance;
	    } /* End For */
	} /* End If */

      OutFile = new ofstream;
      
      OutFile->open("DataSelectInputsAB.dat");

      if (OneInputSelected)
	{
	  (*OutFile) << "Best selection :" << endl;
	  (*OutFile) << "Best Mask : " << Last_Min_Mask << endl;
	  (*OutFile) << "Y=B*X+A" << endl;
	  (*OutFile) << "B -- A -- VR -- R2" << endl;
	  (*OutFile) << Last_Min_b << " -- " << Last_Min_a << " -- " << Last_Min_VR << " -- " << 1 - Last_Min_VR << endl;
	} /* End If */
      else
	{
	  (*OutFile) << "Best selection :" << endl;
	  (*OutFile) << "No input selected !! " << endl;
	} /* End Else */

      (*OutFile) << "Mask value : ";
      for(j=0; j<InputMask.size(); j++)
	{
	  (*OutFile) << InputMask[j].first;
	} /* End For */
      (*OutFile) << endl;
      
      for(j=0; j<InputMask.size(); j++)
	{
	  switch(InputMask[j].second)
	    {
	    case MASK_0:
	      (*OutFile) << 0;
	      break;
	    case MASK_1:
	      (*OutFile) << 1;
	      break;
	    case MASK_X:
	      (*OutFile) << "X";
	      break;
	    } /* End Switch */
	} /* End For */
      (*OutFile) << endl;

      (*OutFile) << "% Y = A*X.B" << endl;
      (*OutFile) << "% Mask -- A -- B -- VR -- R2" << endl;

      for(i=0; i<ListA.size(); i++)
	{
	  (*OutFile) << ListMask[i] << " -- " << ListB[i] << " -- " << ListA[i] << " -- ";
	  (*OutFile) << ListVR[i]   << " -- " << ListR2[i] << endl;
	} /* End For */

      if (ListOfSelectedInputs.size()>1)
	{
	  if (List_Sensitivity.size()!=0)
	    {
	      cout << "Sensitivity of the output with respect to the inputs: " << endl;
	      
	      for(m=0; m<List_Sensitivity.size(); m++)
		{
		  cout << "Input " << m << " removed : decrease of R2 = " << (List_Sensitivity[m] - Max)/Max * 100.0 << endl;
		} /* End For */
	    } /* End If */
	} /* End If */

      OutFile->close();

      delete OutFile;
    } /* End If */

  if (S_Sensitivity)
    {
      ListA.resize(0);
      ListB.resize(0);
      ListVR.resize(0);
      ListR2.resize(0);

      cout << "Computation of the sensitivity of the inputs" << endl;

      // Computation of the Gamma Test values for the whole input data set

      cout << "Testing the whole set of inputs" << endl;

      cout << "Preparing the data with selected inputs" << endl;
      
      Selected_InputData.resize(0);
      
      for(k=0; k<InputData.size(); k++)
	{
	  Selected_InputData.push_back(InputData[k]);
	} /* End For */
      
      cout << "Performing the Gamma test" << endl;
      
      ListDelta.resize(0);
      ListGamma.resize(0);
      ListAllDelta.resize(0);
      ListAllGamma.resize(0);
      
      GammaTest(Selected_InputData, MeasureData, ListDelta, ListGamma, ListAllDelta, ListAllGamma,
		p, a, b, Proportion, 0, UseVersion2, UseUP, IndexUP, Starting_n, false, MTest_Min, MTest_Max,
		MTest_a, MTest_b, BucketSize);
      
      cout << "Value of A (Variance of the modelisation noise) : " << a << endl;
      cout << "Value of B                                      : " << b << endl;
      cout << "Value of the variance ratio                     : " << fabs(a)/Variance << endl;
      
      ListA.push_back(a);
      ListB.push_back(b);
      ListVR.push_back(fabs(a)/Variance);
      ListR2.push_back(1 - fabs(a)/Variance);
      
      Min_a    = a;
      Min_b    = b;
      Min_VR   = fabs(a)/Variance;

      // Computation of the sensitivity
      
      cout << "Computation of the sensitivity" << endl;
      
      for(i=0; i<InputData[0].size(); i++)
	{
	  cout << "Removing input " << i << " from the set of inputs" << endl;

	  cout << "Preparing the data with selected inputs" << endl;
	  
	  Selected_InputData.resize(0);
	  
	  for(k=0; k<InputData.size(); k++)
	    {
	      DList.resize(0);
	      for(l=0; l<InputData[k].size(); l++)
		{
		  if (l!=i) DList.push_back(InputData[k][l]);
		} /* End For */
	      Selected_InputData.push_back(DList);
	    } /* End For */
	  
	  cout << "Performing the Gamma test" << endl;
	  
	  ListDelta.resize(0);
	  ListGamma.resize(0);
	  ListAllDelta.resize(0);
	  ListAllGamma.resize(0);
	  
	  GammaTest(Selected_InputData, MeasureData, ListDelta, ListGamma, ListAllDelta, ListAllGamma,
		    p, a, b, Proportion, 0, UseVersion2, UseUP, IndexUP, Starting_n, false, MTest_Min, MTest_Max,
		    MTest_a, MTest_b, BucketSize);
	  
	  cout << "Value of A (Variance of the modelisation noise) : " << a << endl;
	  cout << "Value of B                                      : " << b << endl;
	  cout << "Value of the variance ratio                     : " << fabs(a)/Variance << endl;
	  
	  ListA.push_back(a);
	  ListB.push_back(b);
	  ListVR.push_back(fabs(a)/Variance);
	  ListR2.push_back(1 - fabs(a)/Variance);

	  if (fabs(a)/Variance<Min_VR)
	    {
	      Min_a    = a;
	      Min_b    = b;
	      Min_VR   = fabs(a)/Variance;
	    } /* End If */
	} /* End For */
      
      OutFile = new ofstream;
      
      OutFile->open("DataSensitivity.dat");

      (*OutFile) << "% Y = B*X+A" << endl;
      (*OutFile) << "% Input Removed -- B -- A -- VR -- R2" << endl;

      (*OutFile) << "No Input Removed" << " -- " << ListB[0] << " -- " << ListA[0] << " -- ";
      (*OutFile) << ListVR[0]   << " -- " << ListR2[0] << endl;

      for(i=1; i<ListA.size(); i++)
	{
	  (*OutFile) << i << " -- " << ListB[i] << " -- " << ListA[i] << " -- ";
	  (*OutFile) << ListVR[i]   << " -- " << ListR2[i] << endl;
	} /* End For */

      (*OutFile) << "Sensitivity of inputs : " << endl;

      for(i=1; i<ListA.size(); i++)
	{
	  if (!((ListR2[i]<0)||(ListR2[i]>1)))
	    {
	      (*OutFile) << i << " -- " << (ListR2[0] - ListR2[i])/(double)ListR2[0]*100.0 << " % " << endl;
	    } /* End If */
	  else
	    {
	      (*OutFile) << i << " -- 0 % " << endl;
	    } /* End Else */
	} /* End For */

      OutFile->close();

      delete OutFile;
    } /* End If */

  if (S_Estimator)
    {
      GammaTest_Estimator(InputData, MeasureData, Y_Output, p, 0, UseVersion2, Starting_n, BucketSize);      

      OutFile = new ofstream;
      
      OutFile->open("DataYEstim_GT.dat");

      for(i=0; i<Y_Output.size(); i++)
	{
	  (*OutFile) << Y_Output[i] << endl;
	} /* End If */

      OutFile->close();

      delete OutFile;
    } /* End If */

  return 0;
} /* End For n */
