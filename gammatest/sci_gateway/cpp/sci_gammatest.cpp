#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <iostream>
#include <math.h>

#ifdef _MSC_VER
#include <kdtree_static.h>
#endif

#include <GammaTest.h>
#include <tokenize.hpp>
#include <Filter.h>

extern "C" {
#include <stack-c.h>
#include <MALLOC.h>
#include <sciprint.h>
#include <Scierror.h>
}

#include <api_common.h>
#include <api_double.h>
#include <helper.h>

using namespace std;

#define DATA_IN    1
#define N_IN       2
#define P_IN       3
#define PARAM_IN   4
#define GT_A_OUT   Rhs+1
#define GT_B_OUT   Rhs+2
#define REG_D_OUT  Rhs+3
#define REG_G_OUT  Rhs+4
#define SCAT_D_OUT Rhs+5
#define SCAT_G_OUT Rhs+6

// #define LOG 1

extern "C" int sci_gammatest(char * fname)
{
  vector<vector<double> >  InputData, Selected_InputData;
  vector<double>           MeasureData, Selected_MeasureData;
  vector<double>           ListDelta;
  vector<double>           ListGamma;
  vector<double>           ListAllDelta;
  vector<double>           ListAllGamma, MTest_a, MTest_b;
  vector<double>           InputMin, InputMax;
  double                   MeasureMin, MeasureMax;
  double                   a = 0.0, b = 0.0;
  unsigned int             i, j;
  bool                     S_Select = false, S_Rand = false, S_Greedy = false, S_Greedy_Fast = false, UseUP = false;
  bool                     S_Input_Select_Exhaust = false, UseVersion2 = false;
  bool                     S_Min = true, S_MTest = false, S_Norm = false, S_Estimator = false, S_BucketSize = false;
  int                      BucketSize = 5, S_Nb_Select = 0, S_Nb_Rand = 0, S_Nb_Greedy = 0, S_Nb_Greedy_Fast = 0;
  unsigned int             S_Start_Point = 0, S_End_Point = 0;
  unsigned int             Retenue = 0, Index = 0, IndexUP = 0, RandSeed = 0;
  int                      IndexMin, Starting_n = 0, MTest_Min = 0, MTest_Max = 0;
  double                   Variance = 0.0, Mean = 0.0, Max = - numeric_limits<double>::max();
  SciErr _SciErr;

  // [gt_a, gt_b, reg_d, reg_g, scat_d, scat_g] = gammatest(data_in,n,p,params);

  // 'select', 1
  // 'greedy', 1
  // 'greedy_fast', 1
  // 'rand', 1
  // 'starting_n', 1
  // 'mtest_min', -1
  // 'mtest_max', -1
  // 'randseed', 0
  // 'norm', 0
  // 'bucket_size', 5
  // 'start_point', 0
  // 'end_point', 0

  int * data_in_addr = NULL, * n_addr = NULL, * p_addr = NULL;
  int * scat_d_out = NULL, * scat_g_out = NULL;
  int * reg_d_out = NULL, * reg_g_out = NULL;
  int * gt_a = NULL, * gt_b = NULL;
  int scat_d_rows, scat_d_cols;
  int scat_g_rows, scat_g_cols;
  int reg_d_rows,  reg_d_cols;
  int reg_g_rows,  reg_g_cols;
  int gt_a_rows,   gt_a_cols;
  int gt_b_rows,   gt_b_cols;
  int data_in_rows, data_in_cols;
  double * data_in = NULL;
  double n, p;

  _SciErr = getVarAddressFromPosition(pvApiCtx, DATA_IN, &data_in_addr);
  _SciErr = getMatrixOfDouble(pvApiCtx, data_in_addr, &data_in_rows, &data_in_cols, &data_in);

  _SciErr = getVarAddressFromPosition(pvApiCtx, N_IN, &n_addr);
  if (!isEmptyMatrix(pvApiCtx, n_addr))
    {
      getScalarDouble(pvApiCtx, n_addr, &n);
    }
  else
    {
      n = 20;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, P_IN, &p_addr);
  if (!isEmptyMatrix(pvApiCtx, p_addr))
    {
      getScalarDouble(pvApiCtx, p_addr, &p);
    }
  else
    {
      p = 1.0;
    }
	
  int tmp_int, tmp_res;
  double tmp_dbl;
  char * tmp_char = NULL;

  INIT_PARAM(PARAM_IN);

  GET_PARAM_INT("select", tmp_int, 1, tmp_res);
  if (tmp_res!=-1) 
    {
      S_Nb_Select = tmp_int;
      S_Select = true;
    }
  else
    {
      S_Select = false;
    }

  GET_PARAM_INT("greedy", tmp_int, 1, tmp_res);
  if (tmp_res!=-1) 
    {
      S_Nb_Greedy = tmp_int;
      S_Greedy = true;
    }
  else
    {
      S_Greedy = false;
    }
  
  GET_PARAM_INT("greedy_fast", tmp_int, 1, tmp_res);
  if (tmp_res!=-1) 
    {
      S_Nb_Greedy_Fast = tmp_int;
      S_Greedy_Fast = true;
    }
  else
    {
      S_Greedy_Fast = false;
    }

  GET_PARAM_INT("rand", tmp_int, 1, tmp_res);
  if (tmp_res!=-1) 
    {
      S_Nb_Rand = tmp_int;
      S_Rand = true;
    }
  else
    {
      S_Rand = false;
    }

  GET_PARAM_INT("useup", tmp_int, 0, tmp_res);
  if (tmp_res!=-1) 
    {
      IndexUP = tmp_int;
      UseUP   = true;
    }
  else
    {
      UseUP = false;
    }

  GET_PARAM_INT("starting_n", tmp_int, 1, tmp_res);
  Starting_n = tmp_int;

  GET_PARAM_INT("mtest_min", tmp_int, -1, tmp_res);
  if (tmp_res!=-1) 
    {
      MTest_Min = tmp_int;
      S_MTest   = true;
    }
  else
    {
      S_MTest = false;
    }

  GET_PARAM_INT("mtest_max", tmp_int, -1, tmp_res);
  if (tmp_res!=-1) 
    {
      MTest_Max = tmp_int;
      S_MTest   = true;
    }
  else
    {
      S_MTest = S_MTest | false;
    }

  if ((MTest_Max==-1)&&(MTest_Min==-1)) 
    {
      MTest_Max = p;
      MTest_Min = p;
    }
  else if ((MTest_Max==-1)&&(MTest_Min!=-1))
    {
      MTest_Max = MTest_Min;
    }
  else if ((MTest_Max!=-1)&&(MTest_Min==-1))
    {
      MTest_Min = MTest_Max;
    }

  GET_PARAM_INT("norm", tmp_int, 1, tmp_res);
  S_Norm = (tmp_res!=-1)&&(tmp_int);

  GET_PARAM_INT("seed", tmp_int, 0, tmp_res);
  RandSeed = tmp_int;

  GET_PARAM_INT("estimator", tmp_int, 0, tmp_res);
  S_Estimator = (tmp_res!=-1)&&(tmp_int);

  GET_PARAM_INT("userversion2", tmp_int, 0, tmp_res);
  UseVersion2 = (tmp_res!=-1)&&(tmp_int);

  GET_PARAM_INT("bucket_size", tmp_int, 5, tmp_res);
  BucketSize   = tmp_int;
  S_BucketSize = (tmp_res!=-1);

  GET_PARAM_INT("start_point", tmp_int, 0, tmp_res);
  S_Start_Point = tmp_int;

  GET_PARAM_INT("end_point", tmp_int, 0, tmp_res);
  S_End_Point = tmp_int;

  // Start Gamma test pre-processing
  
  srand(RandSeed);

  ListDelta.resize(0);
  ListGamma.resize(0);
  ListAllDelta.resize(0);
  ListAllGamma.resize(0);

  InputMin.resize(0);
  InputMax.resize(0);

  MTest_a.resize(0);
  MTest_b.resize(0);

  // We store data_in into InputData

  InputData.resize(data_in_rows, vector<double>(data_in_cols, 0.0));
  InputMin.resize(data_in_cols-1, numeric_limits<double>::max());
  InputMax.resize(data_in_cols-1, numeric_limits<double>::min());
  MeasureData.resize(data_in_rows, 0.0);
  MeasureMin = numeric_limits<double>::max();
  MeasureMax = numeric_limits<double>::min();

  for(i=0;i<data_in_rows;i++)
    {
      for(j=0;j<data_in_cols-1; j++)
	{
	  InputData[i][j] = *(data_in + i + j*data_in_rows);
	  if (InputMin[j]>InputData[i][j]) InputMin[j] = InputData[i][j];
	  if (InputMax[j]<InputData[i][j]) InputMax[j] = InputData[i][j];
	}
      MeasureData[i] = *(data_in + i + (data_in_cols-1)*data_in_rows);
      if (MeasureMin>MeasureData[i]) MeasureMin = MeasureData[i];
      if (MeasureMax<MeasureData[i]) MeasureMax = MeasureData[i];
    }

  // Retrieving the S_Start_Point - S_End_Point interval
  
  if ((S_Start_Point==0)&&(S_End_Point==0))
    {
      S_Start_Point = 0;
      S_End_Point   = InputData.size();
    }
  else if (S_End_Point==0)
    {
      S_End_Point = InputData.size();
    }
  else if (S_End_Point<S_Start_Point)
    {
      Scierror(999,"%s: Error - S_Start_Point must be smaller than S_End_Point\n", fname);
      return 0;
    }

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
      for(i=0; i<InputData.size(); i++)
	{
	  for(j=0; j<InputData[i].size(); j++)
	    {
	      InputData[i][j] = (InputData[i][j] - InputMin[j])/(double)(InputMax[j] - InputMin[j]);
	    }
	  MeasureData[i] = (MeasureData[i] - MeasureMin)/(double)(MeasureMax - MeasureMin);
	}

      for(i=0; i<InputMin.size(); i++)
	{
	  InputMin[i] = 0;
	  InputMax[i] = 1;
	}

      MeasureMin = 0;
      MeasureMax = 1;
    }

  if (S_Select)      Filter_Select(InputData, MeasureData, S_Nb_Select);
  if (S_Greedy)      Filter_EchantGreedy(InputData, MeasureData, InputMin, InputMax, S_Nb_Greedy);
  if (S_Greedy_Fast) Filter_EchantGreedy_Fast(InputData, MeasureData, InputMin, InputMax,
					      (InputData.size() - S_Nb_Greedy_Fast)/100, S_Nb_Greedy_Fast);
  if (S_Rand)        Filter_EchantRand(InputData, MeasureData, S_Nb_Rand);

  // Computation of the Variance Ratio

  for(i=0; i<MeasureData.size(); i++) Mean += MeasureData[i];
  Mean /= (double)MeasureData.size();

  for(i=0; i<MeasureData.size(); i++) Variance += (MeasureData[i] - Mean) * (MeasureData[i] - Mean);
  Variance /= (double)MeasureData.size();

  // Launch the Gamma Test
  GammaTest(InputData, MeasureData, ListDelta, ListGamma, ListAllDelta, ListAllGamma,
	    (int)n, a, b, p, UseVersion2, UseUP, IndexUP, Starting_n, S_MTest, MTest_Min, MTest_Max,
	    MTest_a, MTest_b, BucketSize);
  
  if (S_MTest)
    {
      gt_a_rows = MTest_a.size(); gt_a_cols = 1;
      gt_b_rows = MTest_a.size(); gt_b_cols = 1;
      double * gt_a = (double *)MALLOC(gt_a_rows*sizeof(double));
      double * gt_b = (double *)MALLOC(gt_b_rows*sizeof(double));

      for(i=0; i<MTest_a.size(); i++)
	{
	  gt_a[i] = MTest_a[i];
	  gt_b[i] = MTest_b[i];
	} /* End For */
      _SciErr = createMatrixOfDouble(pvApiCtx, GT_A_OUT, gt_a_rows, gt_a_cols, gt_a);
      _SciErr = createMatrixOfDouble(pvApiCtx, GT_B_OUT, gt_b_rows, gt_b_cols, gt_b);

      FREE(gt_a);
      FREE(gt_b);
    } /* End If */
  else
    {
      gt_a_rows = 1; gt_a_cols = 1;
      gt_b_rows = 1; gt_b_cols = 1;
      double * gt_a = (double *)MALLOC(gt_a_rows*sizeof(double));
      double * gt_b = (double *)MALLOC(gt_b_rows*sizeof(double));

      gt_a[0] = a;
      gt_b[0] = b;

      _SciErr = createMatrixOfDouble(pvApiCtx, GT_A_OUT, gt_a_rows, gt_a_cols, gt_a);
      _SciErr = createMatrixOfDouble(pvApiCtx, GT_B_OUT, gt_b_rows, gt_b_cols, gt_b);

      FREE(gt_a);
      FREE(gt_b);
    }

  LhsVar(1) = GT_A_OUT;
  LhsVar(2) = GT_B_OUT;

  if (Lhs>=3)
    {
      // Regression points
      reg_d_rows = ListDelta.size(); reg_d_cols = 1;
      reg_g_rows = ListDelta.size(); reg_g_cols = 1;
      double * reg_d = (double *)MALLOC(reg_d_rows*sizeof(double));
      double * reg_g = (double *)MALLOC(reg_g_rows*sizeof(double));
      
      for(i=0; i<ListDelta.size(); i++)
	{
	  reg_d[i] = ListDelta[i];
	  reg_g[i] = ListGamma[i];
	} /* End For */
      
      _SciErr = createMatrixOfDouble(pvApiCtx, REG_D_OUT, reg_d_rows, reg_d_cols, reg_d);
      _SciErr = createMatrixOfDouble(pvApiCtx, REG_G_OUT, reg_g_rows, reg_g_cols, reg_g);
      
      FREE(reg_d);
      FREE(reg_g);

      LhsVar(3) = REG_D_OUT;
      LhsVar(4) = REG_G_OUT;
    }

  if (Lhs>=5)
    {
      // All the points for the scatterplot
      scat_d_rows = ListAllDelta.size(); scat_d_cols = 1;
      scat_g_rows = ListAllDelta.size(); scat_g_cols = 1;
      double * scat_d = (double *)MALLOC(scat_d_rows*sizeof(double));
      double * scat_g = (double *)MALLOC(scat_g_rows*sizeof(double));
      
      for(i=0; i<ListAllDelta.size(); i++)
	{
	  scat_d[i] = ListAllDelta[i];
	  scat_g[i] = ListAllGamma[i];
	} /* End For */
      
      _SciErr = createMatrixOfDouble(pvApiCtx, SCAT_D_OUT, scat_d_rows, scat_d_cols, scat_d);
      _SciErr = createMatrixOfDouble(pvApiCtx, SCAT_G_OUT, scat_g_rows, scat_g_cols, scat_g);
      
      FREE(scat_d);
      FREE(scat_g);

      LhsVar(5) = SCAT_D_OUT;
      LhsVar(6) = SCAT_G_OUT;
    }

  // (*OutFile) << "Gamma = G * Delta + A" << endl;
  // (*OutFile) << "G = " << b << " A = " << a << endl;
  
  return 0;
}
