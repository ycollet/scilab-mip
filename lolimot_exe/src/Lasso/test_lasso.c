#include <stdlib.h>
#include <math.h>

#include <lasso.h>

#define NUMOFPARAM       10
#define NUMOFPOINTS      100
#define NUMOFVALIDPOINTS 10
#define NBPROBLEM        30

//#define DONTDISPLAYDATA
#define RANDSEED 1234
#define MAX(A,B) ((A<B)?B:A)
#define FILENAME "data_"
#define RANDINITPARAMETERS

double ComputeResidual_Validation(double * XValid, double * YValid, double * ModelParam, int IndexPb)
{
  unsigned int i, j;
  double       Residual = 0.0, YHatValid;

  for(i=0; i<NUMOFVALIDPOINTS; i++)
    {
      YHatValid = 0.0;
      for(j=0; j<NUMOFPARAM; j++)
	{
	  YHatValid += ModelParam[IndexPb*NUMOFPARAM+j]*XValid[i+j*NUMOFVALIDPOINTS];
	} /* End For */
      Residual += (YHatValid - YValid[i]) * (YHatValid - YValid[i]);
    } /* End For */

  return Residual;
}

int main()
{
  double * X = NULL, * XValid = NULL;
  double * Y = NULL, * YHat = NULL, * YValid = NULL, * YHatValid = NULL;
  double * Residual = NULL, * Bound = NULL, * Lagrangian = NULL, * ValidResidual = NULL;
  double * Param = NULL, Sum = 0.0, Sum_L2 = 0.0, Sum_L1 = 0.0;
  double * SelectionParam = NULL, * ListOfResValid = NULL;
  int      Verbose = 1, PasSub = 0, PSuc, i, j, NbProblem = NBPROBLEM;
  int      NumOfParam = NUMOFPARAM, NumOfPoints = NUMOFPOINTS;
  double   Step, ModelParam[NUMOFPARAM], ResValid;
  char   * Filename;
  FILE   * fid;

  printf("*********\n");
  printf("* lasso *\n");
  printf("*********\n");

  X             = (double *)malloc(NUMOFPARAM*NUMOFPOINTS*sizeof(double));
  Y             = (double *)malloc(NUMOFPOINTS*sizeof(double));
  YHat          = (double *)malloc(NUMOFPOINTS*sizeof(double));
  XValid        = (double *)malloc(NUMOFPARAM*NUMOFVALIDPOINTS*sizeof(double));
  YValid        = (double *)malloc(NUMOFVALIDPOINTS*sizeof(double));
  YHatValid     = (double *)malloc(NUMOFVALIDPOINTS*sizeof(double));
  Residual      = (double *)malloc(NUMOFPOINTS*sizeof(double));
  Param         = (double *)malloc(NUMOFPARAM*sizeof(double));
  Bound         = (double *)malloc(1*sizeof(double));
  Lagrangian    = (double *)malloc(1*sizeof(double));
  ValidResidual = (double *)malloc(1*sizeof(double));

  // Random selection of parameters
#ifdef RANDINITPARAMETERS
  srand(RANDSEED);
  printf("List of the parameters of the model\n");
  ModelParam[0] = rand()/(double)RAND_MAX;
  printf("ModelParam[0] = %f\n", ModelParam[0]);
  for(i=1; i<NUMOFPARAM; i++) 
    {
      ModelParam[i] = rand()/(double)RAND_MAX;
      printf("ModelParam[%d] = %f\n", i, ModelParam[i]);
    } /* End For */
#else
  for(i=0; i<NUMOFPARAM; i++) Param[i] = 0.0;
#endif

  // Generation of the learning data set
  for(i=0; i<NUMOFPOINTS; i++)
    {
      X[i+0*NUMOFPOINTS] = 1.0;
      for(j=1; j<NUMOFPARAM; j++) X[i+j*NUMOFPOINTS] = rand()/(double)RAND_MAX;
      Y[i] = ModelParam[0] * 1.0;
      for(j=1; j<NUMOFPARAM; j++) Y[i] += ModelParam[j]*X[i+j*NUMOFPOINTS]*(1+20.0*X[i+j*NUMOFPOINTS]);
    } /* End For */

  // Generation of the validation data set
  for(i=0; i<NUMOFVALIDPOINTS; i++)
    {
      XValid[i+0*NUMOFVALIDPOINTS] = 1.0;
      for(j=1; j<NUMOFPARAM; j++) 
	XValid[i+j*NUMOFVALIDPOINTS] = rand()/(double)RAND_MAX;
      YValid[i] = ModelParam[0] * 1.0;
      for(j=1; j<NUMOFPARAM; j++) 
	YValid[i] += ModelParam[j]*XValid[i+j*NUMOFVALIDPOINTS]*(1+20.0*XValid[i+j*NUMOFVALIDPOINTS]);
    } /* End For */

  Bound[0] = 0.0;
  for(i=0; i<NUMOFPOINTS; i++) 
    {
      Sum = 0.0;
      for(j=0; j<NUMOFPARAM; j++)
	{
	  Sum += (Y[i]/MAX(X[i+j*NUMOFPOINTS],0.001));
	} /* End For */
      Sum /=(double)NUMOFPARAM;
      Bound[0] += Sum;
    } /* End For */
  Bound[0] /= (double)NUMOFPOINTS;

  Lagrangian[0] = 1.0;

  lasso(X, &NumOfPoints, &NumOfParam, Bound, Param, Y, YHat, Residual, Lagrangian, &PSuc,  &Verbose, &PasSub);

#ifndef DONTDISPLAYDATA
  for(i=0; i<NUMOFPOINTS; i++) printf("Y = %f / YHat = %f / Residual = %f\n", Y[i], YHat[i], Residual[i]);
#endif

  for(i=0; i<NUMOFPARAM; i++)  printf("Param[%d] = %f\n", i, Param[i]);

  Filename = (char *)malloc(100*sizeof(char));

  sprintf(Filename, "%s%s.dat",FILENAME, "lasso");
  printf("Saving file %s \n");
  fid = fopen(Filename, "w");
  for(i=0; i<NUMOFPOINTS; i++) fprintf(fid, "%f %f\n", Y[i], YHat[i]);
  fclose(fid);

  Sum = 0.0;
  for(i=0; i<NUMOFPOINTS; i++) Sum += Residual[i];
  Sum /= (double)NUMOFPOINTS;

  Sum_L2 = 0.0;
  for(i=0; i<NUMOFPOINTS; i++) Sum_L2 += Residual[i]*Residual[i];
  Sum_L2 /= (double)NUMOFPOINTS;

  ResValid = ComputeResidual_Validation(XValid, YValid, Param, 0);

  printf("PSuc       = %d\n", PSuc);
  printf("Residual : Mean = %f - Std = %f\n", Sum, sqrt(Sum_L2));
  printf("Lagrangian = %f\n", Lagrangian);
  printf("PasSub     = %d\n", PasSub);
  printf("Verbose    = %d\n", Verbose);
  printf("Validation Residual = %f\n", ResValid);

  // We compute the L1 norm of the parameters. The result is then used as a bound for mult_lasso
  Sum_L1 = 0.0;
  for(i=0;i<NUMOFPARAM; i++) Sum_L1 += fabs(Param[i]);

  free(YHat);
  free(Residual);
  free(Bound);
  free(Param);
  free(Lagrangian);
  free(ValidResidual);

  printf("**************\n");
  printf("* mult_lasso *\n");
  printf("**************\n");

  YHat           = (double *)malloc(NUMOFPOINTS*NBPROBLEM*sizeof(double));
  Residual       = (double *)malloc(NUMOFPOINTS*NBPROBLEM*sizeof(double));
  Param          = (double *)malloc(NUMOFPARAM*NBPROBLEM*sizeof(double));
  Bound          = (double *)malloc(NBPROBLEM*sizeof(double));
  Lagrangian     = (double *)malloc(NBPROBLEM*sizeof(double));
  SelectionParam = (double *)malloc(NBPROBLEM*sizeof(double));
  ValidResidual  = (double *)malloc(NBPROBLEM*sizeof(double));
  ListOfResValid = (double *)malloc(NBPROBLEM*sizeof(double));

  for(i=0; i<NBPROBLEM*NUMOFPARAM; i++) Param[i] = 0.0;

  printf("L1 norm of the parameters = %f\n", Sum_L1);

  // We divide by NbProblem so to have the last Bound != 0
  for(i=0, Step = Sum_L1; i<NbProblem; i++, Step -= Sum_L1 / (double)(NbProblem)) 
    {
      Bound[i] = Step;
      printf("Bound[%d] = %f\n", i, Step);
      Lagrangian[i] = 1.0;
    } /* End For */

  mult_lasso(X, &NumOfPoints, &NumOfParam, Bound, &NbProblem, Param, Y, YHat, Residual, Lagrangian, &PSuc,  &Verbose);

  for(i=0; i<NBPROBLEM; i++)
    {
      printf("Problem %d:\n", i);

#ifndef DONTDISPLAYDATA
      for(j=0; j<NUMOFPOINTS; j++)
	{
	  printf("Y = %f / YHat = %f / Residual = %f\n", Y[j], YHat[i*NUMOFPOINTS+j], Residual[i*NUMOFPOINTS+j]);
	} /* End For */
#endif

      ListOfResValid[i] = ComputeResidual_Validation(XValid, YValid, Param, i);

      sprintf(Filename, "%s%s%d.dat",FILENAME, "mult_lasso_", i);
      fid = fopen(Filename, "w");
      for(j=0; j<NUMOFPOINTS; j++) fprintf(fid, "%f %f\n", Y[j], YHat[i*NUMOFPOINTS+j]);
      fclose(fid);

      for(j=0; j<NUMOFPARAM; j++)
	{
	  printf("Param[%d] = %f\n", j, Param[i*NUMOFPARAM+j]);
	} /* End For */
      printf("Lagrangian[%d] = %f\n", i, Lagrangian[i]);

      Sum = 0.0;
      for(j=0; j<NUMOFPOINTS; j++) Sum += Residual[i*NUMOFPOINTS+j];
      Sum /= (double)NUMOFPOINTS;
      
      Sum_L2 = 0.0;
      for(j=0; j<NUMOFPOINTS; j++) Sum_L2 += Residual[i*NUMOFPOINTS+j]*Residual[i*NUMOFPOINTS+j];
      Sum_L2 /= (double)NUMOFPOINTS;

      printf("Residual[%d] : Mean = %f - Std = %f\n", i, Sum, sqrt(Sum_L2));
      printf("Lagrangian[%d] = %f - L1 = %f\n", i, Lagrangian[i], Sum_L1);
      //SelectionParam[i] = Sum + Lagrangian[i]*Sum_L1;
      SelectionParam[i] = Sum*(double)NUMOFPOINTS + 2.0*Sum_L1; // AIC = Resi + 2*K
    } /* End For */

  printf("PSuc = %d\n", PSuc);
  printf("Residual = %f\n", Sum);
  printf("PasSub = %d\n", PasSub);
  printf("Verbose = %d\n", Verbose);

  sprintf(Filename, "%s%s%s.dat",FILENAME, "mult_lasso_", "selparam");
  printf("Saving file %s\n", Filename);
  fid = fopen(Filename, "w");
  for(j=0; j<NBPROBLEM; j++) fprintf(fid, "%f %f %f\n", Bound[j], SelectionParam[j], ListOfResValid[j]);
  fclose(fid);

  sprintf(Filename, "%s%s%s.dat",FILENAME, "mult_lasso_", "listofparam");
  printf("Saving file %s\n", Filename);
  fid = fopen(Filename, "w");
  for(i=0; i<NBPROBLEM; i++)
    {
      for(j=0; j<NUMOFPARAM; j++)
	{
	  fprintf(fid, "%f ", Param[i*NUMOFPARAM+j]);
	} /* End For */
      fprintf(fid, "\n");
    } /* End For */
  fclose(fid);

  free(X);
  free(Y);
  free(YHat);
  free(Residual);
  free(Param);
  free(Bound);
  free(Lagrangian);
  free(ValidResidual);

  free(Filename);

  return 0;
}
