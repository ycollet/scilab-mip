/*
 * Scilab interface for the FANN library
 * Author: Dirk Gorissen <dirk.gorissen@ua.ac.be>
 * Author: Yann COLLETTE <ycollet@freesurf.fr>
 * Licence: GPL version 2 or later
 */

#include <helperFann.h>

#include <stdio.h>
#include <math.h>
#include <string.h>

#include <stack-c.h>
#include <api_scilab.h>
#include <MALLOC.h>
#include <Scierror.h>

//--------------------------------------------------------------------------------------------------------
// Function code adapted from:
// http://leenissen.dk/fann/forum/viewtopic.php?p=719&sid=1661ac359e28908e704231faa6310518 
//
struct fann_train_data *read_from_array(const double *din,
					const double *dout,
					const unsigned int num_data,
					const unsigned int num_input,
					const unsigned int num_output) 
{  
  unsigned int i, j;
  struct fann_train_data * data = (struct fann_train_data *)MALLOC(sizeof(struct fann_train_data));
  if (data == NULL) 
    {
      fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
      return NULL;
    }
 
  fann_init_error_data((struct fann_error *) data);
 
  data->num_data   = num_data;
  data->num_input  = num_input;
  data->num_output = num_output;

  data->input = (fann_type **)MALLOC(num_data * sizeof(fann_type *));
  if (data->input == NULL)
    {
      fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
      fann_destroy_train(data);
      return NULL;
    }
 
  data->output = (fann_type **)MALLOC(num_data * sizeof(fann_type *));
  if (data->output == NULL) 
    {
      fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
      fann_destroy_train(data);
      return NULL;
    }

  for(i=0; i<num_data; i++)
    {
      data->input[i] = (fann_type *)MALLOC(num_input*sizeof(fann_type));
      if (data->input[i] == NULL) 
	{
	  fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
	  fann_destroy_train(data);
	  return NULL;
	}

      data->output[i] = (fann_type *)MALLOC(num_output*sizeof(fann_type));
      if (data->output[i] == NULL) 
	{
	  fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
	  fann_destroy_train(data);
	  return NULL;
	}
    }
 
  //Code changed to support the way scilab passes arrays
  for(i=0; i<num_data; i++) 
    {
      for(j=0; j<num_input; j++) 
	{
	  data->input[i][j] = din[(j*num_data)+i];
	}
      
      for(j=0; j!=num_output; j++) 
	{
	  data->output[i][j] = dout[(j*num_data)+i];
	}
    }

  return data;
}

//--------------------------------------------------------------------------------------------------------
//Evaluate the ann on an array of samples
void evaluateNetwork(struct fann *ann, 
		     const double *input, 
		     double* output,
		     const unsigned int numData)
{
  int i,j;
  unsigned int numInputs  = ann->num_input;
  unsigned int numOutputs = ann->num_output;
  fann_type * out = NULL;
  fann_type * in = (fann_type *)MALLOC(numInputs * sizeof(fann_type));
  if (in == NULL) 
    {
      fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
    }
  
  for(i=0;i<(int)numData;++i)
    {
      for(j=0;j<(int)numInputs;j++) 
	{
	  in[j] = input[(j*numData)+i];
	}
      
      out = fann_run(ann,in);
      
      for(j=0;j<(int)numOutputs;j++) 
	{
	  output[(j*numData)+i] = out[j];
	}
    }
  
  FREE(in);
}

//--------------------------------------------------------------------------------------------------------
int detect_fanntraindatalist(int StackPos)
{
  int m_param, n_param, * pi_param_addr;
  int type, * pi_len = NULL, i;
  char ** LabelList = NULL;
  SciErr _sciErr;

  _sciErr = getVarAddressFromPosition(pvApiCtx, StackPos, &pi_param_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }
  _sciErr = getVarType(pvApiCtx, pi_param_addr, &type);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }

  if (type!=sci_mlist)
    {
      Scierror(999,"Argument %d is not a mlist\n",StackPos);
      return 0;
    }

  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, NULL, NULL);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }

  pi_len = (int*)MALLOC(sizeof(int) * m_param * n_param);
  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, pi_len, NULL);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }
  
  LabelList = (char**)MALLOC(sizeof(char*) * m_param * n_param);
  for(i=0; i<m_param*n_param; i++)
    {
      LabelList[i] = (char*)MALLOC(sizeof(char) * (pi_len[i] + 1)); // + 1 for null termination
    }
  
  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, pi_len, LabelList);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }
  
  if (strcmp(LabelList[0],"fanntraindatalist") != 0) 
    {
      Scierror(999,"Argument %d is not a fanntraindatalist\n",StackPos);
      return 0;
    }
  
  return 0;
}

//--------------------------------------------------------------------------------------------------------
int detect_fannlist(int StackPos)
{
  int m_param, n_param, * pi_param_addr;
  int type, * pi_len = NULL, i;
  char ** LabelList = NULL;
  SciErr _sciErr;

  _sciErr = getVarAddressFromPosition(pvApiCtx, StackPos, &pi_param_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }
  _sciErr = getVarType(pvApiCtx, pi_param_addr, &type);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }

  if (type!=sci_mlist)
    {
      Scierror(999,"Argument %d is not a mlist\n",StackPos);
      return 0;
    }

  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, NULL, NULL);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }

  pi_len = (int*)MALLOC(sizeof(int) * m_param * n_param);
  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, pi_len, NULL);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }
  
  LabelList = (char**)MALLOC(sizeof(char*) * m_param * n_param);
  for(i=0; i<m_param*n_param; i++)
    {
      LabelList[i] = (char*)MALLOC(sizeof(char) * (pi_len[i] + 1)); // + 1 for null termination
    }
  
  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, pi_len, LabelList);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }

  if (strcmp(LabelList[0],"fannlist") != 0) 
    {
      Scierror(999,"Argument %d is not a fannlist\n",StackPos);
      return 0;
    }
  
  return 0;
}

//--------------------------------------------------------------------------------------------------------
int detect_fannerrorlist(int StackPos)
{
  int m_param, n_param, * pi_param_addr;
  int type, * pi_len = NULL, i;
  char ** LabelList = NULL;
  SciErr _sciErr;

  _sciErr = getVarAddressFromPosition(pvApiCtx, StackPos, &pi_param_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }
  _sciErr = getVarType(pvApiCtx, pi_param_addr, &type);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }

  if (type!=sci_mlist)
    {
      Scierror(999,"Argument %d is not a mlist\n",StackPos);
      return 0;
    }

  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, NULL, NULL);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }

  pi_len = (int*)MALLOC(sizeof(int) * m_param * n_param);
  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, pi_len, NULL);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }
  
  LabelList = (char**)MALLOC(sizeof(char*) * m_param * n_param);
  for(i=0; i<m_param*n_param; i++)
    {
      LabelList[i] = (char*)MALLOC(sizeof(char) * (pi_len[i] + 1)); // + 1 for null termination
    }
  
  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, pi_len, LabelList);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }

  if (strcmp(LabelList[0],"fannerrorlist") != 0) 
    {
      Scierror(999,"Argument %d is not a fannerrorlist\n",StackPos);
      return 0;
    }
  
  return 0;
}

//--------------------------------------------------------------------------------------------------------
int createScilabFannStructFromCFannStruct(struct fann* ann,unsigned int StackPos)
{
  int * pi_extra_addr = NULL;
  SciErr _sciErr;

  //The struct field names
  static char * fieldnames[] = {"fannlist", "fann"};

  _sciErr = createMList(pvApiCtx, StackPos, 2, &pi_extra_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }

  _sciErr = createMatrixOfStringInList(pvApiCtx, StackPos, pi_extra_addr, 1, 2, 1, fieldnames);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }

  _sciErr = createPointerInList(pvApiCtx, StackPos, pi_extra_addr, 2, ann);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }

  return 0;
}

//--------------------------------------------------------------------------------------------------------
int createScilabFannTrainDataStructFromCFannTrainDataStruct(struct fann_train_data* ann_data, unsigned int StackPos)
{
  int * pi_extra_addr = NULL;
  SciErr _sciErr;

  //The struct field names
  static char * fieldnames[] = {"fanntraindatalist", "fann_train_data"};

  _sciErr = createMList(pvApiCtx, StackPos, 2, &pi_extra_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, -1);
      return 0;
    }

  _sciErr = createMatrixOfStringInList(pvApiCtx, StackPos, pi_extra_addr, 1, 2, 1, fieldnames);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, -1);
      return 0;
    }

  _sciErr = createPointerInList(pvApiCtx, StackPos, pi_extra_addr, 2, ann_data);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, -1);
      return 0;
    }

  return 0;
}

//--------------------------------------------------------------------------------------------------------
int createScilabFannErrorStructFromCFannErrorStruct(struct fann_error* ann, FILE * log_file, unsigned int StackPos)
{
  int * pi_extra_addr = NULL;
  SciErr _sciErr;

  //The struct field names
  static char * fieldnames[] = {"fannerrorlist", "file"};

  _sciErr = createMList(pvApiCtx, StackPos, 2, &pi_extra_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }

  _sciErr = createMatrixOfStringInList(pvApiCtx, StackPos, pi_extra_addr, 1, 2, 1, fieldnames);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }

  _sciErr = createPointerInList(pvApiCtx, StackPos, pi_extra_addr, 2, log_file);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return -1;
    }

  return 0;
}

//--------------------------------------------------------------------------------------------------------
struct fann * createCFannStructFromScilabFannStruct(unsigned int StackPos, int * res)
{
  int m_param, n_param, * pi_param_addr;
  int type, * pi_len = NULL, i;
  char ** LabelList = NULL;
  struct fann * result_ann = NULL;
  SciErr _sciErr;

  *res = -1;

  _sciErr = getVarAddressFromPosition(pvApiCtx, StackPos, &pi_param_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }
  _sciErr = getVarType(pvApiCtx, pi_param_addr, &type);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }

  if (type!=sci_mlist)
    {
      Scierror(999,"Argument %d is not a mlist\n",StackPos);
      return NULL;
    }

  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, NULL, NULL);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }

  pi_len = (int*)MALLOC(sizeof(int) * m_param * n_param);
  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, pi_len, NULL);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }
  
  LabelList = (char**)MALLOC(sizeof(char*) * m_param * n_param);
  for(i=0; i<m_param*n_param; i++)
    {
      LabelList[i] = (char*)MALLOC(sizeof(char) * (pi_len[i] + 1)); // + 1 for null termination
    }
  
  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, pi_len, LabelList);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }

  if (strcmp(LabelList[0],"fannlist") != 0) 
    {
      Scierror(999,"Argument 1 is not a fannlist\r\n");
      return NULL;
    }

  _sciErr = getPointerInList(pvApiCtx, pi_param_addr, 2, (void **)&result_ann);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }

  *res = 0;

  return result_ann;
}


//--------------------------------------------------------------------------------------------------------
struct fann_train_data * createCFannTrainDataStructFromScilabFannTrainDataStruct(unsigned int StackPos, int * res)
{
  int m_param, n_param, * pi_param_addr;
  int type, * pi_len = NULL, i;
  char ** LabelList = NULL;
  struct fann_train_data * result_ann_data = NULL;
  SciErr _sciErr;

  *res = -1;

  _sciErr = getVarAddressFromPosition(pvApiCtx, StackPos, &pi_param_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }
  _sciErr = getVarType(pvApiCtx, pi_param_addr, &type);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }

  if (type!=sci_mlist)
    {
      Scierror(999,"Argument %d is not a mlist\n",StackPos);
      return NULL;
    }

  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, NULL, NULL);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }

  pi_len = (int*)MALLOC(sizeof(int) * m_param * n_param);
  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, pi_len, NULL);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }
  
  LabelList = (char**)MALLOC(sizeof(char*) * m_param * n_param);
  for(i=0; i<m_param*n_param; i++)
    {
      LabelList[i] = (char*)MALLOC(sizeof(char) * (pi_len[i] + 1)); // + 1 for null termination
    }
  
  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, pi_len, LabelList);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }

  if (strcmp(LabelList[0],"fanntraindatalist") != 0) 
    {
      Scierror(999,"Argument 1 is not a fanntraindatalist\r\n");
      return NULL;
    }

  _sciErr = getPointerInList(pvApiCtx, pi_param_addr, 2, (void **)&result_ann_data);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }

  *res = 0;

  return result_ann_data;
}

//--------------------------------------------------------------------------------------------------------
struct fann_error * createCFannErrorStructFromScilabFannErrorStruct(unsigned int StackPos, int * res)
{
  int m_param, n_param, * pi_param_addr;
  int type, * pi_len = NULL, i;
  char ** LabelList = NULL;
  struct fann_error * result_ann_error = NULL;
  SciErr _sciErr;

  *res = -1;

  _sciErr = getVarAddressFromPosition(pvApiCtx, StackPos, &pi_param_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }
  _sciErr = getVarType(pvApiCtx, pi_param_addr, &type);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }

  if (type!=sci_mlist)
    {
      Scierror(999,"Argument %d is not a mlist\n",StackPos);
      return NULL;
    }

  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, NULL, NULL);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }

  pi_len = (int*)MALLOC(sizeof(int) * m_param * n_param);
  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, pi_len, NULL);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }
  
  LabelList = (char**)MALLOC(sizeof(char*) * m_param * n_param);
  for(i=0; i<m_param*n_param; i++)
    {
      LabelList[i] = (char*)MALLOC(sizeof(char) * (pi_len[i] + 1)); // + 1 for null termination
    }
  
  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, pi_len, LabelList);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }

  if (strcmp(LabelList[0],"fannerrorlist") != 0) 
    {
      Scierror(999,"Argument 1 is not a fannerrorlist\r\n");
      return NULL;
    }

  _sciErr = getPointerInList(pvApiCtx, pi_param_addr, 2, (void **)&result_ann_error);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }

  *res = 0;

  return result_ann_error;
}

//--------------------------------------------------------------------------------------------------------
FILE * createFILEFromScilabFannErrorStruct(unsigned int StackPos, int * res)
{
  int m_param, n_param, * pi_param_addr;
  int type, * pi_len = NULL, i;
  char ** LabelList = NULL;
  FILE * log_file = NULL;
  SciErr _sciErr;

  *res = -1;

  _sciErr = getVarAddressFromPosition(pvApiCtx, StackPos, &pi_param_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }
  _sciErr = getVarType(pvApiCtx, pi_param_addr, &type);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }

  if (type!=sci_mlist)
    {
      Scierror(999,"Argument %d is not a mlist\n",StackPos);
      return NULL;
    }

  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, NULL, NULL);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }

  pi_len = (int*)MALLOC(sizeof(int) * m_param * n_param);
  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, pi_len, NULL);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }
  
  LabelList = (char**)MALLOC(sizeof(char*) * m_param * n_param);
  for(i=0; i<m_param*n_param; i++)
    {
      LabelList[i] = (char*)MALLOC(sizeof(char) * (pi_len[i] + 1)); // + 1 for null termination
    }
  
  _sciErr = getMatrixOfStringInList(pvApiCtx, pi_param_addr, 1, &m_param, &n_param, pi_len, LabelList);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }

  if (strcmp(LabelList[0],"fannerrorlist") != 0) 
    {
      Scierror(999,"Argument 1 is not a fannerrorlist\r\n");
      return NULL;
    }

  _sciErr = getPointerInList(pvApiCtx, pi_param_addr, 2, (void **)&log_file);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return NULL;
    }

  *res = 0;

  return log_file;
}
