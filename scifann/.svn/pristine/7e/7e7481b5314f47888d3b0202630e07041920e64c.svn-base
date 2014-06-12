/*
 * Scilab interface for the FANN library
 * Author: Dirk Gorissen <dirk.gorissen@ua.ac.be>
 * Author: Yann COLLETTE <ycollet@fresurf.fr>
 * Licence: GPL version 2 or later
 */

#include <helperFann.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include <stack-c.h>
#include <api_scilab.h>
#include <MALLOC.h>
#include <Scierror.h>

#include <fann_error.h>

// fann_set_error_log   Change where errors are logged to.
int sci_fann_set_error_log(char * fname)
{
  int * pi_name_addr = NULL;
  int res;
  char * Name     = NULL;
  FILE * log_file = NULL;
  struct fann_error * result_error = (struct fann_error *)MALLOC(1*sizeof(struct fann_error));
  SciErr _sciErr;

  // Initialisation of the structure
  result_error->errstr    = NULL;
  result_error->errno_f   = FANN_E_NO_ERROR;
  result_error->error_log = fann_default_error_log;
  
  if (Rhs==0)
    {
      fann_set_error_log(result_error,NULL);
      
      fann_reset_errno(result_error);
      fann_reset_errstr(result_error);

      res = createScilabFannErrorStructFromCFannErrorStruct(result_error,NULL,Rhs + 1);
      
      LhsVar(1) = Rhs + 1;
    }
  else
    {
      if ((Rhs!=1)&&(Lhs!=1))
	{
	  Scierror(999,"%s: usage log_out = %(filename).\n", fname, fname);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 1, &pi_name_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      
      getAllocatedSingleString(pvApiCtx,  pi_name_addr, &Name);

      log_file = fopen(Name,"w+");
      if (log_file == NULL)
	{
	  Scierror(999,"%s: unable to open the file %s for writing.\n",fname,Name);
	  freeAllocatedSingleString(Name);
	  return 0;
	}
      freeAllocatedSingleString(Name);

      fann_set_error_log(result_error,log_file);
      fann_reset_errno(result_error);
      fann_reset_errstr(result_error);

      res = createScilabFannErrorStructFromCFannErrorStruct(result_error, log_file, Rhs + 1);
      
      LhsVar(1) = Rhs + 1;
    }

  if (result_error==NULL)
    {
      Scierror(999,"%s: unable to create a fann_error structure\n",fname);
      return 0;
    }

  return 0;
}

// usage: value = fann_get_errno(ann_error);
int sci_fann_get_errno(char * fname)
{
  int res;
  double d_errno = 0.0;
  struct fann_error * result_error = NULL;

  if ((Rhs!=1)&&(Lhs!=1))
    {
      Scierror(999,"%s: usage result = %(error_fann_in).\n", fname, fname);
      return 0;
    }
  
  // Get the fann structure
  res = detect_fannerrorlist(1);
  if (res==-1) return 0;

  result_error = createCFannErrorStructFromScilabFannErrorStruct(1,&res);
  if (res==-1) return 0;

  if (result_error)
    d_errno = fann_get_errno(result_error);
  else
    d_errno = FANN_E_NO_ERROR;

  createScalarDouble(pvApiCtx, Rhs + 1, d_errno);

  LhsVar(1) = Rhs + 1;

  return 0;
}

// usage: str = fann_get_errstr(ann_error);
int sci_fann_get_errstr(char * fname)
{
  int res;
  char * ErrStr = NULL;
  struct fann_error * result_error = NULL;

  if ((Rhs!=1)&&(Lhs!=1))
    {
      Scierror(999,"%s: usage result = %(error_fann_in).\n", fname, fname);
      return 0;
    }
  
  // Get the fann structure
  res = detect_fannerrorlist(1);
  if (res==-1) return 0;

  result_error = createCFannErrorStructFromScilabFannErrorStruct(1,&res);
  if (res==-1) return 0;

  if (result_error)
    {
      if (fann_get_errno(result_error)!=FANN_E_NO_ERROR)
	{
	  ErrStr = strdup(fann_get_errstr(result_error));
	  createSingleString(pvApiCtx, Rhs + 1, ErrStr);
	  FREE(ErrStr);
	}
    }
  else
    {
      ErrStr = strdup("no error");
      createSingleString(pvApiCtx, Rhs + 1, ErrStr);
      FREE(ErrStr);
    }

  LhsVar(1) = Rhs + 1;

  return 0;
}

int sci_fann_reset_errno(char * fname)
{
  int res;
  struct fann_error * result_error = NULL;

  if (Rhs!=1)
    {
      Scierror(999,"%s: usage error_fann_out = %(error_fann_in).\n", fname, fname);
      return 0;
    }
  
  // Get the fann structure
  res = detect_fannerrorlist(1);
  if (res==-1) return 0;

  result_error = createCFannErrorStructFromScilabFannErrorStruct(1,&res);
  if (res==-1) return 0;

  if (result_error) fann_reset_errno(result_error);

  LhsVar(1) = 0;

  return 0;
}

int sci_fann_reset_errstr(char * fname)
{
  int res;
  struct fann_error * result_error = NULL;

  if (Rhs!=1)
    {
      Scierror(999,"%s: usage %(error_fann_in).\n", fname, fname);
      return 0;
    }

  // Get the fann structure
  res = detect_fannerrorlist(1);
  if (res==-1) return 0;

  result_error = createCFannErrorStructFromScilabFannErrorStruct(1,&res);
  if (res==-1) return 0;

  if (result_error) fann_reset_errstr(result_error);

  return 0;
}

int sci_fann_destroy_error_log(char * fname)
{
  int res;
  struct fann_error * result_error = NULL;
  FILE * log_file = NULL;

  if (Rhs!=1)
    {
      Scierror(999,"%s: usage %(error_fann_in).\n", fname, fname);
      return 0;
    }

  // Get the fann structure
  res = detect_fannerrorlist(1);
  if (res==-1) return 0;

  result_error = createCFannErrorStructFromScilabFannErrorStruct(1,&res);
  if (res==-1) return 0;

  if (result_error==NULL)
    {
      Scierror(999,"%s: unable to create the fann_error structure\n",fname);
      return 0;
    }

  log_file = createFILEFromScilabFannErrorStruct(1,&res);

  fclose(log_file);
  FREE(result_error);

  return 0;
}
