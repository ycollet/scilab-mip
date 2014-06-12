/*
 * Scilab interface for the FANN library
 * Author: Dirk Gorissen <dirk.gorissen@ua.ac.be>
 * Author: Yann COLLETTE <ycollet@freesurf.fr>
 * Licence: GPL version 2 or later
 */

#include <helperFann.h>

#include <stdio.h>
#include <string.h>

#include <stack-c.h>
#include <api_scilab.h>
#include <Scierror.h>

int sci_fann_read(char * fname)
{
  int * pi_name_addr = NULL;
  int res;
  char * Type = NULL;
  struct fann * result_fann = NULL;
  SciErr _sciErr;

  if (Rhs!=1)
    {
      Scierror(999,"%s: usage fann_out = %s(filename)\n", fname, fname);
      return 0;
    }

  _sciErr = getVarAddressFromPosition(pvApiCtx, 1, &pi_name_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }

  getAllocatedSingleString(pvApiCtx,  pi_name_addr, &Type);

  result_fann = fann_create_from_file(Type);
  if (result_fann==NULL)
    {
      Scierror(999,"%s: file error: %s\n",fname,Type);
      freeAllocatedSingleString(Type);
      return 0;
    }
  freeAllocatedSingleString(Type);

  //Create the struct representing this ann in Scilab
  res = createScilabFannStructFromCFannStruct(result_fann, Rhs + 1);
  if (res==-1) return 0;

  LhsVar(1) = Rhs + 1;

  return 0;
}

int sci_fann_save(char * fname)
{
  int * pi_name_addr = NULL;
  int res;
  char * Type = NULL;
  struct fann* result_fann = NULL;
  SciErr _sciErr;

  if (Rhs!=2)
    {
      Scierror(999,"%s usage: %s(ann,filename)\n", fname, fname);
      return 0;
    }

  res = detect_fannlist(1);
  if (res==-1) return 0;

  result_fann  = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (result_fann==NULL)
    {
      Scierror(999,"%s: problem while creating the fann structure\n",fname);
      return 0;
    }
  
  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_name_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }

  getAllocatedSingleString(pvApiCtx,  pi_name_addr, &Type);

  res = fann_save(result_fann, Type);
  if (res==-1)
    {
      Scierror(999,"%s: problem while saving the fann structure\n",fname);
      freeAllocatedSingleString(Type);
      return 0;
    }
  freeAllocatedSingleString(Type);

  return 0;
}

int sci_fann_savetofixed(char * fname)
{
  int * pi_name_addr = NULL;
  int res;
  char * Type = NULL;
  struct fann* result_fann = NULL;
  SciErr _sciErr;

  if (Rhs!=2)
    {
      Scierror(999,"%s usage: %s(ann,filename)\n", fname, fname);
      return 0;
    }

  res = detect_fannlist(1);
  if (res==-1) return 0;

  result_fann  = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (result_fann==NULL)
    {
      Scierror(999,"%s: Problem while creating the fann scilab structure\n",fname);
      return 0;
    }
  
  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_name_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }

  getAllocatedSingleString(pvApiCtx,  pi_name_addr, &Type);

  res = fann_save_to_fixed(result_fann, Type);
  if (res==-1)
    {
      Scierror(999,"%s: Problem while saving the fann structure\n",fname);
      freeAllocatedSingleString(Type);
      return 0;
    }
  freeAllocatedSingleString(Type);

  return 0;
}
