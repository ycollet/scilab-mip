#include <serialize.h>

#include <api_scilab.h>
#include <stack-c.h>
#include <Scierror.h>
#include <sciprint.h>

static struct serialize_list internal_variable;

int sci_serialize_set(char * fname)
{
  int * p_address = NULL;
  SciErr _SciErr;

  CheckRhs(1,1);

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_address);

  serialize(p_address, &internal_variable);
  
  return 0;
}

int sci_serialize_get(char * fname)
{
  int * p_address = NULL;
  SciErr _SciErr;

  CheckLhs(1,1);
  CheckRhs(0,0);

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_address);

  unserialize(1, &internal_variable, 1, 0, 0, NULL);
  
  LhsVar(1) = 1;

  return 0;
}
