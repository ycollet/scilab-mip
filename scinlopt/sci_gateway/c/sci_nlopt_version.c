#include <stack-c.h>
#include <api_scilab.h>

#include <nlopt.h>

int sci_nlopt_version(char * fname)
{
  int major, minor, bugfix;
  SciErr _SciErr;

  CheckRhs(0, 0);
  CheckLhs(1, 3);

  nlopt_version(&major, &minor, &bugfix);

  createScalarDouble(pvApiCtx, 1, major);
  createScalarDouble(pvApiCtx, 2, minor);
  createScalarDouble(pvApiCtx, 3, bugfix);

  LhsVar(1) = 1;
  if (Lhs>=2) 
    {
      LhsVar(2) = 2;
    }
  if (Lhs>=3)
    {
      LhsVar(3) = 3;
    }

  return 0;
}
