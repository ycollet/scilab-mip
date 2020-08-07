#include <math.h>
#include "f2c.h"

#ifdef USEDOUBLE
doublereal d_lg10(doublereal *x)
{
  return( log10(*x) );
}
#else
doublereal d_lg10(doublereal *x)
{
  return( log10f(*x) );
}
#endif
