#include <math.h>
#include "f2c.h"

#ifdef USEDOUBLE
doublereal d_lg(doublereal *x)
{
  return( log(*x) );
}
#else
doublereal d_lg(doublereal *x)
{
  return( logf(*x) );
}
#endif
