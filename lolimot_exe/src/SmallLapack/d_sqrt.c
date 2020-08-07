#include <math.h>
#include "f2c.h"

#ifdef USEDOUBLE
doublereal d_sqrt(doublereal *x)
{
  return( sqrt(*x) );
}
#else
doublereal d_sqrt(doublereal *x)
{
  return( sqrtf(*x) );
}
#endif
