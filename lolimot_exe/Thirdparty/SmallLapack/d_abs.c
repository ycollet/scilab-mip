#include <math.h>
#include "f2c.h"

#ifdef USEDOUBLE
doublereal d_abs(doublereal x)
{
  return( fabs(x) );
}
#else
doublereal d_abs(doublereal x)
{
  return( fabsf(x) );
}
#endif
