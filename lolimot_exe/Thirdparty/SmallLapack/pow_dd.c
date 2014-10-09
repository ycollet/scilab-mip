#include <math.h>
#include "f2c.h"

#ifdef USEDOUBLE
doublereal pow_dd(doublereal *ap, doublereal *bp)
{
  return(pow(*ap, *bp) );
}
#else
doublereal pow_dd(doublereal *ap, doublereal *bp)
{
  return(powf(*ap, *bp) );
}
#endif
