#include "f2c.h"

doublereal d_sign(doublereal *a, doublereal *b)
{
  doublereal x;
  x = (*a >= 0 ? *a : - *a);
  return( *b >= 0 ? x : -x);
}
