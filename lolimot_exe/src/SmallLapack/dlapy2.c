#include "blaswrap.h"
#include "f2c.h"

doublereal dlapy2_(doublereal *x, doublereal *y)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary   
    overflow.   

    Arguments   
    =========   

    X       (input) DOUBLE PRECISION   
    Y       (input) DOUBLE PRECISION   
            X and Y specify the values x and y.   

    ===================================================================== */
    /* System generated locals */
    doublereal ret_val, d__1;
    /* Builtin functions */
    doublereal d_sqrt(doublereal *);
    doublereal d_abs(doublereal);

    /* Local variables */
    static doublereal xabs, yabs, w, z__;
    doublereal AuxRes = 0.0;


    xabs = d_abs(*x);
    yabs = d_abs(*y);
    w = max(xabs,yabs);
    z__ = min(xabs,yabs);
    if (z__ == 0.) {
	ret_val = w;
    } else {
/* Computing 2nd power */
	d__1 = z__ / w;
	AuxRes = d__1 * d__1 + 1.;
	ret_val = w * d_sqrt(&AuxRes);
    }
    return ret_val;

/*     End of DLAPY2 */

} /* dlapy2_ */

