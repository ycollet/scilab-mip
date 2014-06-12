/* conmin.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

#ifdef WIN32
#define Extern __declspec(dllexport)
#else
#define Extern
#endif

/* Common Block Declarations */

struct {
    doublereal delfun, dabfun, fdch, fdchm, ct, ctmin, ctl, ctlmin, alphax, 
	    abobj1, theta, obj;
    integer ndv, ncon, nside, iprint, nfdg, nscal, linobj, itmax, itrm, 
	    icndir, igoto, nac, info, infog, iter;
} Extern cnmn1_;

#define cnmn1_1 cnmn1_ 

struct {
    doublereal dm1, dm2, dm3, dm4, dm5, dm6, dm7, dm8, dm9, dm10, dm11, dm12, 
	    dct, dctl, phi, abobj, cta, ctam, ctbm, obj1, slope, dx, dx1, fi, 
	    xi, dftdf1, alp, fff, a1, a2, a3, a4, f1, f2, f3, f4, cv1, cv2, 
	    cv3, cv4, app, alpca, alpfes, alpln, alpmin, alpnc, alpsav, 
	    alpsid, alptot, rspace;
    integer idm1, idm2, idm3, jdir, iobj, kobj, kcount, ncal[2], nfeas, mscal,
	     ncobj, nvc, kount, icount, igood1, igood2, igood3, igood4, ibest,
	     iii, nlnc, jgoto, ispace[2];
} consav_;

#define consav_1 consav_

/* Table of constant values */

static integer c__1 = 1;

/* ----- CONMIN double precision version */
/* Subroutine */ int conmin_(doublereal *x, doublereal *vlb, doublereal *vub, 
	doublereal *g, doublereal *scal, doublereal *df, doublereal *a, 
	doublereal *s, doublereal *g1, doublereal *g2, doublereal *b, 
	doublereal *c__, integer *isc, integer *ic, integer *ms1, integer *n1,
	 integer *n2, integer *n3, integer *n4, integer *n5)
{
    /* Format strings */
    static char fmt_980[] = "(///5x,\002CONMIN HAS ACHIEVED A SOLUTION OF OB"
	    "J LESS THAN -1.0E+40\002/5x,\002SOLUTION APPEARS TOABE UNBOUNDE"
	    "D\002/5x,\002OPTIMIZATION IS TERMINATED\002)";
    static char fmt_1220[] = "(\0021\002,////12x,27(\002* \002)/12x,\002*"
	    "\002,51x,\002*\002/12x,\002*\002,20x,\002C O N M I N\002,20x,"
	    "\002*\002/12x,\002*\002,51x,\002*\002/12x,\002*\002,15x,\002 FOR"
	    "TRAN PROGRAM FOR \002,15x,\002*\002/12x,\002*\002,51x,\002*\002/"
	    "12x,\002*\002,9x,\002CONSTRAINED FUNCTION MINIMIZATION\002,9x"
	    ",\002*\002/12x,\002*\002,51x,\002*\002/12x,27(\002* \002))";
    static char fmt_970[] = "(///5x,\002A COMPLETELY UNCONSTRAINED FUNCTION "
	    "WITH A LINEAR OBJECTIVE IS SPECIFIED\002//10x,\002LINOBJ =\002,i"
	    "5/10x,\002NCON   =\002,i5/10x,\002NSIDE  =\002,i5//5x,\002CONTRO"
	    "L RETURNED TO CALLING PROGRAM\002)";
    static char fmt_1120[] = "(///5x,\002* * CONMIN DETECTS VLB(I).GT.VUB(I"
	    ")\002/5x,\002FIX IS SET X(I)=VLB(I)=VUB(I) = .5*(VLB(I)+VUB(I) F"
	    "OR I =\002,i5)";
    static char fmt_1130[] = "(///5x,\002* * CONMIN DETECTS INITIAL X(I).LT."
	    "VLB(I)\002/5x,\002X(I) =\002,e12.4,2x,\002VLB(I) =\002,e12.4/5x"
	    ",\002X(I) IS SET EQUAL TO VLB(I) FOR I =\002,i5)";
    static char fmt_1140[] = "(///5x,\002* * CONMIN DETECTS INITIAL X(I).GT."
	    "VUB(I)\002/5x,\002X(I) =\002,e12.4,2x,\002VUB(I) =\002,e12.4/5x"
	    ",\002X(I) IS SET EQUAL TO VUB(I) FOR I =\002,i5)";
    static char fmt_1290[] = "(////5x,\002UNCONSTRAINED FUNCTION MINIMIZAT"
	    "ION\002//5x,\002CONTROL PARAMETERS\002)";
    static char fmt_1230[] = "(////5x,\002CONSTRAINED FUNCTION MINIMIZATIO"
	    "N\002//5x,\002CONTROL PARAMETERS\002)";
    static char fmt_1240[] = "(/5x,\002IPRINT  NDV    ITMAX    NCON    NSIDE"
	    "  ICNDIR   NSCAL   NFDG\002/8i8//5x,\002LINOBJ  ITRM\002,5x,\002"
	    "N1\002,6x,\002N2\002,6x,\002N3\002,6x,\002N4\002,6x,\002N5\002/8"
	    "i8)";
    static char fmt_1260[] = "(/9x,\002CT\002,14x,\002CTMIN\002,11x,\002CT"
	    "L\002,13x,\002CTLMIN\002/1x,4(2x,e14.5)//9x,\002THETA\002,11x"
	    ",\002PHI\002,13x,\002DELFUN\002,10x,\002DABFUN\002/1x,4(2x,e14.5"
	    "))";
    static char fmt_1250[] = "(/9x,\002FDCH\002,12x,\002FDCHM\002,11x,\002AL"
	    "PHAX\002,10x,\002ABOBJ1\002/1x,4(2x,e14.5))";
    static char fmt_1270[] = "(/5x,\002LOWER BOUNDS ON DECISION VARIABLES (V"
	    "LB)\002)";
    static char fmt_1010[] = "(3x,i5,\002)\002,2x,6e13.5)";
    static char fmt_1280[] = "(/5x,\002UPPER BOUNDS ON DECISION VARIABLES (V"
	    "UB)\002)";
    static char fmt_1300[] = "(/5x,\002SCALING VECTOR (SCAL)\002)";
    static char fmt_1460[] = "(3x,7e13.4)";
    static char fmt_1020[] = "(/5x,\002LINEAR CONSTRAINT IDENTIFIERS (ISC"
	    ")\002/5x,\002NON-ZERO INDICATES LINEAR CONSTRAINT\002)";
    static char fmt_1030[] = "(3x,i5,\002)\002,2x,15i5)";
    static char fmt_1040[] = "(/5x,\002ALL CONSTRAINTS ARE LINEAR\002)";
    static char fmt_1050[] = "(/5x,\002ALL CONSTRAINTS ARE NON-LINEAR\002)";
    static char fmt_1440[] = "(//5x,\002INITIAL FUNCTION INFORMATION\002//"
	    "5x,\002OBJ =\002,e15.6)";
    static char fmt_1450[] = "(/5x,\002DECISION VARIABLES (X-VECTOR)\002)";
    static char fmt_1470[] = "(/5x,\002CONSTRAINT VALUES (G-VECTOR)\002)";
    static char fmt_1360[] = "(\002 \002)";
    static char fmt_1310[] = "(////5x,\002BEGIN ITERATION NUMBER\002,i5)";
    static char fmt_1320[] = "(/5x,\002CT =\002,e14.5,5x,\002CTL =\002,e14.5"
	    ",5x,\002PHI =\002,e14.5)";
    static char fmt_1330[] = "(/5x,\002NEW SCALING VECTOR (SCAL)\002)";
    static char fmt_1060[] = "(/5x,\002THERE ARE\002,i5,\002 ACTIVE CONSTRAI"
	    "NTS\002)";
    static char fmt_1070[] = "(5x,\002CONSTRAINT NUMBERS ARE\002)";
    static char fmt_1480[] = "(5x,15i5)";
    static char fmt_1080[] = "(/5x,\002THERE ARE\002,i5,\002 VIOLATED CONSTR"
	    "AINTS\002)";
    static char fmt_1090[] = "(/5x,\002THERE ARE\002,i5,\002 ACTIVE SIDE CON"
	    "STRAINTS\002)";
    static char fmt_1100[] = "(5x,\002DECISION VARIABLES AT LOWER OR UPPER B"
	    "OUNDS\002,\002 (MINUS INDICATES LOWER BOUND)\002)";
    static char fmt_1340[] = "(/5x,\002GRADIENT OF OBJ\002)";
    static char fmt_1350[] = "(/5x,\002GRADIENTS OF ACTIVE AND VIOLATED CONS"
	    "TRAINTS\002)";
    static char fmt_990[] = "(5x,\002CONSTRAINT NUMBER\002,i5)";
    static char fmt_1000[] = "(5x,\002SIDE CONSTRAINT ON VARIABLE\002,i5)";
    static char fmt_1370[] = "(/5x,\002PUSH-OFF FACTORS, (THETA(I), I=1,NAC"
	    ")\002)";
    static char fmt_1210[] = "(/5x,\002CONSTRAINT PARAMETER, BETA =\002,e14."
	    "5)";
    static char fmt_1380[] = "(/5x,\002SEARCH DIRECTION (S-VECTOR)\002)";
    static char fmt_1110[] = "(/5x,\002ONE-DIMENSIONAL SEARCH\002/5x,\002INI"
	    "TIAL SLOPE =\002,e12.4,2x,\002PROPOSED ALPHA =\002,e12.4)";
    static char fmt_1390[] = "(/5x,\002CALCULATED ALPHA =\002,e14.5)";
    static char fmt_1400[] = "(////5x,\002ITER =\002,i5,5x,\002OBJ =\002,e14"
	    ".5,5x,\002NO CHANGE IN OBJ\002)";
    static char fmt_1410[] = "(/5x,\002OBJ =\002,e15.6,5x,\002NO CHANGE ON O"
	    "BJ\002)";
    static char fmt_1420[] = "(/5x,\002OBJ =\002,e15.6)";
    static char fmt_1430[] = "(////5x,\002ITER =\002,i5,5x,\002OBJ =\002,e14"
	    ".5)";
    static char fmt_1490[] = "(/5x,\002THE NUMBER OF ACTIVE AND VIOLATED CON"
	    "STRAINTS EXCEEDS N3-1.\002/5x,\002DIMENSIONED SIZE OF MATRICES A"
	    " AND B AND VECTOR IC IS INSUFFICIENT\002/5x,\002OPTIMIZATION TER"
	    "MINATED AND CONTROL RETURNED TO MAIN PROGRAM.\002)";
    static char fmt_1500[] = "(\0021\002,////4x,\002FINAL OPTIMIZATION INFOR"
	    "MATION\002)";
    static char fmt_1150[] = "(/5x,\002TERMINATION CRITERION\002)";
    static char fmt_1160[] = "(10x,\002ITER EQUALS ITMAX\002)";
    static char fmt_1170[] = "(10x,\002NFEASCT CONSECUTIVE ITERATIONS FAILED"
	    " TO PRODUCE A   FEASIBLE DESIGN\002)";
    static char fmt_1180[] = "(10x,\002ABS(1-OBJ(I-1)/OBJ(I)) LESS THAN DELF"
	    "UN FOR\002,i3,\002 ITERATIONS\002)";
    static char fmt_1190[] = "(10x,\002ABS(OBJ(I)-OBJ(I-1))   LESS THAN DABF"
	    "UN FOR\002,i3,\002 ITERATIONS\002)";
    static char fmt_1200[] = "(/5x,\002NUMBER OF ITERATIONS =\002,i5)";
    static char fmt_1510[] = "(/5x,\002OBJECTIVE FUNCTION WAS EVALUATED\002,"
	    "8x,i5,2x,\002TIMES\002)";
    static char fmt_1520[] = "(/5x,\002CONSTRAINT FUNCTIONS WERE EVALUATE"
	    "D\002,i10,2x,\002TIMES\002)";
    static char fmt_1530[] = "(/5x,\002GRADIENT OF OBJECTIVE WAS CALCULATE"
	    "D\002,i9,2x,\002TIMES\002)";
    static char fmt_1540[] = "(/5x,\002GRADIENTS OF CONSTRAINTS WERE CALCULA"
	    "TED\002,i5,2x,\002TIMES\002)";

    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal c1;
    static integer m1, m2, m3;
    static doublereal x1, gi;
    static integer ii;
    static doublereal x12, si, xx, ff1, ct1;
    static integer nic;
    static doublereal sib, scj, xid;
    static integer nci;
    static doublereal ctc;
    static integer mcn1;
    static doublereal alp1;
    static integer ndv1, ndv2;
    static doublereal objb;
    static integer nnac;
    static doublereal alp11, alp12, objd;
    extern /* Subroutine */ int cnmn01_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , integer *, integer *, integer *, integer *), cnmn02_(integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    , cnmn05_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *), cnmn03_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *), cnmn06_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static integer nfeasct;

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 6, 0, fmt_980, 0 };
    static cilist io___5 = { 0, 6, 0, fmt_1220, 0 };
    static cilist io___6 = { 0, 6, 0, fmt_970, 0 };
    static cilist io___9 = { 0, 6, 0, fmt_1120, 0 };
    static cilist io___10 = { 0, 6, 0, fmt_1130, 0 };
    static cilist io___11 = { 0, 6, 0, fmt_1140, 0 };
    static cilist io___13 = { 0, 6, 0, fmt_1290, 0 };
    static cilist io___14 = { 0, 6, 0, fmt_1230, 0 };
    static cilist io___15 = { 0, 6, 0, fmt_1240, 0 };
    static cilist io___16 = { 0, 6, 0, fmt_1260, 0 };
    static cilist io___17 = { 0, 6, 0, fmt_1250, 0 };
    static cilist io___18 = { 0, 6, 0, fmt_1270, 0 };
    static cilist io___20 = { 0, 6, 0, fmt_1010, 0 };
    static cilist io___22 = { 0, 6, 0, fmt_1280, 0 };
    static cilist io___23 = { 0, 6, 0, fmt_1010, 0 };
    static cilist io___24 = { 0, 6, 0, fmt_1300, 0 };
    static cilist io___25 = { 0, 6, 0, fmt_1460, 0 };
    static cilist io___26 = { 0, 6, 0, fmt_1020, 0 };
    static cilist io___27 = { 0, 6, 0, fmt_1030, 0 };
    static cilist io___28 = { 0, 6, 0, fmt_1040, 0 };
    static cilist io___29 = { 0, 6, 0, fmt_1050, 0 };
    static cilist io___30 = { 0, 6, 0, fmt_1440, 0 };
    static cilist io___31 = { 0, 6, 0, fmt_1450, 0 };
    static cilist io___33 = { 0, 6, 0, fmt_1010, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_1470, 0 };
    static cilist io___35 = { 0, 6, 0, fmt_1010, 0 };
    static cilist io___36 = { 0, 6, 0, fmt_1360, 0 };
    static cilist io___37 = { 0, 6, 0, fmt_1310, 0 };
    static cilist io___38 = { 0, 6, 0, fmt_1320, 0 };
    static cilist io___44 = { 0, 6, 0, fmt_1330, 0 };
    static cilist io___45 = { 0, 6, 0, fmt_1460, 0 };
    static cilist io___51 = { 0, 6, 0, fmt_1060, 0 };
    static cilist io___52 = { 0, 6, 0, fmt_1070, 0 };
    static cilist io___53 = { 0, 6, 0, fmt_1480, 0 };
    static cilist io___54 = { 0, 6, 0, fmt_1080, 0 };
    static cilist io___55 = { 0, 6, 0, fmt_1070, 0 };
    static cilist io___56 = { 0, 6, 0, fmt_1480, 0 };
    static cilist io___60 = { 0, 6, 0, fmt_1090, 0 };
    static cilist io___61 = { 0, 6, 0, fmt_1100, 0 };
    static cilist io___62 = { 0, 6, 0, fmt_1480, 0 };
    static cilist io___63 = { 0, 6, 0, fmt_1340, 0 };
    static cilist io___64 = { 0, 6, 0, fmt_1010, 0 };
    static cilist io___65 = { 0, 6, 0, fmt_1350, 0 };
    static cilist io___66 = { 0, 6, 0, fmt_990, 0 };
    static cilist io___67 = { 0, 6, 0, fmt_1000, 0 };
    static cilist io___68 = { 0, 6, 0, fmt_1010, 0 };
    static cilist io___69 = { 0, 6, 0, fmt_1360, 0 };
    static cilist io___70 = { 0, 6, 0, fmt_1370, 0 };
    static cilist io___71 = { 0, 6, 0, fmt_1010, 0 };
    static cilist io___72 = { 0, 6, 0, fmt_1210, 0 };
    static cilist io___77 = { 0, 6, 0, fmt_1380, 0 };
    static cilist io___78 = { 0, 6, 0, fmt_1010, 0 };
    static cilist io___79 = { 0, 6, 0, fmt_1110, 0 };
    static cilist io___83 = { 0, 6, 0, fmt_1390, 0 };
    static cilist io___84 = { 0, 6, 0, fmt_1400, 0 };
    static cilist io___85 = { 0, 6, 0, fmt_1410, 0 };
    static cilist io___86 = { 0, 6, 0, fmt_1420, 0 };
    static cilist io___87 = { 0, 6, 0, fmt_1430, 0 };
    static cilist io___88 = { 0, 6, 0, fmt_1450, 0 };
    static cilist io___90 = { 0, 6, 0, fmt_1010, 0 };
    static cilist io___91 = { 0, 6, 0, fmt_1470, 0 };
    static cilist io___92 = { 0, 6, 0, fmt_1010, 0 };
    static cilist io___94 = { 0, 6, 0, fmt_1490, 0 };
    static cilist io___95 = { 0, 6, 0, fmt_1500, 0 };
    static cilist io___96 = { 0, 6, 0, fmt_1420, 0 };
    static cilist io___97 = { 0, 6, 0, fmt_1450, 0 };
    static cilist io___98 = { 0, 6, 0, fmt_1010, 0 };
    static cilist io___99 = { 0, 6, 0, fmt_1470, 0 };
    static cilist io___100 = { 0, 6, 0, fmt_1010, 0 };
    static cilist io___101 = { 0, 6, 0, fmt_1060, 0 };
    static cilist io___102 = { 0, 6, 0, fmt_1070, 0 };
    static cilist io___103 = { 0, 6, 0, fmt_1480, 0 };
    static cilist io___104 = { 0, 6, 0, fmt_1080, 0 };
    static cilist io___105 = { 0, 6, 0, fmt_1070, 0 };
    static cilist io___106 = { 0, 6, 0, fmt_1480, 0 };
    static cilist io___107 = { 0, 6, 0, fmt_1090, 0 };
    static cilist io___108 = { 0, 6, 0, fmt_1100, 0 };
    static cilist io___109 = { 0, 6, 0, fmt_1480, 0 };
    static cilist io___110 = { 0, 6, 0, fmt_1150, 0 };
    static cilist io___111 = { 0, 6, 0, fmt_1160, 0 };
    static cilist io___112 = { 0, 6, 0, fmt_1170, 0 };
    static cilist io___113 = { 0, 6, 0, fmt_1180, 0 };
    static cilist io___114 = { 0, 6, 0, fmt_1190, 0 };
    static cilist io___115 = { 0, 6, 0, fmt_1200, 0 };
    static cilist io___116 = { 0, 6, 0, fmt_1510, 0 };
    static cilist io___117 = { 0, 6, 0, fmt_1520, 0 };
    static cilist io___118 = { 0, 6, 0, fmt_1530, 0 };
    static cilist io___119 = { 0, 6, 0, fmt_1540, 0 };


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*      revision history */
/*       double precision version for workstations */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*      INCLUDE "stack.h" */
/* ----- Doesn't work under cygwin (include ???) */


/*     ROUTINE TO SOLVE CONSTRAINED OR UNCONSTRAINED FUNCTION */
/*     MINIMIZATION. */
/*     BY G. N. VANDERPLAATS                          APRIL, 1972. */
/*     * * * * * * * * * * *   JUNE, 1979 VERSION   * * * * * * * * * * * */
/*     NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF. */
/*     REFERENCE;  CONMIN - A FORTRAN PROGRAM FOR CONSTRAINED FUNCTION */
/*         MINIMIZATION:  USER'S MANUAL,  BY G. N. VANDERPLAATS, */
/*         NASA TM X-62,282, AUGUST, 1973. */
/*     STORAGE REQUIREMENTS: */
/*         PROGRAM - 7000 DECIMAL WORDS (CDC COMPUTER) */
/*         ARRAYS  - APPROX. 2*(NDV**2)+26*NDV+4*NCON, */
/*               WHERE N3 = NDV+2. */
/*     RE-SCALE VARIABLES IF REQUIRED. */
    /* Parameter adjustments */
    --s;
    --df;
    --scal;
    --vub;
    --vlb;
    --x;
    --isc;
    --g2;
    --g1;
    --g;
    --ic;
    b_dim1 = *n3;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    a_dim1 = *n1;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --c__;
    --ms1;

    /* Function Body */
    if (cnmn1_1.nscal == 0 || cnmn1_1.igoto == 0) {
	goto L20;
    }
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	x[i__] = c__[i__];
    }
L20:
/*     CONSTANTS. */
    ndv1 = cnmn1_1.ndv + 1;
    ndv2 = cnmn1_1.ndv + 2;
    if (cnmn1_1.igoto == 0) {
	goto L40;
    }
/*     ------------------------------------------------------------------ */
/*                     CHECK FOR UNBOUNDED SOLUTION */
/*     ------------------------------------------------------------------ */
/*     STOP IF OBJ IS LESS THAN -1.0D+40 */
    if (cnmn1_1.obj > -1e40) {
	goto L30;
    }
    s_wsfe(&io___4);
    e_wsfe();
    goto L810;
L30:
    switch (cnmn1_1.igoto) {
	case 1:  goto L160;
	case 2:  goto L390;
	case 3:  goto L380;
	case 4:  goto L670;
	case 5:  goto L690;
    }
/*     ------------------------------------------------------------------ */
/*                      SAVE INPUT CONTROL PARAMETERS */
/*     ------------------------------------------------------------------ */
L40:
    if (cnmn1_1.iprint > 0) {
	s_wsfe(&io___5);
	e_wsfe();
    }
    if (cnmn1_1.linobj == 0 || (cnmn1_1.ncon > 0 || cnmn1_1.nside > 0)) {
	goto L50;
    }
/*     TOTALLY UNCONSTRAINED FUNCTION WITH LINEAR OBJECTIVE. */
/*     SOLUTION IS UNBOUNDED. */
    s_wsfe(&io___6);
    do_fio(&c__1, (char *)&cnmn1_1.linobj, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&cnmn1_1.ncon, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&cnmn1_1.nside, (ftnlen)sizeof(integer));
    e_wsfe();
    return 0;
L50:
    consav_1.idm1 = cnmn1_1.itrm;
    consav_1.idm2 = cnmn1_1.itmax;
    consav_1.idm3 = cnmn1_1.icndir;
    consav_1.dm1 = cnmn1_1.delfun;
    consav_1.dm2 = cnmn1_1.dabfun;
    consav_1.dm3 = cnmn1_1.ct;
    consav_1.dm4 = cnmn1_1.ctmin;
    consav_1.dm5 = cnmn1_1.ctl;
    consav_1.dm6 = cnmn1_1.ctlmin;
    consav_1.dm7 = cnmn1_1.theta;
    consav_1.dm8 = consav_1.phi;
    consav_1.dm9 = cnmn1_1.fdch;
    consav_1.dm10 = cnmn1_1.fdchm;
    consav_1.dm11 = cnmn1_1.abobj1;
    consav_1.dm12 = cnmn1_1.alphax;
/*     ------------------------------------------------------------------ */
/*                                DEFAULTS */
/*     ------------------------------------------------------------------ */
    if (cnmn1_1.itrm <= 0) {
	cnmn1_1.itrm = 3;
    }
    if (cnmn1_1.itmax <= 0) {
	cnmn1_1.itmax = 20;
    }
    ndv1 = cnmn1_1.ndv + 1;
    if (cnmn1_1.icndir == 0) {
	cnmn1_1.icndir = ndv1;
    }
    if (cnmn1_1.delfun <= 0.f) {
	cnmn1_1.delfun = 1e-4f;
    }
    cnmn1_1.ct = -abs(cnmn1_1.ct);
    if (cnmn1_1.ct >= 0.f) {
	cnmn1_1.ct = -.1f;
    }
    cnmn1_1.ctmin = abs(cnmn1_1.ctmin);
    if (cnmn1_1.ctmin <= 0.f) {
	cnmn1_1.ctmin = .004f;
    }
    cnmn1_1.ctl = -abs(cnmn1_1.ctl);
    if (cnmn1_1.ctl >= 0.f) {
	cnmn1_1.ctl = -.01f;
    }
    cnmn1_1.ctlmin = abs(cnmn1_1.ctlmin);
    if (cnmn1_1.ctlmin <= 0.f) {
	cnmn1_1.ctlmin = .001f;
    }
    if (cnmn1_1.theta <= 0.f) {
	cnmn1_1.theta = 1.f;
    }
    if (cnmn1_1.abobj1 <= 0.f) {
	cnmn1_1.abobj1 = .1f;
    }
    if (cnmn1_1.alphax <= 0.f) {
	cnmn1_1.alphax = .1f;
    }
    if (cnmn1_1.fdch <= 0.f) {
	cnmn1_1.fdch = .01f;
    }
    if (cnmn1_1.fdchm <= 0.f) {
	cnmn1_1.fdchm = .01f;
    }
/*     ------------------------------------------------------------------ */
/*                     INITIALIZE INTERNAL PARAMETERS */
/*     ------------------------------------------------------------------ */
    cnmn1_1.infog = 0;
    cnmn1_1.iter = 0;
    consav_1.jdir = 0;
    consav_1.iobj = 0;
    consav_1.kobj = 0;
    ndv2 = cnmn1_1.ndv + 2;
    consav_1.kcount = 0;
    consav_1.ncal[0] = 0;
    consav_1.ncal[1] = 0;
    cnmn1_1.nac = 0;
    consav_1.nfeas = 0;
    consav_1.mscal = cnmn1_1.nscal;
    ct1 = (doublereal) cnmn1_1.itrm;
    ct1 = 1.f / ct1;
    d__1 = cnmn1_1.ctmin / abs(cnmn1_1.ct);
    consav_1.dct = pow_dd(&d__1, &ct1);
    d__1 = cnmn1_1.ctlmin / abs(cnmn1_1.ctl);
    consav_1.dctl = pow_dd(&d__1, &ct1);
    consav_1.phi = 5.f;
    consav_1.abobj = cnmn1_1.abobj1;
    consav_1.ncobj = 0;
    consav_1.ctam = abs(cnmn1_1.ctmin);
    consav_1.ctbm = abs(cnmn1_1.ctlmin);
/*     CALCULATE NUMBER OF LINEAR CONSTRAINTS, NLNC. */
    consav_1.nlnc = 0;
    if (cnmn1_1.ncon == 0) {
	goto L70;
    }
    i__1 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (isc[i__] > 0) {
	    ++consav_1.nlnc;
	}
/* L60: */
    }
L70:
/*     ------------------------------------------------------------------ */
/*          CHECK TO BE SURE THAT SIDE CONSTRAINTS ARE SATISFIED */
/*     ------------------------------------------------------------------ */
    if (cnmn1_1.nside == 0) {
	goto L110;
    }
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (vlb[i__] <= vub[i__]) {
	    goto L80;
	}
	xx = (vlb[i__] + vub[i__]) * .5f;
	x[i__] = xx;
	vlb[i__] = xx;
	vub[i__] = xx;
	s_wsfe(&io___9);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfe();
L80:
	xx = x[i__] - vlb[i__];
	if (xx >= 0.f) {
	    goto L90;
	}
/*     LOWER BOUND VIOLATED. */
	s_wsfe(&io___10);
	do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&vlb[i__], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfe();
	x[i__] = vlb[i__];
	goto L100;
L90:
	xx = vub[i__] - x[i__];
	if (xx >= 0.f) {
	    goto L100;
	}
	s_wsfe(&io___11);
	do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&vub[i__], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfe();
	x[i__] = vub[i__];
L100:
	;
    }
L110:
/*     ------------------------------------------------------------------ */
/*                        INITIALIZE SCALING VECTOR, SCAL */
/*     ------------------------------------------------------------------ */
    if (cnmn1_1.nscal == 0) {
	goto L150;
    }
    if (cnmn1_1.nscal < 0) {
	goto L130;
    }
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L120: */
	scal[i__] = 1.f;
    }
    goto L150;
L130:
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	si = (d__1 = scal[i__], abs(d__1));
	if (si < 1e-20) {
	    si = 1e-5;
	}
	scal[i__] = si;
	si = 1.f / si;
	x[i__] *= si;
	if (cnmn1_1.nside == 0) {
	    goto L140;
	}
	vlb[i__] *= si;
	vub[i__] *= si;
L140:
	;
    }
L150:
/*     ------------------------------------------------------------------ */
/*     ***** CALCULATE INITIAL FUNCTION AND CONSTRAINT VALUES  ***** */
/*     ------------------------------------------------------------------ */
    cnmn1_1.info = 1;
    consav_1.ncal[0] = 1;
    cnmn1_1.igoto = 1;
    goto L950;
L160:
    consav_1.obj1 = cnmn1_1.obj;
    if (cnmn1_1.dabfun <= 0.f) {
	cnmn1_1.dabfun = abs(cnmn1_1.obj) * .001f;
    }
    if (cnmn1_1.dabfun < 1e-10) {
	cnmn1_1.dabfun = 1e-10;
    }
    if (cnmn1_1.iprint <= 0) {
	goto L270;
    }
/*     ------------------------------------------------------------------ */
/*                    PRINT INITIAL DESIGN INFORMATION */
/*     ------------------------------------------------------------------ */
    if (cnmn1_1.iprint <= 1) {
	goto L230;
    }
    if (cnmn1_1.nside == 0 && cnmn1_1.ncon == 0) {
	s_wsfe(&io___13);
	e_wsfe();
    }
    if (cnmn1_1.nside != 0 || cnmn1_1.ncon > 0) {
	s_wsfe(&io___14);
	e_wsfe();
    }
    s_wsfe(&io___15);
    do_fio(&c__1, (char *)&cnmn1_1.iprint, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&cnmn1_1.ndv, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&cnmn1_1.itmax, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&cnmn1_1.ncon, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&cnmn1_1.nside, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&cnmn1_1.icndir, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&cnmn1_1.nscal, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&cnmn1_1.nfdg, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&cnmn1_1.linobj, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&cnmn1_1.itrm, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*n1), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*n2), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*n3), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*n4), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*n5), (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___16);
    do_fio(&c__1, (char *)&cnmn1_1.ct, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&cnmn1_1.ctmin, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&cnmn1_1.ctl, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&cnmn1_1.ctlmin, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&cnmn1_1.theta, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&consav_1.phi, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&cnmn1_1.delfun, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&cnmn1_1.dabfun, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___17);
    do_fio(&c__1, (char *)&cnmn1_1.fdch, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&cnmn1_1.fdchm, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&cnmn1_1.alphax, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&cnmn1_1.abobj1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    if (cnmn1_1.nside == 0) {
	goto L190;
    }
    s_wsfe(&io___18);
    e_wsfe();
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; i__ += 6) {
/* Computing MIN */
	i__2 = cnmn1_1.ndv, i__3 = i__ + 5;
	m1 = min(i__2,i__3);
/* L170: */
	s_wsfe(&io___20);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	i__2 = m1;
	for (j = i__; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&vlb[j], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    s_wsfe(&io___22);
    e_wsfe();
    i__2 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__2; i__ += 6) {
/* Computing MIN */
	i__1 = cnmn1_1.ndv, i__3 = i__ + 5;
	m1 = min(i__1,i__3);
/* L180: */
	s_wsfe(&io___23);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	i__1 = m1;
	for (j = i__; j <= i__1; ++j) {
	    do_fio(&c__1, (char *)&vub[j], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
L190:
    if (cnmn1_1.nscal >= 0) {
	goto L200;
    }
    s_wsfe(&io___24);
    e_wsfe();
    s_wsfe(&io___25);
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&scal[i__], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
L200:
    if (cnmn1_1.ncon == 0) {
	goto L230;
    }
    if (consav_1.nlnc == 0 || consav_1.nlnc == cnmn1_1.ncon) {
	goto L220;
    }
    s_wsfe(&io___26);
    e_wsfe();
    i__1 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__1; i__ += 15) {
/* Computing MIN */
	i__2 = cnmn1_1.ncon, i__3 = i__ + 14;
	m1 = min(i__2,i__3);
/* L210: */
	s_wsfe(&io___27);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	i__2 = m1;
	for (j = i__; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&isc[j], (ftnlen)sizeof(integer));
	}
	e_wsfe();
    }
    goto L230;
L220:
    if (consav_1.nlnc == cnmn1_1.ncon) {
	s_wsfe(&io___28);
	e_wsfe();
    }
    if (consav_1.nlnc == 0) {
	s_wsfe(&io___29);
	e_wsfe();
    }
L230:
    s_wsfe(&io___30);
    do_fio(&c__1, (char *)&cnmn1_1.obj, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___31);
    e_wsfe();
    i__2 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__2; ++i__) {
	x1 = 1.f;
	if (cnmn1_1.nscal != 0) {
	    x1 = scal[i__];
	}
/* L240: */
	g1[i__] = x[i__] * x1;
    }
    i__2 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__2; i__ += 6) {
/* Computing MIN */
	i__1 = cnmn1_1.ndv, i__3 = i__ + 5;
	m1 = min(i__1,i__3);
/* L250: */
	s_wsfe(&io___33);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	i__1 = m1;
	for (j = i__; j <= i__1; ++j) {
	    do_fio(&c__1, (char *)&g1[j], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    if (cnmn1_1.ncon == 0) {
	goto L270;
    }
    s_wsfe(&io___34);
    e_wsfe();
    i__1 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__1; i__ += 6) {
/* Computing MIN */
	i__2 = cnmn1_1.ncon, i__3 = i__ + 5;
	m1 = min(i__2,i__3);
/* L260: */
	s_wsfe(&io___35);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	i__2 = m1;
	for (j = i__; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&g[j], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
L270:
    if (cnmn1_1.iprint > 1) {
	s_wsfe(&io___36);
	e_wsfe();
    }
/*     ------------------------------------------------------------------ */
/*     ********************  BEGIN MINIMIZATION  ************************ */
/*     ------------------------------------------------------------------ */
L280:
    ++cnmn1_1.iter;
    if (cnmn1_1.abobj1 < 1e-4f) {
	cnmn1_1.abobj1 = 1e-4f;
    }
    if (cnmn1_1.abobj1 > .2f) {
	cnmn1_1.abobj1 = .2f;
    }
    if (cnmn1_1.alphax > 1.f) {
	cnmn1_1.alphax = 1.f;
    }
    if (cnmn1_1.alphax < .001f) {
	cnmn1_1.alphax = .001f;
    }

/*  THE FOLLOWING TWO LINES OF CODE WERE COMMENTED OUT ON 3/5/81 */

/*     NFEAS=NFEAS+1 */
/*     IF (NFEAS.GT.10) GO TO 810 */
    if (cnmn1_1.iprint > 2) {
	s_wsfe(&io___37);
	do_fio(&c__1, (char *)&cnmn1_1.iter, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (cnmn1_1.iprint > 3 && cnmn1_1.ncon > 0) {
	s_wsfe(&io___38);
	do_fio(&c__1, (char *)&cnmn1_1.ct, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&cnmn1_1.ctl, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&consav_1.phi, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    consav_1.cta = abs(cnmn1_1.ct);
    if (consav_1.ncobj == 0) {
	goto L340;
    }
/*     ------------------------------------------------------------------ */
/*     NO MOVE ON LAST ITERATION.  DELETE CONSTRAINTS THAT ARE NO */
/*     LONGER ACTIVE. */
/*     ------------------------------------------------------------------ */
    nnac = cnmn1_1.nac;
    i__2 = nnac;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (ic[i__] > cnmn1_1.ncon) {
	    --cnmn1_1.nac;
	}
/* L290: */
    }
    if (cnmn1_1.nac <= 0) {
	goto L420;
    }
    nnac = cnmn1_1.nac;
    i__2 = nnac;
    for (i__ = 1; i__ <= i__2; ++i__) {
L300:
	nic = ic[i__];
	ct1 = cnmn1_1.ct;
	if (isc[nic] > 0) {
	    ct1 = cnmn1_1.ctl;
	}
	if (g[nic] > ct1) {
	    goto L330;
	}
	--cnmn1_1.nac;
	if (i__ > cnmn1_1.nac) {
	    goto L420;
	}
	i__1 = cnmn1_1.nac;
	for (k = i__; k <= i__1; ++k) {
	    ii = k + 1;
	    i__3 = ndv2;
	    for (j = 1; j <= i__3; ++j) {
/* L310: */
		a[j + k * a_dim1] = a[j + ii * a_dim1];
	    }
/* L320: */
	    ic[k] = ic[ii];
	}
	goto L300;
L330:
	;
    }
    goto L420;
L340:
    if (consav_1.mscal < cnmn1_1.nscal || cnmn1_1.nscal == 0) {
	goto L360;
    }
    if (cnmn1_1.nscal < 0 && consav_1.kcount < cnmn1_1.icndir) {
	goto L360;
    }
    consav_1.mscal = 0;
    consav_1.kcount = 0;
/*     ------------------------------------------------------------------ */
/*                          SCALE VARIABLES */
/*     ------------------------------------------------------------------ */
    i__2 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__2; ++i__) {
	si = scal[i__];
	consav_1.xi = si * x[i__];
	sib = si;
	if (cnmn1_1.nscal > 0) {
	    si = abs(consav_1.xi);
	}
	if (si < 1e-10) {
	    goto L350;
	}
	scal[i__] = si;
	si = 1.f / si;
	x[i__] = consav_1.xi * si;
	if (cnmn1_1.nside == 0) {
	    goto L350;
	}
	vlb[i__] = sib * si * vlb[i__];
	vub[i__] = sib * si * vub[i__];
L350:
	;
    }
    if (cnmn1_1.iprint < 4 || cnmn1_1.nscal < 0 && cnmn1_1.iter > 1) {
	goto L360;
    }
    s_wsfe(&io___44);
    e_wsfe();
    s_wsfe(&io___45);
    i__2 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__2; ++i__) {
	do_fio(&c__1, (char *)&scal[i__], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
L360:
    ++consav_1.mscal;
    cnmn1_1.nac = 0;
/*     ------------------------------------------------------------------ */
/*          OBTAIN GRADIENTS OF OBJECTIVE AND ACTIVE CONSTRAINTS */
/*     ------------------------------------------------------------------ */
    cnmn1_1.info = 2;
    ++consav_1.ncal[1];
    if (cnmn1_1.nfdg != 1) {
	goto L370;
    }
    cnmn1_1.igoto = 2;
    goto L950;
L370:
    consav_1.jgoto = 0;
L380:
    cnmn01_(&consav_1.jgoto, &x[1], &df[1], &g[1], &isc[1], &ic[1], &a[
	    a_offset], &g1[1], &vlb[1], &vub[1], &scal[1], &c__[1], 
	    consav_1.ncal, &consav_1.dx, &consav_1.dx1, &consav_1.fi, &
	    consav_1.xi, &consav_1.iii, n1, n2, n3, n4);
    cnmn1_1.igoto = 3;
    if (consav_1.jgoto > 0) {
	goto L950;
    }
L390:
    cnmn1_1.info = 1;
    if (cnmn1_1.nac >= *n3) {
	goto L810;
    }
    if (cnmn1_1.nscal == 0 || cnmn1_1.nfdg == 0) {
	goto L420;
    }
/*     ------------------------------------------------------------------ */
/*                              SCALE GRADIENTS */
/*     ------------------------------------------------------------------ */
/*     SCALE GRADIENT OF OBJECTIVE FUNCTION. */
    i__2 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L400: */
	df[i__] *= scal[i__];
    }
    if (cnmn1_1.nfdg == 2 || cnmn1_1.nac == 0) {
	goto L420;
    }
/*     SCALE GRADIENTS OF ACTIVE CONSTRAINTS. */
    i__2 = cnmn1_1.ndv;
    for (j = 1; j <= i__2; ++j) {
	scj = scal[j];
	i__1 = cnmn1_1.nac;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L410: */
	    a[j + i__ * a_dim1] *= scj;
	}
    }
L420:
    if (cnmn1_1.iprint < 3 || cnmn1_1.ncon == 0) {
	goto L470;
    }
/*     ------------------------------------------------------------------ */
/*                                   PRINT */
/*     ------------------------------------------------------------------ */
/*     PRINT ACTIVE AND VIOLATED CONSTRAINT NUMBERS. */
    m1 = 0;
    m2 = *n3;
    if (cnmn1_1.nac == 0) {
	goto L450;
    }
    i__1 = cnmn1_1.nac;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = ic[i__];
	if (j > cnmn1_1.ncon) {
	    goto L440;
	}
	gi = g[j];
	c1 = consav_1.ctam;
	if (isc[j] > 0) {
	    c1 = consav_1.ctbm;
	}
	gi -= c1;
	if (gi > 0.f) {
	    goto L430;
	}
/*     ACTIVE CONSTRAINT. */
	++m1;
	ms1[m1] = j;
	goto L440;
L430:
	++m2;
/*     VIOLATED CONSTRAINT. */
	ms1[m2] = j;
L440:
	;
    }
L450:
    m3 = m2 - *n3;
    s_wsfe(&io___51);
    do_fio(&c__1, (char *)&m1, (ftnlen)sizeof(integer));
    e_wsfe();
    if (m1 == 0) {
	goto L460;
    }
    s_wsfe(&io___52);
    e_wsfe();
    s_wsfe(&io___53);
    i__1 = m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&ms1[i__], (ftnlen)sizeof(integer));
    }
    e_wsfe();
L460:
    s_wsfe(&io___54);
    do_fio(&c__1, (char *)&m3, (ftnlen)sizeof(integer));
    e_wsfe();
    if (m3 == 0) {
	goto L470;
    }
    s_wsfe(&io___55);
    e_wsfe();
    m3 = *n3 + 1;
    s_wsfe(&io___56);
    i__1 = m2;
    for (i__ = m3; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&ms1[i__], (ftnlen)sizeof(integer));
    }
    e_wsfe();
L470:
/*     ------------------------------------------------------------------ */
/*            CALCULATE GRADIENTS OF ACTIVE SIDE CONSTRAINTS */
/*     ------------------------------------------------------------------ */
    if (cnmn1_1.nside == 0) {
	goto L530;
    }
    mcn1 = cnmn1_1.ncon;
    m1 = 0;
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*     LOWER BOUND. */
	consav_1.xi = x[i__];
	xid = vlb[i__];
	x12 = abs(xid);
	if (x12 < 1.f) {
	    x12 = 1.f;
	}
	gi = (xid - consav_1.xi) / x12;
	if (gi < -1e-6) {
	    goto L490;
	}
	++m1;
	ms1[m1] = -i__;
	++cnmn1_1.nac;
	if (cnmn1_1.nac >= *n3) {
	    goto L810;
	}
	++mcn1;
	i__2 = cnmn1_1.ndv;
	for (j = 1; j <= i__2; ++j) {
/* L480: */
	    a[j + cnmn1_1.nac * a_dim1] = 0.f;
	}
	a[i__ + cnmn1_1.nac * a_dim1] = -1.f;
	ic[cnmn1_1.nac] = mcn1;
	g[mcn1] = gi;
	isc[mcn1] = 1;
/*     UPPER BOUND. */
L490:
	xid = vub[i__];
	x12 = abs(xid);
	if (x12 < 1.f) {
	    x12 = 1.f;
	}
	gi = (consav_1.xi - xid) / x12;
	if (gi < -1e-6) {
	    goto L510;
	}
	++m1;
	ms1[m1] = i__;
	++cnmn1_1.nac;
	if (cnmn1_1.nac >= *n3) {
	    goto L810;
	}
	++mcn1;
	i__2 = cnmn1_1.ndv;
	for (j = 1; j <= i__2; ++j) {
/* L500: */
	    a[j + cnmn1_1.nac * a_dim1] = 0.f;
	}
	a[i__ + cnmn1_1.nac * a_dim1] = 1.f;
	ic[cnmn1_1.nac] = mcn1;
	g[mcn1] = gi;
	isc[mcn1] = 1;
L510:
	;
    }
/*     ------------------------------------------------------------------ */
/*                                  PRINT */
/*     ------------------------------------------------------------------ */
/*     PRINT ACTIVE SIDE CONSTRAINT NUMBERS. */
    if (cnmn1_1.iprint < 3) {
	goto L530;
    }
    s_wsfe(&io___60);
    do_fio(&c__1, (char *)&m1, (ftnlen)sizeof(integer));
    e_wsfe();
    if (m1 == 0) {
	goto L530;
    }
    s_wsfe(&io___61);
    e_wsfe();
    s_wsfe(&io___62);
    i__1 = m1;
    for (j = 1; j <= i__1; ++j) {
	do_fio(&c__1, (char *)&ms1[j], (ftnlen)sizeof(integer));
    }
    e_wsfe();
L530:
/*     PRINT GRADIENTS OF ACTIVE AND VIOLATED CONSTRAINTS. */
    if (cnmn1_1.iprint < 4) {
	goto L570;
    }
    s_wsfe(&io___63);
    e_wsfe();
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; i__ += 6) {
/* Computing MIN */
	i__2 = cnmn1_1.ndv, i__3 = i__ + 5;
	m1 = min(i__2,i__3);
/* L540: */
	s_wsfe(&io___64);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	i__2 = m1;
	for (j = i__; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&df[j], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    if (cnmn1_1.nac == 0) {
	goto L570;
    }
    s_wsfe(&io___65);
    e_wsfe();
    i__2 = cnmn1_1.nac;
    for (i__ = 1; i__ <= i__2; ++i__) {
	m1 = ic[i__];
	m2 = m1 - cnmn1_1.ncon;
	m3 = 0;
	if (m2 > 0) {
	    m3 = (i__1 = ms1[m2], abs(i__1));
	}
	if (m2 <= 0) {
	    s_wsfe(&io___66);
	    do_fio(&c__1, (char *)&m1, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	if (m2 > 0) {
	    s_wsfe(&io___67);
	    do_fio(&c__1, (char *)&m3, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	i__1 = cnmn1_1.ndv;
	for (k = 1; k <= i__1; k += 6) {
/* Computing MIN */
	    i__3 = cnmn1_1.ndv, i__4 = k + 5;
	    m1 = min(i__3,i__4);
/* L550: */
	    s_wsfe(&io___68);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    i__3 = m1;
	    for (j = k; j <= i__3; ++j) {
		do_fio(&c__1, (char *)&a[j + i__ * a_dim1], (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsfe();
	}
/* L560: */
	s_wsfe(&io___69);
	e_wsfe();
    }
L570:
/*     ------------------------------------------------------------------ */
/*     ******************  DETERMINE SEARCH DIRECTION ******************* */
/*     ------------------------------------------------------------------ */
    consav_1.alp = 1e20;
    if (cnmn1_1.nac > 0) {
	goto L580;
    }
/*     ------------------------------------------------------------------ */
/*                        UNCONSTRAINED FUNCTION */
/*     ------------------------------------------------------------------ */
/*     FIND DIRECTION OF STEEPEST DESCENT OR CONJUGATE DIRECTION. */

/*  S. N. 575 ADDED ON 2/25/81 */

L575:
    consav_1.nvc = 0;
    consav_1.nfeas = 0;
    ++consav_1.kcount;
/*     IF KCOUNT.GT.ICNDIR  RESTART CONJUGATE DIRECTION ALGORITHM. */
    if (consav_1.kcount > cnmn1_1.icndir || consav_1.iobj == 2) {
	consav_1.kcount = 1;
    }
    if (consav_1.kcount == 1) {
	consav_1.jdir = 0;
    }
/*     IF JDIR = 0 FIND DIRECTION OF STEEPEST DESCENT. */
    cnmn02_(&consav_1.jdir, &consav_1.slope, &consav_1.dftdf1, &df[1], &s[1], 
	    n1);
    goto L630;
L580:
/*     ------------------------------------------------------------------ */
/*                          CONSTRAINED FUNCTION */
/*     ------------------------------------------------------------------ */
/*     FIND USABLE-FEASIBLE DIRECTION. */
    consav_1.kcount = 0;
    consav_1.jdir = 0;
    consav_1.phi *= 10.f;
    if (consav_1.phi > 1e3f) {
	consav_1.phi = 1e3f;
    }

/*  THE FOLLOWING LINE OF CODE WAS COMMENTED OUT ON 3/5/81 */

/*     IF (NFEAS.EQ.1) PHI=5. */
/*     CALCULATE DIRECTION, S. */
    cnmn05_(&g[1], &df[1], &a[a_offset], &s[1], &b[b_offset], &c__[1], &
	    consav_1.slope, &consav_1.phi, &isc[1], &ic[1], &ms1[1], &
	    consav_1.nvc, n1, n2, n3, n4, n5);

/*  THE FOLLOWING LINE WAS ADDED ON 2/25/81 */

    if (cnmn1_1.nac == 0) {
	goto L575;
    }

/*  THE FOLLOWING FIVE LINES WERE COMMENTED OUT ON 3/5/81 */
/*  REASON : THEY WERE NOT IN G. VANDERPLAATS LISTING */

/*     IF THIS DESIGN IS FEASIBLE AND LAST ITERATION WAS INFEASIBLE, */
/*     SET ABOBJ1=.05 (5 PERCENT). */
/*     IF (NVC.EQ.0.AND.NFEAS.GT.1) ABOBJ1=.05 */
/*     IF (NVC.EQ.0) NFEAS=0 */
    if (cnmn1_1.iprint < 3) {
	goto L600;
    }
    s_wsfe(&io___70);
    e_wsfe();
    i__2 = cnmn1_1.nac;
    for (i__ = 1; i__ <= i__2; i__ += 6) {
/* Computing MIN */
	i__3 = cnmn1_1.nac, i__1 = i__ + 5;
	m1 = min(i__3,i__1);
/* L590: */
	s_wsfe(&io___71);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	i__3 = m1;
	for (j = i__; j <= i__3; ++j) {
	    do_fio(&c__1, (char *)&a[ndv1 + j * a_dim1], (ftnlen)sizeof(
		    doublereal));
	}
	e_wsfe();
    }
    s_wsfe(&io___72);
    do_fio(&c__1, (char *)&s[ndv1], (ftnlen)sizeof(doublereal));
    e_wsfe();
L600:
/*     ------------------------------------------------------------------ */
/*     ****************** ONE-DIMENSIONAL SEARCH ************************ */
/*     ------------------------------------------------------------------ */
    if (s[ndv1] < 1e-6 && consav_1.nvc == 0) {
	goto L710;
    }
/*     ------------------------------------------------------------------ */
/*                 FIND ALPHA TO OBTAIN A FEASIBLE DESIGN */
/*     ------------------------------------------------------------------ */
    if (consav_1.nvc == 0) {
	goto L630;
    }
    consav_1.alp = -1.f;
    i__3 = cnmn1_1.nac;
    for (i__ = 1; i__ <= i__3; ++i__) {
	nci = ic[i__];
	c1 = g[nci];
	ctc = consav_1.ctam;
	if (isc[nci] > 0) {
	    ctc = consav_1.ctbm;
	}
	if (c1 <= ctc) {
	    goto L620;
	}
	alp1 = 0.f;
	i__2 = cnmn1_1.ndv;
	for (j = 1; j <= i__2; ++j) {
/* L610: */
	    alp1 += s[j] * a[j + i__ * a_dim1];
	}
	alp1 *= a[ndv2 + i__ * a_dim1];
	if (abs(alp1) < 1e-20) {
	    goto L620;
	}
	alp1 = -c1 / alp1;
	if (alp1 > consav_1.alp) {
	    consav_1.alp = alp1;
	}
L620:
	;
    }
L630:
/*     ------------------------------------------------------------------ */
/*                       LIMIT CHANCE TO ABOBJ1*OBJ */
/*     ------------------------------------------------------------------ */
    alp1 = 1e20;
    si = abs(cnmn1_1.obj);
    if (si < .01f) {
	si = .01f;
    }
    if (abs(consav_1.slope) > 1e-20) {
	alp1 = cnmn1_1.abobj1 * si / consav_1.slope;
    }
    alp1 = abs(alp1);
    if (consav_1.nvc > 0) {
	alp1 *= 10.f;
    }
    if (alp1 < consav_1.alp) {
	consav_1.alp = alp1;
    }
/*     ------------------------------------------------------------------ */
/*                   LIMIT CHANGE IN VARIABLE TO ALPHAX */
/*     ------------------------------------------------------------------ */
    alp11 = 1e20;
    i__3 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__3; ++i__) {
	si = (d__1 = s[i__], abs(d__1));
	consav_1.xi = (d__1 = x[i__], abs(d__1));
	if (si < 1e-10 || consav_1.xi < .1f) {
	    goto L640;
	}
	alp1 = cnmn1_1.alphax * consav_1.xi / si;
	if (alp1 < alp11) {
	    alp11 = alp1;
	}
L640:
	;
    }
    if (consav_1.nvc > 0) {
	alp11 *= 10.f;
    }
    if (alp11 < consav_1.alp) {
	consav_1.alp = alp11;
    }
    if (consav_1.alp > 1e20) {
	consav_1.alp = 1e20;
    }
    if (consav_1.alp <= 1e-20) {
	consav_1.alp = 1e-20;
    }
    if (cnmn1_1.iprint < 3) {
	goto L660;
    }
    s_wsfe(&io___77);
    e_wsfe();
    i__3 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__3; i__ += 6) {
/* Computing MIN */
	i__2 = cnmn1_1.ndv, i__1 = i__ + 5;
	m1 = min(i__2,i__1);
/* L650: */
	s_wsfe(&io___78);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	i__2 = m1;
	for (j = i__; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&s[j], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    s_wsfe(&io___79);
    do_fio(&c__1, (char *)&consav_1.slope, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&consav_1.alp, (ftnlen)sizeof(doublereal));
    e_wsfe();
L660:
    if (cnmn1_1.ncon > 0 || cnmn1_1.nside > 0) {
	goto L680;
    }
/*     ------------------------------------------------------------------ */
/*           DO ONE-DIMENSIONAL SEARCH FOR UNCONSTRAINED FUNCTION */
/*     ------------------------------------------------------------------ */
    consav_1.jgoto = 0;
L670:
    cnmn03_(&x[1], &s[1], &consav_1.slope, &consav_1.alp, &consav_1.fff, &
	    consav_1.a1, &consav_1.a2, &consav_1.a3, &consav_1.a4, &
	    consav_1.f1, &consav_1.f2, &consav_1.f3, &consav_1.f4, &
	    consav_1.app, n1, consav_1.ncal, &consav_1.kount, &consav_1.jgoto)
	    ;
    cnmn1_1.igoto = 4;
    if (consav_1.jgoto > 0) {
	goto L950;
    }
    consav_1.jdir = 1;
/*     PROCEED TO CONVERGENCE CHECK. */
    goto L700;
/*     ------------------------------------------------------------------ */
/*       SOLVE ONE-DIMENSIONAL SEARCH PROBLEM FOR CONSTRAINED FUNCTION */
/*     ------------------------------------------------------------------ */
L680:
    consav_1.jgoto = 0;
L690:
    cnmn06_(&x[1], &vlb[1], &vub[1], &g[1], &scal[1], &df[1], &s[1], &g1[1], &
	    g2[1], &consav_1.ctam, &consav_1.ctbm, &consav_1.slope, &
	    consav_1.alp, &consav_1.a2, &consav_1.a3, &consav_1.a4, &
	    consav_1.f1, &consav_1.f2, &consav_1.f3, &consav_1.cv1, &
	    consav_1.cv2, &consav_1.cv3, &consav_1.cv4, &consav_1.alpca, &
	    consav_1.alpfes, &consav_1.alpln, &consav_1.alpmin, &
	    consav_1.alpnc, &consav_1.alpsav, &consav_1.alpsid, &
	    consav_1.alptot, &isc[1], n1, n2, consav_1.ncal, &consav_1.nvc, &
	    consav_1.icount, &consav_1.igood1, &consav_1.igood2, &
	    consav_1.igood3, &consav_1.igood4, &consav_1.ibest, &consav_1.iii,
	     &consav_1.nlnc, &consav_1.jgoto);
    cnmn1_1.igoto = 5;
    if (consav_1.jgoto > 0) {
	goto L950;
    }
    if (cnmn1_1.nac == 0) {
	consav_1.jdir = 1;
    }
/*     ------------------------------------------------------------------ */
/*     *******************     UPDATE ALPHAX   ************************** */
/*     ------------------------------------------------------------------ */
L700:
L710:
    if (consav_1.alp > 1e19) {
	consav_1.alp = 0.f;
    }
/*     UPDATE ALPHAX TO BE AVERAGE OF MAXIMUM CHANGE IN X(I) */
/*     AND ALHPAX. */
    alp11 = 0.f;
    i__2 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__2; ++i__) {
	si = (d__1 = s[i__], abs(d__1));
	consav_1.xi = (d__1 = x[i__], abs(d__1));
	if (consav_1.xi < 1e-10) {
	    goto L720;
	}
	alp1 = consav_1.alp * si / consav_1.xi;
	if (alp1 > alp11) {
	    alp11 = alp1;
	}
L720:
	;
    }
    alp11 = (alp11 + cnmn1_1.alphax) * .5f;
    alp12 = cnmn1_1.alphax * 5.f;
    if (alp11 > alp12) {
	alp11 = alp12;
    }
    cnmn1_1.alphax = alp11;
    ++consav_1.ncobj;
/*     ABSOLUTE CHANGE IN OBJECTIVE. */
    objd = consav_1.obj1 - cnmn1_1.obj;
    objb = abs(objd);
    if (objb < 1e-10) {
	objb = 0.f;
    }
    if (cnmn1_1.nac == 0 || objb > 0.f) {
	consav_1.ncobj = 0;
    }
    if (consav_1.ncobj > 1) {
	consav_1.ncobj = 0;
    }
/*     ------------------------------------------------------------------ */
/*                                  PRINT */
/*     ------------------------------------------------------------------ */
/*     PRINT MOVE PARAMETER, NEW X-VECTOR AND CONSTRAINTS. */
    if (cnmn1_1.iprint < 3) {
	goto L730;
    }
    s_wsfe(&io___83);
    do_fio(&c__1, (char *)&consav_1.alp, (ftnlen)sizeof(doublereal));
    e_wsfe();
L730:
    if (cnmn1_1.iprint < 2) {
	goto L800;
    }
    if (objb > 0.f) {
	goto L740;
    }
    if (cnmn1_1.iprint == 2) {
	s_wsfe(&io___84);
	do_fio(&c__1, (char *)&cnmn1_1.iter, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&cnmn1_1.obj, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (cnmn1_1.iprint > 2) {
	s_wsfe(&io___85);
	do_fio(&c__1, (char *)&cnmn1_1.obj, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    goto L760;
L740:
    if (cnmn1_1.iprint == 2) {
	goto L750;
    }
    s_wsfe(&io___86);
    do_fio(&c__1, (char *)&cnmn1_1.obj, (ftnlen)sizeof(doublereal));
    e_wsfe();
    goto L760;
L750:
    s_wsfe(&io___87);
    do_fio(&c__1, (char *)&cnmn1_1.iter, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&cnmn1_1.obj, (ftnlen)sizeof(doublereal));
    e_wsfe();
L760:
    s_wsfe(&io___88);
    e_wsfe();
    i__2 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ff1 = 1.f;
	if (cnmn1_1.nscal != 0) {
	    ff1 = scal[i__];
	}
/* L770: */
	g1[i__] = ff1 * x[i__];
    }
    i__2 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__2; i__ += 6) {
/* Computing MIN */
	i__3 = cnmn1_1.ndv, i__1 = i__ + 5;
	m1 = min(i__3,i__1);
/* L780: */
	s_wsfe(&io___90);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	i__3 = m1;
	for (j = i__; j <= i__3; ++j) {
	    do_fio(&c__1, (char *)&g1[j], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    if (cnmn1_1.ncon == 0) {
	goto L800;
    }
    s_wsfe(&io___91);
    e_wsfe();
    i__3 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__3; i__ += 6) {
/* Computing MIN */
	i__2 = cnmn1_1.ncon, i__1 = i__ + 5;
	m1 = min(i__2,i__1);
/* L790: */
	s_wsfe(&io___92);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	i__2 = m1;
	for (j = i__; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&g[j], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
L800:

/*  THE FOLLOWING CODE WAS ADDED ON 3/5/81 */

/*  IT HAD NOT BEEN REPORTED AS A FIX TO MAOB */
/*  BUT WAS SENT TO JEFF STROUD A YEAR AGO */
/*  SEE OTHER COMMENTS IN CONMIN SUBROUTINE FOR DELETIONS OF CODE */
/*  ON 3/5/81 PERTAINING TO THIS FIX */


/*                   CHECK FEASIBILITY */

    if (cnmn1_1.ncon <= 0) {
	goto L808;
    }
    nfeasct = 10;
/*  added by slp 11/17/94 */
    i__2 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__2; ++i__) {
	c1 = consav_1.ctam;
	if (isc[i__] > 0) {
	    c1 = consav_1.ctbm;
	}
	if (g[i__] <= c1) {
	    goto L804;
	}
	++consav_1.nfeas;
	goto L806;
L804:
	;
    }
    if (consav_1.nfeas > 0) {
	cnmn1_1.abobj1 = .05f;
    }
/* cc */
    consav_1.nfeas = 0;
    consav_1.phi = 5.f;
L806:
    if (consav_1.nfeas >= nfeasct) {
	goto L810;
    }
L808:

/*  END OF INSERTED FIX */

/*     ------------------------------------------------------------------ */
/*                          CHECK CONVERGENCE */
/*     ------------------------------------------------------------------ */
/*     STOP IF ITER EQUALS ITMAX. */
    if (cnmn1_1.iter >= cnmn1_1.itmax) {
	goto L810;
    }
/*     ------------------------------------------------------------------ */
/*                     ABSOLUTE CHANGE IN OBJECTIVE */
/*     ------------------------------------------------------------------ */
    objb = abs(objd);
    ++consav_1.kobj;
    if (objb >= cnmn1_1.dabfun || consav_1.nfeas > 0) {
	consav_1.kobj = 0;
    }
/*     ------------------------------------------------------------------ */
/*                     RELATIVE CHANGE IN OBJECTIVE */
/*     ------------------------------------------------------------------ */
    if (abs(consav_1.obj1) > 1e-10) {
	objd /= abs(consav_1.obj1);
    }
    cnmn1_1.abobj1 = (abs(consav_1.abobj) + abs(objd)) * .5f;
    consav_1.abobj = abs(objd);
    ++consav_1.iobj;
    if (consav_1.nvc > 0 || objd >= cnmn1_1.delfun) {
	consav_1.iobj = 0;
    }
    if (consav_1.iobj >= cnmn1_1.itrm || consav_1.kobj >= cnmn1_1.itrm) {
	goto L810;
    }
    consav_1.obj1 = cnmn1_1.obj;
/*     ------------------------------------------------------------------ */
/*           REDUCE CT IF OBJECTIVE FUNCTION IS CHANGING SLOWLY */
/*     ------------------------------------------------------------------ */
    if (consav_1.iobj < 1 || cnmn1_1.nac == 0) {
	goto L280;
    }
    cnmn1_1.ct = consav_1.dct * cnmn1_1.ct;
    cnmn1_1.ctl *= consav_1.dctl;
    if (abs(cnmn1_1.ct) < cnmn1_1.ctmin) {
	cnmn1_1.ct = -cnmn1_1.ctmin;
    }
    if (abs(cnmn1_1.ctl) < cnmn1_1.ctlmin) {
	cnmn1_1.ctl = -cnmn1_1.ctlmin;
    }
    goto L280;
L810:
    if (cnmn1_1.nac >= *n3) {
	s_wsfe(&io___94);
	e_wsfe();
    }
/*     ------------------------------------------------------------------ */
/*     ****************  FINAL FUNCTION INFORMATION  ******************** */
/*     ------------------------------------------------------------------ */
    if (cnmn1_1.nscal == 0) {
	goto L830;
    }
/*     UN-SCALE THE DESIGN VARIABLES. */
    i__2 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__2; ++i__) {
	consav_1.xi = scal[i__];
	if (cnmn1_1.nside == 0) {
	    goto L820;
	}
	vlb[i__] = consav_1.xi * vlb[i__];
	vub[i__] = consav_1.xi * vub[i__];
L820:
	x[i__] = consav_1.xi * x[i__];
    }
/*     ------------------------------------------------------------------ */
/*                           PRINT FINAL RESULTS */
/*     ------------------------------------------------------------------ */
L830:
    if (cnmn1_1.iprint == 0 || cnmn1_1.nac >= *n3) {
	goto L940;
    }
    s_wsfe(&io___95);
    e_wsfe();
    s_wsfe(&io___96);
    do_fio(&c__1, (char *)&cnmn1_1.obj, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___97);
    e_wsfe();
    i__2 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__2; i__ += 6) {
/* Computing MIN */
	i__3 = cnmn1_1.ndv, i__1 = i__ + 5;
	m1 = min(i__3,i__1);
/* L840: */
	s_wsfe(&io___98);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	i__3 = m1;
	for (j = i__; j <= i__3; ++j) {
	    do_fio(&c__1, (char *)&x[j], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    if (cnmn1_1.ncon == 0) {
	goto L900;
    }
    s_wsfe(&io___99);
    e_wsfe();
    i__3 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__3; i__ += 6) {
/* Computing MIN */
	i__2 = cnmn1_1.ncon, i__1 = i__ + 5;
	m1 = min(i__2,i__1);
/* L850: */
	s_wsfe(&io___100);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	i__2 = m1;
	for (j = i__; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&g[j], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
/*     DETERMINE WHICH CONSTRAINTS ARE ACTIVE AND PRINT. */
    cnmn1_1.nac = 0;
    consav_1.nvc = 0;
    i__2 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__2; ++i__) {
	consav_1.cta = consav_1.ctam;
	if (isc[i__] > 0) {
	    consav_1.cta = consav_1.ctbm;
	}
	gi = g[i__];
	if (gi > consav_1.cta) {
	    goto L860;
	}
	if (gi < cnmn1_1.ct && isc[i__] == 0) {
	    goto L870;
	}
	if (gi < cnmn1_1.ctl && isc[i__] > 0) {
	    goto L870;
	}
	++cnmn1_1.nac;
	ic[cnmn1_1.nac] = i__;
	goto L870;
L860:
	++consav_1.nvc;
	ms1[consav_1.nvc] = i__;
L870:
	;
    }
    s_wsfe(&io___101);
    do_fio(&c__1, (char *)&cnmn1_1.nac, (ftnlen)sizeof(integer));
    e_wsfe();
    if (cnmn1_1.nac == 0) {
	goto L880;
    }
    s_wsfe(&io___102);
    e_wsfe();
    s_wsfe(&io___103);
    i__2 = cnmn1_1.nac;
    for (j = 1; j <= i__2; ++j) {
	do_fio(&c__1, (char *)&ic[j], (ftnlen)sizeof(integer));
    }
    e_wsfe();
L880:
    s_wsfe(&io___104);
    do_fio(&c__1, (char *)&consav_1.nvc, (ftnlen)sizeof(integer));
    e_wsfe();
    if (consav_1.nvc == 0) {
	goto L890;
    }
    s_wsfe(&io___105);
    e_wsfe();
    s_wsfe(&io___106);
    i__2 = consav_1.nvc;
    for (j = 1; j <= i__2; ++j) {
	do_fio(&c__1, (char *)&ms1[j], (ftnlen)sizeof(integer));
    }
    e_wsfe();
L890:
L900:
    if (cnmn1_1.nside == 0) {
	goto L930;
    }
/*     DETERMINE WHICH SIDE CONSTRAINTS ARE ACTIVE AND PRINT. */
    cnmn1_1.nac = 0;
    i__2 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__2; ++i__) {
	consav_1.xi = x[i__];
	xid = vlb[i__];
	x12 = abs(xid);
	if (x12 < 1.f) {
	    x12 = 1.f;
	}
	gi = (xid - consav_1.xi) / x12;
	if (gi < -1e-6) {
	    goto L910;
	}
	++cnmn1_1.nac;
	ms1[cnmn1_1.nac] = -i__;
L910:
	xid = vub[i__];
	x12 = abs(xid);
	if (x12 < 1.f) {
	    x12 = 1.f;
	}
	gi = (consav_1.xi - xid) / x12;
	if (gi < -1e-6) {
	    goto L920;
	}
	++cnmn1_1.nac;
	ms1[cnmn1_1.nac] = i__;
L920:
	;
    }
    s_wsfe(&io___107);
    do_fio(&c__1, (char *)&cnmn1_1.nac, (ftnlen)sizeof(integer));
    e_wsfe();
    if (cnmn1_1.nac == 0) {
	goto L930;
    }
    s_wsfe(&io___108);
    e_wsfe();
    s_wsfe(&io___109);
    i__2 = cnmn1_1.nac;
    for (j = 1; j <= i__2; ++j) {
	do_fio(&c__1, (char *)&ms1[j], (ftnlen)sizeof(integer));
    }
    e_wsfe();
L930:
    s_wsfe(&io___110);
    e_wsfe();
    if (cnmn1_1.iter >= cnmn1_1.itmax) {
	s_wsfe(&io___111);
	e_wsfe();
    }
    if (consav_1.nfeas >= nfeasct) {
	s_wsfe(&io___112);
	e_wsfe();
    }
    if (consav_1.iobj >= cnmn1_1.itrm) {
	s_wsfe(&io___113);
	do_fio(&c__1, (char *)&cnmn1_1.itrm, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (consav_1.kobj >= cnmn1_1.itrm) {
	s_wsfe(&io___114);
	do_fio(&c__1, (char *)&cnmn1_1.itrm, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    s_wsfe(&io___115);
    do_fio(&c__1, (char *)&cnmn1_1.iter, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___116);
    do_fio(&c__1, (char *)&consav_1.ncal[0], (ftnlen)sizeof(integer));
    e_wsfe();
    if (cnmn1_1.ncon > 0) {
	s_wsfe(&io___117);
	do_fio(&c__1, (char *)&consav_1.ncal[0], (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (cnmn1_1.nfdg != 0) {
	s_wsfe(&io___118);
	do_fio(&c__1, (char *)&consav_1.ncal[1], (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (cnmn1_1.ncon > 0 && cnmn1_1.nfdg == 1) {
	s_wsfe(&io___119);
	do_fio(&c__1, (char *)&consav_1.ncal[1], (ftnlen)sizeof(integer));
	e_wsfe();
    }
/*     ------------------------------------------------------------------ */
/*                   RE-SET BASIC PARAMETERS TO INPUT VALUES */
/*     ------------------------------------------------------------------ */
L940:
    cnmn1_1.itrm = consav_1.idm1;
    cnmn1_1.itmax = consav_1.idm2;
    cnmn1_1.icndir = consav_1.idm3;
    cnmn1_1.delfun = consav_1.dm1;
    cnmn1_1.dabfun = consav_1.dm2;
    cnmn1_1.ct = consav_1.dm3;
    cnmn1_1.ctmin = consav_1.dm4;
    cnmn1_1.ctl = consav_1.dm5;
    cnmn1_1.ctlmin = consav_1.dm6;
    cnmn1_1.theta = consav_1.dm7;
    consav_1.phi = consav_1.dm8;
    cnmn1_1.fdch = consav_1.dm9;
    cnmn1_1.fdchm = consav_1.dm10;
    cnmn1_1.abobj1 = consav_1.dm11;
    cnmn1_1.alphax = consav_1.dm12;
    cnmn1_1.igoto = 0;
L950:
    if (cnmn1_1.nscal == 0 || cnmn1_1.igoto == 0) {
	return 0;
    }
/*     UN-SCALE VARIABLES. */
    i__2 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__2; ++i__) {
	c__[i__] = x[i__];
/* L960: */
	x[i__] *= scal[i__];
    }
    return 0;
/*     ------------------------------------------------------------------ */
/*                                FORMATS */
/*     ------------------------------------------------------------------ */


} /* conmin_ */

/* ----- CNMN01 */
/* Subroutine */ int cnmn01_(integer *jgoto, doublereal *x, doublereal *df, 
	doublereal *g, integer *isc, integer *ic, doublereal *a, doublereal *
	g1, doublereal *vlb, doublereal *vub, doublereal *scal, doublereal *
	c__, integer *ncal, doublereal *dx, doublereal *dx1, doublereal *fi, 
	doublereal *xi, integer *iii, integer *n1, integer *n2, integer *n3, 
	integer *n4)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i__, j, i1;
    static doublereal x1;
    static integer inf;
    static doublereal fdch1;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*      revision history */
/*     double precision version for workstations */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     ROUTINE TO CALCULATE GRADIENT INFORMATION BY FINITE DIFFERENCE. */
/*     BY G. N. VANDERPLAATS                         JUNE, 1972. */
/*     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF. */
    /* Parameter adjustments */
    --ncal;
    --scal;
    --vub;
    --vlb;
    --df;
    --x;
    --g1;
    --isc;
    --g;
    a_dim1 = *n1;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ic;
    --c__;

    /* Function Body */
    if (*jgoto == 1) {
	goto L10;
    }
    if (*jgoto == 2) {
	goto L70;
    }
    cnmn1_1.infog = 0;
    inf = cnmn1_1.info;
    cnmn1_1.nac = 0;
    if (cnmn1_1.linobj != 0 && cnmn1_1.iter > 1) {
	goto L10;
    }
/*     ------------------------------------------------------------------ */
/*                    GRADIENT OF LINEAR OBJECTIVE */
/*     ------------------------------------------------------------------ */
    if (cnmn1_1.nfdg == 2) {
	*jgoto = 1;
    }
    if (cnmn1_1.nfdg == 2) {
	return 0;
    }
L10:
    *jgoto = 0;
    if (cnmn1_1.nfdg == 2 && cnmn1_1.ncon == 0) {
	return 0;
    }
    if (cnmn1_1.ncon == 0) {
	goto L40;
    }
/*     ------------------------------------------------------------------ */
/*       * * * DETERMINE WHICH CONSTRAINTS ARE ACTIVE OR VIOLATED * * * */
/*     ------------------------------------------------------------------ */
    i__1 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (g[i__] < cnmn1_1.ct) {
	    goto L20;
	}
	if (isc[i__] > 0 && g[i__] < cnmn1_1.ctl) {
	    goto L20;
	}
	++cnmn1_1.nac;
	if (cnmn1_1.nac >= *n3) {
	    return 0;
	}
	ic[cnmn1_1.nac] = i__;
L20:
	;
    }
    if (cnmn1_1.nfdg == 2 && cnmn1_1.nac == 0) {
	return 0;
    }
    if (cnmn1_1.linobj > 0 && cnmn1_1.iter > 1 && cnmn1_1.nac == 0) {
	return 0;
    }
/*     ------------------------------------------------------------------ */
/*                  STORE VALUES OF CONSTRAINTS IN G1 */
/*     ------------------------------------------------------------------ */
    i__1 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
	g1[i__] = g[i__];
    }
L40:
    *jgoto = 0;
    if (cnmn1_1.nac == 0 && cnmn1_1.nfdg == 2) {
	return 0;
    }
/*     ------------------------------------------------------------------ */
/*                            CALCULATE GRADIENTS */
/*     ------------------------------------------------------------------ */
    cnmn1_1.infog = 1;
    cnmn1_1.info = 1;
    *fi = cnmn1_1.obj;
    *iii = 0;
L50:
    ++(*iii);
    *xi = x[*iii];
    *dx = cnmn1_1.fdch * *xi;
    *dx = abs(*dx);
    fdch1 = cnmn1_1.fdchm;
    if (cnmn1_1.nscal != 0) {
	fdch1 = cnmn1_1.fdchm / scal[*iii];
    }
    if (*dx < fdch1) {
	*dx = fdch1;
    }
    x1 = *xi + *dx;
    if (cnmn1_1.nside == 0) {
	goto L60;
    }
    if (x1 > vub[*iii]) {
	*dx = -(*dx);
    }
L60:
    *dx1 = 1.f / *dx;
    x[*iii] = *xi + *dx;
    ++ncal[1];
/*     ------------------------------------------------------------------ */
/*                         FUNCTION EVALUATION */
/*     ------------------------------------------------------------------ */
    *jgoto = 2;
    return 0;
L70:
    x[*iii] = *xi;
    if (cnmn1_1.nfdg == 0) {
	df[*iii] = *dx1 * (cnmn1_1.obj - *fi);
    }
    if (cnmn1_1.nac == 0) {
	goto L90;
    }
/*     ------------------------------------------------------------------ */
/*             DETERMINE GRADIENT COMPONENTS OF ACTIVE CONSTRAINTS */
/*     ------------------------------------------------------------------ */
    i__1 = cnmn1_1.nac;
    for (j = 1; j <= i__1; ++j) {
	i1 = ic[j];
/* L80: */
	a[*iii + j * a_dim1] = *dx1 * (g[i1] - g1[i1]);
    }
L90:
    if (*iii < cnmn1_1.ndv) {
	goto L50;
    }
    cnmn1_1.infog = 0;
    cnmn1_1.info = inf;
    *jgoto = 0;
    cnmn1_1.obj = *fi;
    if (cnmn1_1.ncon == 0) {
	return 0;
    }
/*     ------------------------------------------------------------------ */
/*             STORE CURRENT CONSTRAINT VALUES BACK IN G-VECTOR */
/*     ------------------------------------------------------------------ */
    i__1 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L100: */
	g[i__] = g1[i__];
    }
    return 0;
} /* cnmn01_ */

/* ----- CNMN02 */
/* Subroutine */ int cnmn02_(integer *ncalc, doublereal *slope, doublereal *
	dftdf1, doublereal *df, doublereal *s, integer *n1)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static doublereal s1, s2, si, dfi, beta, dftdf;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*      revision history */
/*     double precision version for workstations */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     ROUTINE TO DETERMINE CONJUGATE DIRECTION VECTOR OR DIRECTION */
/*     OF STEEPEST DESCENT FOR UNCONSTRAINED FUNCTION MINIMIZATION. */
/*     BY G. N. VANDERPLAATS                       APRIL, 1972. */
/*     NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF. */
/*     NCALC = CALCULATION CONTROL. */
/*         NCALC = 0,     S = STEEPEST DESCENT. */
/*         NCALC = 1,     S = CONJUGATE DIRECTION. */
/*     CONJUGATE DIRECTION IS FOUND BY FLETCHER-REEVES ALGORITHM. */
/*     ------------------------------------------------------------------ */
/*                   CALCULATE NORM OF GRADIENT VECTOR */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --s;
    --df;

    /* Function Body */
    dftdf = 0.f;
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dfi = df[i__];
/* L10: */
	dftdf += dfi * dfi;
    }
/*     ------------------------------------------------------------------ */
/*     **********                FIND DIRECTION S              ********** */
/*     ------------------------------------------------------------------ */
    if (*ncalc != 1) {
	goto L30;
    }
    if (*dftdf1 < 1e-20) {
	goto L30;
    }
/*     ------------------------------------------------------------------ */
/*                 FIND FLETCHER-REEVES CONJUGATE DIRECTION */
/*     ------------------------------------------------------------------ */
    beta = dftdf / *dftdf1;
    *slope = 0.f;
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dfi = df[i__];
	si = beta * s[i__] - dfi;
	*slope += si * dfi;
/* L20: */
	s[i__] = si;
    }
    goto L50;
L30:
    *ncalc = 0;
/*     ------------------------------------------------------------------ */
/*                  CALCULATE DIRECTION OF STEEPEST DESCENT */
/*     ------------------------------------------------------------------ */
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L40: */
	s[i__] = -df[i__];
    }
    *slope = -dftdf;
L50:
/*     ------------------------------------------------------------------ */
/*                  NORMALIZE S TO MAX ABS VALUE OF UNITY */
/*     ------------------------------------------------------------------ */
    s1 = 0.f;
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s2 = (d__1 = s[i__], abs(d__1));
	if (s2 > s1) {
	    s1 = s2;
	}
/* L60: */
    }
    if (s1 < 1e-20) {
	s1 = 1e-20;
    }
    s1 = 1.f / s1;
    *dftdf1 = dftdf * s1;
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L70: */
	s[i__] = s1 * s[i__];
    }
    *slope = s1 * *slope;
    return 0;
} /* cnmn02_ */

/* ----- CNMN03 */
/* Subroutine */ int cnmn03_(doublereal *x, doublereal *s, doublereal *slope, 
	doublereal *alp, doublereal *fff, doublereal *a1, doublereal *a2, 
	doublereal *a3, doublereal *a4, doublereal *f1, doublereal *f2, 
	doublereal *f3, doublereal *f4, doublereal *app, integer *n1, integer 
	*ncal, integer *kount, integer *jgoto)
{
    /* Format strings */
    static char fmt_360[] = "(/////5x,\002* * * UNCONSTRAINED ONE-DIMENSIONA"
	    "L SEARCH INFORMATION * * *\002)";
    static char fmt_370[] = "(/5x,\002ALPHA =\002,e14.5/5x,\002X-VECTOR\002)";
    static char fmt_380[] = "(5x,6e13.5)";
    static char fmt_390[] = "(/5x,\002OBJ =\002,e14.5)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__;
    static doublereal aa, ab, ff, ap;
    static integer ii, jj;
    static doublereal ab2, ab3, ap1, zro;
    extern /* Subroutine */ int cnmn04_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);

    /* Fortran I/O blocks */
    static cilist io___134 = { 0, 6, 0, fmt_360, 0 };
    static cilist io___138 = { 0, 6, 0, fmt_370, 0 };
    static cilist io___139 = { 0, 6, 0, fmt_380, 0 };
    static cilist io___140 = { 0, 6, 0, fmt_390, 0 };
    static cilist io___142 = { 0, 6, 0, fmt_370, 0 };
    static cilist io___143 = { 0, 6, 0, fmt_380, 0 };
    static cilist io___144 = { 0, 6, 0, fmt_390, 0 };
    static cilist io___147 = { 0, 6, 0, fmt_370, 0 };
    static cilist io___148 = { 0, 6, 0, fmt_380, 0 };
    static cilist io___149 = { 0, 6, 0, fmt_390, 0 };
    static cilist io___150 = { 0, 6, 0, fmt_370, 0 };
    static cilist io___151 = { 0, 6, 0, fmt_380, 0 };
    static cilist io___152 = { 0, 6, 0, fmt_390, 0 };
    static cilist io___153 = { 0, 6, 0, fmt_370, 0 };
    static cilist io___154 = { 0, 6, 0, fmt_380, 0 };
    static cilist io___155 = { 0, 6, 0, fmt_390, 0 };
    static cilist io___160 = { 0, 6, 0, fmt_370, 0 };
    static cilist io___161 = { 0, 6, 0, fmt_380, 0 };
    static cilist io___162 = { 0, 6, 0, fmt_390, 0 };
    static cilist io___163 = { 0, 6, 0, fmt_370, 0 };
    static cilist io___164 = { 0, 6, 0, fmt_380, 0 };
    static cilist io___165 = { 0, 6, 0, fmt_390, 0 };


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*      revision history */
/*     double precision version for workstations */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     ROUTINE TO SOLVE ONE-DIMENSIONAL SEARCH IN UNCONSTRAINED */
/*     MINIMIZATION USING 2-POINT QUADRATIC INTERPOLATION, 3-POINT */
/*     CUBIC INTERPOLATION AND 4-POINT CUBIC INTERPOLATION. */
/*     BY G. N. VANDERPLAATS                         APRIL, 1972. */
/*     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF. */
/*     ALP = PROPOSED MOVE PARAMETER. */
/*     SLOPE = INITIAL FUNCTION SLOPE = S-TRANSPOSE TIMES DF. */
/*     SLOPE MUST BE NEGATIVE. */
/*     OBJ = INITIAL FUNCTION VALUE. */
    /* Parameter adjustments */
    --s;
    --x;
    --ncal;

    /* Function Body */
    zro = 0.f;
    if (*jgoto == 0) {
	goto L10;
    }
    switch (*jgoto) {
	case 1:  goto L50;
	case 2:  goto L80;
	case 3:  goto L110;
	case 4:  goto L140;
	case 5:  goto L180;
	case 6:  goto L220;
	case 7:  goto L270;
    }
/*     ------------------------------------------------------------------ */
/*                     INITIAL INFORMATION  (ALPHA=0) */
/*     ------------------------------------------------------------------ */
L10:
    if (*slope < 0.f) {
	goto L20;
    }
    *alp = 0.f;
    return 0;
L20:
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___134);
	e_wsfe();
    }
    *fff = cnmn1_1.obj;
    ap1 = 0.f;
    *a1 = 0.f;
    *f1 = cnmn1_1.obj;
    *a2 = *alp;
    *a3 = 0.f;
    *f3 = 0.f;
    ap = *a2;
    *kount = 0;
/*     ------------------------------------------------------------------ */
/*            MOVE A DISTANCE AP*S AND UPDATE FUNCTION VALUE */
/*     ------------------------------------------------------------------ */
L30:
    ++(*kount);
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L40: */
	x[i__] += ap * s[i__];
    }
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___138);
	do_fio(&c__1, (char *)&ap, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___139);
	i__1 = cnmn1_1.ndv;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    ++ncal[1];
    *jgoto = 1;
    return 0;
L50:
    *f2 = cnmn1_1.obj;
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___140);
	do_fio(&c__1, (char *)&(*f2), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (*f2 < *f1) {
	goto L120;
    }
/*     ------------------------------------------------------------------ */
/*                     CHECK FOR ILL-CONDITIONING */
/*     ------------------------------------------------------------------ */
    if (*kount > 5) {
	goto L60;
    }
    ff = abs(*f1) * 2.f;
    if (*f2 < ff) {
	goto L90;
    }
    ff = abs(*f1) * 5.f;
    if (*f2 < ff) {
	goto L60;
    }
    *a2 *= .5f;
    ap = -(*a2);
    *alp = *a2;
    goto L30;
L60:
    *f3 = *f2;
    *a3 = *a2;
    *a2 *= .5f;
/*     ------------------------------------------------------------------ */
/*                 UPDATE DESIGN VECTOR AND FUNCTION VALUE */
/*     ------------------------------------------------------------------ */
    ap = *a2 - *alp;
    *alp = *a2;
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L70: */
	x[i__] += ap * s[i__];
    }
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___142);
	do_fio(&c__1, (char *)&(*a2), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___143);
	i__1 = cnmn1_1.ndv;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    ++ncal[1];
    *jgoto = 2;
    return 0;
L80:
    *f2 = cnmn1_1.obj;
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___144);
	do_fio(&c__1, (char *)&(*f2), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
/*     PROCEED TO CUBIC INTERPOLATION. */
    goto L160;
L90:
/*     ------------------------------------------------------------------ */
/*     **********        2-POINT QUADRATIC INTERPOLATION       ********** */
/*     ------------------------------------------------------------------ */
    jj = 1;
    ii = 1;
    cnmn04_(&ii, app, &zro, a1, f1, slope, a2, f2, &zro, &zro, &zro, &zro);
    if (*app < zro || *app > *a2) {
	goto L120;
    }
    *f3 = *f2;
    *a3 = *a2;
    *a2 = *app;
    jj = 0;
/*     ------------------------------------------------------------------ */
/*                  UPDATE DESIGN VECTOR AND FUNCTION VALUE */
/*     ------------------------------------------------------------------ */
    ap = *a2 - *alp;
    *alp = *a2;
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L100: */
	x[i__] += ap * s[i__];
    }
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___147);
	do_fio(&c__1, (char *)&(*a2), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___148);
	i__1 = cnmn1_1.ndv;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    ++ncal[1];
    *jgoto = 3;
    return 0;
L110:
    *f2 = cnmn1_1.obj;
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___149);
	do_fio(&c__1, (char *)&(*f2), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    goto L150;
L120:
    *a3 = *a2 * 2.f;
/*     ------------------------------------------------------------------ */
/*                  UPDATE DESIGN VECTOR AND FUNCTION VALUE */
/*     ------------------------------------------------------------------ */
    ap = *a3 - *alp;
    *alp = *a3;
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L130: */
	x[i__] += ap * s[i__];
    }
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___150);
	do_fio(&c__1, (char *)&(*a3), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___151);
	i__1 = cnmn1_1.ndv;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    ++ncal[1];
    *jgoto = 4;
    return 0;
L140:
    *f3 = cnmn1_1.obj;
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___152);
	do_fio(&c__1, (char *)&(*f3), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
L150:
    if (*f3 < *f2) {
	goto L190;
    }
L160:
/*     ------------------------------------------------------------------ */
/*     **********       3-POINT CUBIC INTERPOLATION      ********** */
/*     ------------------------------------------------------------------ */
    ii = 3;
    cnmn04_(&ii, app, &zro, a1, f1, slope, a2, f2, a3, f3, &zro, &zro);
    if (*app < zro || *app > *a3) {
	goto L190;
    }
/*     ------------------------------------------------------------------ */
/*     UPDATE DESIGN VECTOR AND FUNCTION VALUE. */
/*     ------------------------------------------------------------------ */
    ap1 = *app;
    ap = *app - *alp;
    *alp = *app;
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L170: */
	x[i__] += ap * s[i__];
    }
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___153);
	do_fio(&c__1, (char *)&(*alp), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___154);
	i__1 = cnmn1_1.ndv;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    ++ncal[1];
    *jgoto = 5;
    return 0;
L180:
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___155);
	do_fio(&c__1, (char *)&cnmn1_1.obj, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
/*     ------------------------------------------------------------------ */
/*                         CHECK CONVERGENCE */
/*     ------------------------------------------------------------------ */
    aa = 1.f - *app / *a2;
    ab2 = abs(*f2);
    ab3 = abs(cnmn1_1.obj);
    ab = ab2;
    if (ab3 > ab) {
	ab = ab3;
    }
    if (ab < 1e-15) {
	ab = 1e-15;
    }
    ab = (ab2 - ab3) / ab;
    if (abs(ab) < 1e-15 && abs(aa) < .001f) {
	goto L330;
    }
    *a4 = *a3;
    *f4 = *f3;
    *a3 = *app;
    *f3 = cnmn1_1.obj;
    if (*a3 > *a2) {
	goto L230;
    }
    *a3 = *a2;
    *f3 = *f2;
    *a2 = *app;
    *f2 = cnmn1_1.obj;
    goto L230;
L190:
/*     ------------------------------------------------------------------ */
/*     **********        4-POINT CUBIC INTERPOLATION       ********** */
/*     ------------------------------------------------------------------ */
L200:
    *a4 = *a3 * 2.f;
/*     UPDATE DESIGN VECTOR AND FUNCTION VALUE. */
    ap = *a4 - *alp;
    *alp = *a4;
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L210: */
	x[i__] += ap * s[i__];
    }
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___160);
	do_fio(&c__1, (char *)&(*alp), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___161);
	i__1 = cnmn1_1.ndv;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    ++ncal[1];
    *jgoto = 6;
    return 0;
L220:
    *f4 = cnmn1_1.obj;
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___162);
	do_fio(&c__1, (char *)&(*f4), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (*f4 > *f3) {
	goto L230;
    }
    *a1 = *a2;
    *f1 = *f2;
    *a2 = *a3;
    *f2 = *f3;
    *a3 = *a4;
    *f3 = *f4;
    goto L200;
L230:
    ii = 4;
    cnmn04_(&ii, app, a1, a1, f1, slope, a2, f2, a3, f3, a4, f4);
    if (*app > *a1) {
	goto L250;
    }
    ap = *a1 - *alp;
    *alp = *a1;
    cnmn1_1.obj = *f1;
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L240: */
	x[i__] += ap * s[i__];
    }
    goto L280;
L250:
/*     ------------------------------------------------------------------ */
/*                 UPDATE DESIGN VECTOR AND FUNCTION VALUE */
/*     ------------------------------------------------------------------ */
    ap = *app - *alp;
    *alp = *app;
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L260: */
	x[i__] += ap * s[i__];
    }
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___163);
	do_fio(&c__1, (char *)&(*alp), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___164);
	i__1 = cnmn1_1.ndv;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    ++ncal[1];
    *jgoto = 7;
    return 0;
L270:
    if (cnmn1_1.iprint > 4) {
	s_wsfe(&io___165);
	do_fio(&c__1, (char *)&cnmn1_1.obj, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
L280:
/*     ------------------------------------------------------------------ */
/*                    CHECK FOR ILL-CONDITIONING */
/*     ------------------------------------------------------------------ */
    if (cnmn1_1.obj > *f2 || cnmn1_1.obj > *f3) {
	goto L290;
    }
    if (cnmn1_1.obj <= *f1) {
	goto L330;
    }
    ap = *a1 - *alp;
    *alp = *a1;
    cnmn1_1.obj = *f1;
    goto L310;
L290:
    if (*f2 < *f3) {
	goto L300;
    }
    cnmn1_1.obj = *f3;
    ap = *a3 - *alp;
    *alp = *a3;
    goto L310;
L300:
    cnmn1_1.obj = *f2;
    ap = *a2 - *alp;
    *alp = *a2;
L310:
/*     ------------------------------------------------------------------ */
/*                       UPDATE DESIGN VECTOR */
/*     ------------------------------------------------------------------ */
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L320: */
	x[i__] += ap * s[i__];
    }
L330:
/*     ------------------------------------------------------------------ */
/*                     CHECK FOR MULTIPLE MINIMA */
/*     ------------------------------------------------------------------ */
    if (cnmn1_1.obj <= *fff) {
	goto L350;
    }
/*     INITIAL FUNCTION IS MINIMUM. */
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L340: */
	x[i__] -= *alp * s[i__];
    }
    *alp = 0.f;
    cnmn1_1.obj = *fff;
L350:
    *jgoto = 0;
    return 0;
/*     ------------------------------------------------------------------ */
/*                                 FORMATS */
/*     ------------------------------------------------------------------ */


} /* cnmn03_ */

/* ----- CNMN04 */
/* Subroutine */ int cnmn04_(integer *ii, doublereal *xbar, doublereal *eps, 
	doublereal *x1, doublereal *y1, doublereal *slope, doublereal *x2, 
	doublereal *y2, doublereal *x3, doublereal *y3, doublereal *x4, 
	doublereal *y4)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal q1, q2, q3, q4, q5, q6, aa, bb, cc, x11, x21, dx, x31, 
	    x32, x41, x42, x22, qq, x33, x44, x111, x222, bac, dnom, xbar1;
    static integer nslop;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*      revision history */
/*     double precision version for workstations */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     ROUTINE TO FIND FIRST XBAR.GE.EPS CORRESPONDING TO A MINIMUM */
/*     OF A ONE-DIMENSIONAL REAL FUNCTION BY POLYNOMIEL INTERPOLATION. */
/*     BY G. N. VANDERPLAATS                          APRIL, 1972. */
/*     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF. */

/*     II = CALCULATION CONTROL. */
/*          1:  2-POINT QUADRATIC INTERPOLATION, GIVEN X1, Y1, SLOPE, */
/*              X2 AND Y2. */
/*          2:  3-POINT QUADRATIC INTERPOLATION, GIVEN X1, Y1, X2, Y2, */
/*              X3 AND Y3. */
/*          3:  3-POINT CUBIC INTERPOLATION, GIVEN X1, Y1, SLOPE, X2, Y2, */
/*              X3 AND Y3. */
/*          4:  4-POINT CUBIC INTERPOLATION, GIVEN X1, Y1, X2, Y2, X3, */
/*              Y3, X4 AND Y4. */
/*     EPS MAY BE NEGATIVE. */
/*     IF REQUIRED MINIMUM ON Y DOES NOT EXITS, OR THE FUNCTION IS */
/*     ILL-CONDITIONED, XBAR = EPS-1.0 WILL BE RETURNED AS AN ERROR */
/*     INDICATOR. */
/*     IF DESIRED INTERPOLATION IS ILL-CONDITIONED, A LOWER ORDER */
/*     INTERPOLATION, CONSISTANT WITH INPUT DATA, WILL BE ATTEMPTED, */
/*     AND II WILL BE CHANGED ACCORDINGLY. */
    xbar1 = *eps - 1.f;
    *xbar = xbar1;
    x21 = *x2 - *x1;
    if (abs(x21) < 1e-20) {
	return 0;
    }
    nslop = *ii % 2;
    switch (*ii) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L40;
	case 4:  goto L50;
    }
L10:
/*     ------------------------------------------------------------------ */
/*                 II=1: 2-POINT QUADRATIC INTERPOLATION */
/*     ------------------------------------------------------------------ */
    *ii = 1;
    dx = *x1 - *x2;
    if (abs(dx) < 1e-20) {
	return 0;
    }
    aa = (*slope + (*y2 - *y1) / dx) / dx;
    if (aa < 1e-20) {
	return 0;
    }
    bb = *slope - aa * 2.f * *x1;
    *xbar = bb * -.5f / aa;
    if (*xbar < *eps) {
	*xbar = xbar1;
    }
    return 0;
L20:
/*     ------------------------------------------------------------------ */
/*                 II=2: 3-POINT QUADRATIC INTERPOLATION */
/*     ------------------------------------------------------------------ */
    *ii = 2;
    x21 = *x2 - *x1;
    x31 = *x3 - *x1;
    x32 = *x3 - *x2;
    qq = x21 * x31 * x32;
    if (abs(qq) < 1e-20) {
	return 0;
    }
    aa = (*y1 * x32 - *y2 * x31 + *y3 * x21) / qq;
    if (aa < 1e-20) {
	goto L30;
    }
    bb = (*y2 - *y1) / x21 - aa * (*x1 + *x2);
    *xbar = bb * -.5f / aa;
    if (*xbar < *eps) {
	*xbar = xbar1;
    }
    return 0;
L30:
    if (nslop == 0) {
	return 0;
    }
    goto L10;
L40:
/*     ------------------------------------------------------------------ */
/*                   II=3: 3-POINT CUBIC INTERPOLATION */
/*     ------------------------------------------------------------------ */
    *ii = 3;
    x21 = *x2 - *x1;
    x31 = *x3 - *x1;
    x32 = *x3 - *x2;
    qq = x21 * x31 * x32;
    if (abs(qq) < 1e-20) {
	return 0;
    }
    x11 = *x1 * *x1;
    dnom = *x2 * *x2 * x31 - x11 * x32 - *x3 * *x3 * x21;
    if (abs(dnom) < 1e-20) {
	goto L20;
    }
    aa = ((x31 * x31 * (*y2 - *y1) - x21 * x21 * (*y3 - *y1)) / (x31 * x21) - 
	    *slope * x32) / dnom;
    if (abs(aa) < 1e-20) {
	goto L20;
    }
    bb = ((*y2 - *y1) / x21 - *slope - aa * (*x2 * *x2 + *x1 * *x2 - x11 * 
	    2.f)) / x21;
    cc = *slope - aa * 3.f * x11 - bb * 2.f * *x1;
    bac = bb * bb - aa * 3.f * cc;
    if (bac < 0.f) {
	goto L20;
    }
    bac = sqrt(bac);
    *xbar = (bac - bb) / (aa * 3.f);
    if (*xbar < *eps) {
	*xbar = *eps;
    }
    return 0;
L50:
/*     ------------------------------------------------------------------ */
/*                    II=4: 4-POINT CUBIC INTERPOLATION */
/*     ------------------------------------------------------------------ */
    x21 = *x2 - *x1;
    x31 = *x3 - *x1;
    x41 = *x4 - *x1;
    x32 = *x3 - *x2;
    x42 = *x4 - *x2;
    x11 = *x1 * *x1;
    x22 = *x2 * *x2;
    x33 = *x3 * *x3;
    x44 = *x4 * *x4;
    x111 = *x1 * x11;
    x222 = *x2 * x22;
    q2 = x31 * x21 * x32;
    if (abs(q2) < 1e-30) {
	return 0;
    }
    q1 = x111 * x32 - x222 * x31 + *x3 * x33 * x21;
    q4 = x111 * x42 - x222 * x41 + *x4 * x44 * x21;
    q5 = x41 * x21 * x42;
    dnom = q2 * q4 - q1 * q5;
    if (abs(dnom) < 1e-30) {
	goto L60;
    }
    q3 = *y3 * x21 - *y2 * x31 + *y1 * x32;
    q6 = *y4 * x21 - *y2 * x41 + *y1 * x42;
    aa = (q2 * q6 - q3 * q5) / dnom;
    bb = (q3 - q1 * aa) / q2;
    cc = (*y2 - *y1 - aa * (x222 - x111)) / x21 - bb * (*x1 + *x2);
    bac = bb * bb - aa * 3.f * cc;
    if (abs(aa) < 1e-20 || bac < 0.f) {
	goto L60;
    }
    bac = sqrt(bac);
    *xbar = (bac - bb) / (aa * 3.f);
    if (*xbar < *eps) {
	*xbar = xbar1;
    }
    return 0;
L60:
    if (nslop == 1) {
	goto L40;
    }
    goto L20;
} /* cnmn04_ */

/* ----- CNMN05 */
/* Subroutine */ int cnmn05_(doublereal *g, doublereal *df, doublereal *a, 
	doublereal *s, doublereal *b, doublereal *c__, doublereal *slope, 
	doublereal *phi, integer *isc, integer *ic, integer *ms1, integer *
	nvc, integer *n1, integer *n2, integer *n3, integer *n4, integer *n5)
{
    /* Format strings */
    static char fmt_43[] = "(5x,\002** CONSTRAINT\002,i5,\002 HAS ZERO GRADI"
	    "ENT\002/5x,\002DELETED FROM ACTIVE SET\002)";
    static char fmt_180[] = "(//5x,\002* * DIRECTION FINDING PROCESS DID NOT"
	    " CONVERGE\002/5x,\002* * S-VECTOR MAY NOT BE VALID\002)";
    static char fmt_178[] = "(5x,\002** CALCULATED S-VECTOR IS NOT FEASIBL"
	    "E\002/5x,\002BETA IS SET TO ZERO\002)";

    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, j, k;
    static doublereal a1, c1;
    static integer j1;
    static doublereal s1, gg, sg, ct1, ct2;
    static integer ndb;
    static doublereal cta, ctb;
    static integer nci, ncj;
    static doublereal ctd, ctc;
    static integer ner;
    static doublereal tht;
    static integer nac1, ndv1, ndv2;
    static doublereal ctam, ctbm;
    extern /* Subroutine */ int cnmn08_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *);
    static doublereal thmax;

    /* Fortran I/O blocks */
    static cilist io___212 = { 0, 6, 0, fmt_43, 0 };
    static cilist io___217 = { 0, 6, 0, fmt_180, 0 };
    static cilist io___220 = { 0, 6, 0, fmt_178, 0 };


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*      revision history */
/*     double precision version for workstations */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     ROUTINE TO SOLVE DIRECTION FINDING PROBLEM IN MODIFIED METHOD OF */
/*     FEASIBLE DIRECTIONS. */
/*     BY G. N. VANDERPLAATS                            MAY, 1972. */
/*     NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF. */
/*     NORM OF S VECTOR USED HERE IS S-TRANSPOSE TIMES S.LE.1. */
/*     IF NVC = 0 FIND DIRECTION BY ZOUTENDIJK'S METHOD.  OTHERWISE */
/*     FIND MODIFIED DIRECTION. */
/*     ------------------------------------------------------------------ */
/*     ***  NORMALIZE GRADIENTS, CALCULATE THETA'S AND DETERMINE NVC  *** */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    --s;
    --df;
    --isc;
    --g;
    --ic;
    b_dim1 = *n3;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    a_dim1 = *n1;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --c__;
    --ms1;

    /* Function Body */
    ndv1 = cnmn1_1.ndv + 1;
    ndv2 = cnmn1_1.ndv + 2;
    nac1 = cnmn1_1.nac + 1;
    *nvc = 0;
    thmax = 0.f;
    cta = abs(cnmn1_1.ct);
    ct1 = 1.f / cta;
    ctam = abs(cnmn1_1.ctmin);
    ctb = abs(cnmn1_1.ctl);
    ct2 = 1.f / ctb;
    ctbm = abs(cnmn1_1.ctlmin);
    a1 = 1.f;
    i__1 = cnmn1_1.nac;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*     CALCULATE THETA */
	nci = ic[i__];
	ncj = 1;
	if (nci <= cnmn1_1.ncon) {
	    ncj = isc[nci];
	}
	c1 = g[nci];
	ctd = ct1;
	ctc = ctam;
	if (ncj <= 0) {
	    goto L10;
	}
	ctc = ctbm;
	ctd = ct2;
L10:
	if (c1 > ctc) {
	    ++(*nvc);
	}
	tht = 0.f;
	gg = ctd * c1 + 1.f;
	if (ncj == 0 || c1 > ctc) {
	    tht = cnmn1_1.theta * gg * gg;
	}
	if (tht > 50.f) {
	    tht = 50.f;
	}
	if (tht > thmax) {
	    thmax = tht;
	}
	a[ndv1 + i__ * a_dim1] = tht;
/*     ------------------------------------------------------------------ */
/*                    NORMALIZE GRADIENTS OF CONSTRAINTS */
/*     ------------------------------------------------------------------ */
	a[ndv2 + i__ * a_dim1] = 1.f;
	if (nci > cnmn1_1.ncon) {
	    goto L40;
	}
	a1 = 0.f;
	i__2 = cnmn1_1.ndv;
	for (j = 1; j <= i__2; ++j) {
/* Computing 2nd power */
	    d__1 = a[j + i__ * a_dim1];
	    a1 += d__1 * d__1;
/* L20: */
	}
	if (a1 < 1e-20) {
	    a1 = 1e-20;
	}
	a1 = sqrt(a1);
	a[ndv2 + i__ * a_dim1] = a1;
	a1 = 1.f / a1;
	i__2 = cnmn1_1.ndv;
	for (j = 1; j <= i__2; ++j) {
/* L30: */
	    a[j + i__ * a_dim1] = a1 * a[j + i__ * a_dim1];
	}
L40:
	;
    }
/*     ------------------------------------------------------------------ */
/*     CHECK FOR ZERO GRADIENT.  PROGRAM CHANGE-FEB, 1981, GV. */
/*     ------------------------------------------------------------------ */
    i__ = 0;
L41:
    ++i__;
L42:
    if (a[ndv2 + i__ * a_dim1] > 1e-6) {
	goto L45;
    }
/*     ZERO GRADIENT IS FOUND.  WRITE ERROR MESSAGE. */
    if (cnmn1_1.iprint >= 2) {
	s_wsfe(&io___212);
	do_fio(&c__1, (char *)&ic[i__], (ftnlen)sizeof(integer));
	e_wsfe();
    }
/*     REDUCE NAC BY ONE. */
    --cnmn1_1.nac;
/*     SHIFT COLUMNS OF A AND ROWS OF IC IF I.LE.NAC. */
    if (i__ > cnmn1_1.nac) {
	goto L46;
    }
/*     SHIFT. */
    i__1 = cnmn1_1.nac;
    for (j = i__; j <= i__1; ++j) {
	j1 = j + 1;
	ic[j] = ic[j1];
	i__2 = ndv2;
	for (k = 1; k <= i__2; ++k) {
/* L44: */
	    a[k + j * a_dim1] = a[k + j1 * a_dim1];
	}
    }
    if (i__ <= cnmn1_1.nac) {
	goto L42;
    }
L45:
    if (i__ < cnmn1_1.nac) {
	goto L41;
    }
L46:
    if (cnmn1_1.nac <= 0) {
	return 0;
    }
    nac1 = cnmn1_1.nac + 1;
/*     DETERMINE IF CONSTRAINTS ARE VIOLATED. */
    *nvc = 0;
    i__2 = cnmn1_1.nac;
    for (i__ = 1; i__ <= i__2; ++i__) {
	nci = ic[i__];
	ncj = 1;
	if (nci <= cnmn1_1.ncon) {
	    ncj = isc[nci];
	}
	ctc = ctam;
	if (ncj > 0) {
	    ctc = ctbm;
	}
	if (g[nci] > ctc) {
	    ++(*nvc);
	}
/* L47: */
    }
/*     ------------------------------------------------------------------ */
/*     NORMALIZE GRADIENT OF OBJECTIVE FUNCTION AND STORE IN NAC+1 */
/*     COLUMN OF A */
/*     ------------------------------------------------------------------ */
    a1 = 0.f;
    i__2 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing 2nd power */
	d__1 = df[i__];
	a1 += d__1 * d__1;
/* L50: */
    }
    if (a1 < 1e-20) {
	a1 = 1e-20;
    }
    a1 = sqrt(a1);
    a1 = 1.f / a1;
    i__2 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L60: */
	a[i__ + nac1 * a_dim1] = a1 * df[i__];
    }
/*     BUILD C VECTOR. */
    if (*nvc > 0) {
	goto L80;
    }
/*     ------------------------------------------------------------------ */
/*                 BUILD C FOR CLASSICAL METHOD */
/*     ------------------------------------------------------------------ */
    ndb = nac1;
    a[ndv1 + ndb * a_dim1] = 1.f;
    i__2 = ndb;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L70: */
	c__[i__] = -a[ndv1 + i__ * a_dim1];
    }
    goto L110;
L80:
/*     ------------------------------------------------------------------ */
/*                   BUILD C FOR MODIFIED METHOD */
/*     ------------------------------------------------------------------ */
    ndb = cnmn1_1.nac;
    a[ndv1 + nac1 * a_dim1] = -(*phi);
/*     ------------------------------------------------------------------ */
/*           SCALE THETA'S SO THAT MAXIMUM THETA IS UNITY */
/*     ------------------------------------------------------------------ */
    if (thmax > 1e-5f) {
	thmax = 1.f / thmax;
    }
    i__2 = ndb;
    for (i__ = 1; i__ <= i__2; ++i__) {
	a[ndv1 + i__ * a_dim1] *= thmax;
/* L90: */
    }
    i__2 = ndb;
    for (i__ = 1; i__ <= i__2; ++i__) {
	c__[i__] = 0.f;
	i__1 = ndv1;
	for (j = 1; j <= i__1; ++j) {
/* L100: */
	    c__[i__] += a[j + i__ * a_dim1] * a[j + nac1 * a_dim1];
	}
    }
L110:
/*     ------------------------------------------------------------------ */
/*                      BUILD B MATRIX */
/*     ------------------------------------------------------------------ */
    i__1 = ndb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ndb;
	for (j = 1; j <= i__2; ++j) {
	    b[i__ + j * b_dim1] = 0.f;
	    i__3 = ndv1;
	    for (k = 1; k <= i__3; ++k) {
/* L120: */
		b[i__ + j * b_dim1] -= a[k + i__ * a_dim1] * a[k + j * a_dim1]
			;
	    }
	}
    }
/*     ------------------------------------------------------------------ */
/*                    SOLVE SPECIAL L. P. PROBLEM */
/*     ------------------------------------------------------------------ */
    cnmn08_(&ndb, &ner, &c__[1], &ms1[1], &b[b_offset], n3, n4, n5);
    if (cnmn1_1.iprint > 1 && ner > 0) {
	s_wsfe(&io___217);
	e_wsfe();
    }
/*     CALCULATE RESULTING DIRECTION VECTOR, S. */
    *slope = 0.f;
/*     ------------------------------------------------------------------ */
/*                  USABLE-FEASIBLE DIRECTION */
/*     ------------------------------------------------------------------ */
    i__3 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__3; ++i__) {
	s1 = 0.f;
	if (*nvc > 0) {
	    s1 = -a[i__ + nac1 * a_dim1];
	}
	i__2 = ndb;
	for (j = 1; j <= i__2; ++j) {
/* L130: */
	    s1 -= a[i__ + j * a_dim1] * c__[j];
	}
	*slope += s1 * df[i__];
/* L140: */
	s[i__] = s1;
    }
    s[ndv1] = 1.f;
    if (*nvc > 0) {
	s[ndv1] = -a[ndv1 + nac1 * a_dim1];
    }
    i__3 = ndb;
    for (j = 1; j <= i__3; ++j) {
/* L150: */
	s[ndv1] -= a[ndv1 + j * a_dim1] * c__[j];
    }
/*     ------------------------------------------------------------------ */
/*     CHECK TO INSURE THE S-VECTOR IS FEASIBLE. */
/*     PROGRAM MOD-FEB, 1981, GV. */
/*     ------------------------------------------------------------------ */
    i__3 = cnmn1_1.nac;
    for (j = 1; j <= i__3; ++j) {
/*     S DOT DEL(G). */
	sg = 0.f;
	i__2 = cnmn1_1.ndv;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L172: */
	    sg += s[i__] * a[i__ + j * a_dim1];
	}
/*     IF(SG.GT.0.) GO TO 176 */

/*  THIS CHANGE MADE ON 4/8/81 FOR G. VANDERPLAATS */

	if (sg > 1e-4) {
	    goto L176;
	}
/*     FEASIBLE FOR THIS CONSTRAINT.  CONTINUE. */
/* L174: */
    }
    goto L179;
L176:
/*     S-VECTOR IS NOT FEASIBLE DUE TO SOME NUMERICAL PROBLEM. */
    if (cnmn1_1.iprint >= 2) {
	s_wsfe(&io___220);
	e_wsfe();
    }
    s[ndv1] = 0.f;
    *nvc = 0;
    return 0;
L179:
/*     ------------------------------------------------------------------ */
/*                  NORMALIZE S TO MAX ABS OF UNITY */
/*     ------------------------------------------------------------------ */
    s1 = 0.f;
    i__3 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__3; ++i__) {
	a1 = (d__1 = s[i__], abs(d__1));
	if (a1 > s1) {
	    s1 = a1;
	}
/* L160: */
    }
/*     IF (S1.LT.1.0E-10) RETURN */

/*  E-10 CHANGED TO E-04 ON 1/12/81 */

    if (s1 < 1e-4) {
	return 0;
    }
    s1 = 1.f / s1;
    i__3 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L170: */
	s[i__] = s1 * s[i__];
    }
    *slope = s1 * *slope;
    s[ndv1] = s1 * s[ndv1];
    return 0;
/*     ------------------------------------------------------------------ */
/*                           FORMATS */
/*     ------------------------------------------------------------------ */


} /* cnmn05_ */

/* ----- CNMN06 */
/* Subroutine */ int cnmn06_(doublereal *x, doublereal *vlb, doublereal *vub, 
	doublereal *g, doublereal *scal, doublereal *df, doublereal *s, 
	doublereal *g1, doublereal *g2, doublereal *ctam, doublereal *ctbm, 
	doublereal *slope, doublereal *alp, doublereal *a2, doublereal *a3, 
	doublereal *a4, doublereal *f1, doublereal *f2, doublereal *f3, 
	doublereal *cv1, doublereal *cv2, doublereal *cv3, doublereal *cv4, 
	doublereal *alpca, doublereal *alpfes, doublereal *alpln, doublereal *
	alpmin, doublereal *alpnc, doublereal *alpsav, doublereal *alpsid, 
	doublereal *alptot, integer *isc, integer *n1, integer *n2, integer *
	ncal, integer *nvc, integer *icount, integer *igood1, integer *igood2,
	 integer *igood3, integer *igood4, integer *ibest, integer *iii, 
	integer *nlnc, integer *jgoto)
{
    /* Format strings */
    static char fmt_730[] = "(/////\002* * * CONSTRAINED ONE-DIMENSIONAL SEA"
	    "RCH INFORMATION * * *\002)";
    static char fmt_740[] = "(//5x,\002PROPOSED DESIGN\002/5x,\002ALPHA ="
	    "\002,e12.5/5x,\002X-VECTOR\002)";
    static char fmt_750[] = "(1x,8e12.4)";
    static char fmt_760[] = "(/5x,\002OBJ =\002,e13.5)";
    static char fmt_770[] = "(/5x,\002CONSTRAINT VALUES\002)";
    static char fmt_780[] = "(/5x,\002TWO-POINT INTERPOLATION\002)";
    static char fmt_720[] = "(/5x,\002THREE-POINT INTERPOLATION\002)";
    static char fmt_790[] = "(/5x,\002* * * END OF ONE-DIMENSIONAL SEARCH"
	    "\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__;
    static doublereal c1, c2, c3, f4, cc, gi;
    static integer ii;
    static doublereal si, xi, xi1, xi2, zro;
    static integer nvc1;
    static doublereal alpa, alpb;
    static integer ksid;
    extern /* Subroutine */ int cnmn04_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), cnmn07_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer jbest;

    /* Fortran I/O blocks */
    static cilist io___222 = { 0, 6, 0, fmt_730, 0 };
    static cilist io___231 = { 0, 6, 0, fmt_740, 0 };
    static cilist io___232 = { 0, 6, 0, fmt_750, 0 };
    static cilist io___233 = { 0, 6, 0, fmt_750, 0 };
    static cilist io___234 = { 0, 6, 0, fmt_760, 0 };
    static cilist io___235 = { 0, 6, 0, fmt_770, 0 };
    static cilist io___236 = { 0, 6, 0, fmt_750, 0 };
    static cilist io___244 = { 0, 6, 0, fmt_780, 0 };
    static cilist io___245 = { 0, 6, 0, fmt_740, 0 };
    static cilist io___246 = { 0, 6, 0, fmt_750, 0 };
    static cilist io___247 = { 0, 6, 0, fmt_750, 0 };
    static cilist io___248 = { 0, 6, 0, fmt_760, 0 };
    static cilist io___249 = { 0, 6, 0, fmt_770, 0 };
    static cilist io___250 = { 0, 6, 0, fmt_750, 0 };
    static cilist io___252 = { 0, 6, 0, fmt_720, 0 };
    static cilist io___253 = { 0, 6, 0, fmt_740, 0 };
    static cilist io___254 = { 0, 6, 0, fmt_750, 0 };
    static cilist io___255 = { 0, 6, 0, fmt_750, 0 };
    static cilist io___257 = { 0, 6, 0, fmt_760, 0 };
    static cilist io___258 = { 0, 6, 0, fmt_770, 0 };
    static cilist io___259 = { 0, 6, 0, fmt_750, 0 };
    static cilist io___260 = { 0, 6, 0, fmt_790, 0 };


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*      revision history */
/*     double precision version for workstations */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     ROUTINE TO SOLVE ONE-DIMENSIONAL SEARCH PROBLEM FOR CONSTRAINED */
/*     FUNCTION MINIMIZATION. */
/*     BY G. N. VANDERPLAATS                           AUG., 1974. */
/*     NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF. */
/*     OBJ = INITIAL AND FINAL FUNCTION VALUE. */
/*     ALP = MOVE PARAMETER. */
/*     SLOPE = INITIAL SLOPE. */

/*     ALPSID = MOVE TO SIDE CONSTRAINT. */
/*     ALPFES = MOVE TO FEASIBLE REGION. */
/*     ALPNC = MOVE TO NEW NON-LINEAR CONSTRAINT. */
/*     ALPLN = MOVE TO LINEAR CONSTRAINT. */
/*     ALPCA = MOVE TO RE-ENCOUNTER CURRENTLY ACTIVE CONSTRAINT. */
/*     ALPMIN = MOVE TO MINIMIZE FUNCTION. */
/*     ALPTOT = TOTAL MOVE PARAMETER. */
    /* Parameter adjustments */
    --s;
    --df;
    --scal;
    --vub;
    --vlb;
    --x;
    --isc;
    --g2;
    --g1;
    --g;
    --ncal;

    /* Function Body */
    zro = 0.f;
    if (*jgoto == 0) {
	goto L10;
    }
    switch (*jgoto) {
	case 1:  goto L140;
	case 2:  goto L310;
	case 3:  goto L520;
    }
L10:
    if (cnmn1_1.iprint >= 5) {
	s_wsfe(&io___222);
	e_wsfe();
    }
    *alpsav = *alp;
    *icount = 0;
    *alptot = 0.f;
/*     TOLERANCES. */
    *ctam = abs(cnmn1_1.ctmin);
    *ctbm = abs(cnmn1_1.ctlmin);
/*     PROPOSED MOVE. */
L20:
/*     ------------------------------------------------------------------ */
/*     *****  BEGIN SEARCH OR IMPOSE SIDE CONSTRAINT MODIFICATION  ***** */
/*     ------------------------------------------------------------------ */
    *a2 = *alpsav;
    ++(*icount);
    *alpsid = 1e20;
/*     INITIAL ALPHA AND OBJ. */
    *alp = 0.f;
    *f1 = cnmn1_1.obj;
    ksid = 0;
    if (cnmn1_1.nside == 0) {
	goto L70;
    }
/*     ------------------------------------------------------------------ */
/*     FIND MOVE TO SIDE CONSTRAINT AND INSURE AGAINST VIOLATION OF */
/*     SIDE CONSTRAINTS */
/*     ------------------------------------------------------------------ */
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	si = s[i__];
	if (abs(si) > 1e-20) {
	    goto L30;
	}
/*     ITH COMPONENT OF S IS SMALL.  SET TO ZERO. */
	s[i__] = 0.f;
	*slope -= si * df[i__];
	goto L60;
L30:
	xi = x[i__];
	si = 1.f / si;
	if (si > 0.f) {
	    goto L40;
	}
/*     LOWER BOUND. */
	xi2 = vlb[i__];
	xi1 = abs(xi2);
	if (xi1 < 1.f) {
	    xi1 = 1.f;
	}
/*     CONSTRAINT VALUE. */
	gi = (xi2 - xi) / xi1;
	if (gi > -1e-6) {
	    goto L50;
	}
/*     PROPOSED MOVE TO LOWER BOUND. */
	alpa = (xi2 - xi) * si;
	if (alpa < *alpsid) {
	    *alpsid = alpa;
	}
	goto L60;
L40:
/*     UPPER BOUND. */
	xi2 = vub[i__];
	xi1 = abs(xi2);
	if (xi1 < 1.f) {
	    xi1 = 1.f;
	}
/*     CONSTRAINT VALUE. */
	gi = (xi - xi2) / xi1;
	if (gi > -1e-6) {
	    goto L50;
	}
/*     PROPOSED MOVE TO UPPER BOUND. */
	alpa = (xi2 - xi) * si;
	if (alpa < *alpsid) {
	    *alpsid = alpa;
	}
	goto L60;
L50:
/*     MOVE WILL VIOLATE SIDE CONSTRAINT.  SET S(I)=0. */
	*slope -= s[i__] * df[i__];
	s[i__] = 0.f;
	++ksid;
L60:
	;
    }
/*     ALPSID IS UPPER BOUND ON ALPHA. */
    if (*a2 > *alpsid) {
	*a2 = *alpsid;
    }
L70:
/*     ------------------------------------------------------------------ */
/*               CHECK ILL-CONDITIONING */
/*     ------------------------------------------------------------------ */
    if (ksid == cnmn1_1.ndv || *icount > 10) {
	goto L710;
    }
    if (*nvc == 0 && *slope > 0.f) {
	goto L710;
    }
    *alpfes = -1.f;
    *alpmin = -1.f;
    *alpln = *alpsid * 1.1f;
    *alpnc = *alpsid;
    *alpca = *alpsid;
    if (cnmn1_1.ncon == 0) {
	goto L90;
    }
/*     STORE CONSTRAINT VALUES IN G1. */
    i__1 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g1[i__] = g[i__];
/* L80: */
    }
L90:
/*     ------------------------------------------------------------------ */
/*                  MOVE A DISTANCE A2*S */
/*     ------------------------------------------------------------------ */
    *alptot += *a2;
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] += *a2 * s[i__];
/* L100: */
    }
    if (cnmn1_1.iprint < 5) {
	goto L130;
    }
    s_wsfe(&io___231);
    do_fio(&c__1, (char *)&(*a2), (ftnlen)sizeof(doublereal));
    e_wsfe();
    if (cnmn1_1.nscal == 0) {
	goto L120;
    }
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L110: */
	g[i__] = scal[i__] * x[i__];
    }
    s_wsfe(&io___232);
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&g[i__], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    goto L130;
L120:
    s_wsfe(&io___233);
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
/*     ------------------------------------------------------------------ */
/*                   UPDATE FUNCTION AND CONSTRAINT VALUES */
/*     ------------------------------------------------------------------ */
L130:
    ++ncal[1];
    *jgoto = 1;
    return 0;
L140:
    *f2 = cnmn1_1.obj;
    if (cnmn1_1.iprint >= 5) {
	s_wsfe(&io___234);
	do_fio(&c__1, (char *)&(*f2), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (cnmn1_1.iprint < 5 || cnmn1_1.ncon == 0) {
	goto L150;
    }
    s_wsfe(&io___235);
    e_wsfe();
    s_wsfe(&io___236);
    i__1 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&g[i__], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
L150:
/*     ------------------------------------------------------------------ */
/*               IDENTIFY ACCAPTABILITY OF DESIGNS F1 AND F2 */
/*     ------------------------------------------------------------------ */
/*     IGOOD = 0 IS ACCAPTABLE. */
/*     CV = MAXIMUM CONSTRAINT VIOLATION. */
    *igood1 = 0;
    *igood2 = 0;
    *cv1 = 0.f;
    *cv2 = 0.f;
    nvc1 = 0;
    if (cnmn1_1.ncon == 0) {
	goto L170;
    }
    i__1 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cc = *ctam;
	if (isc[i__] > 0) {
	    cc = *ctbm;
	}
	c1 = g1[i__] - cc;
	c2 = g[i__] - cc;
	if (c2 > 0.f) {
	    ++nvc1;
	}
	if (c1 > *cv1) {
	    *cv1 = c1;
	}
	if (c2 > *cv2) {
	    *cv2 = c2;
	}
/* L160: */
    }
    if (*cv1 > 0.f) {
	*igood1 = 1;
    }
    if (*cv2 > 0.f) {
	*igood2 = 1;
    }
L170:
    *alp = *a2;
    cnmn1_1.obj = *f2;
/*     ------------------------------------------------------------------ */
/*     IF F2 VIOLATES FEWER CONSTRAINTS THAN F1 BUT STILL HAS CONSTRAINT */
/*     VIOLATIONS RETURN */
/*     ------------------------------------------------------------------ */
    if (nvc1 < *nvc && nvc1 > 0) {
	goto L710;
    }
/*     ------------------------------------------------------------------ */
/*             IDENTIFY BEST OF DESIGNS F1 ANF F2 */
/*     ------------------------------------------------------------------ */
/*     IBEST CORRESPONDS TO MINIMUM VALUE DESIGN. */
/*     IF CONSTRAINTS ARE VIOLATED, IBEST CORRESPONDS TO MINIMUM */
/*     CONSTRAINT VIOLATION. */
    if (*igood1 == 0 && *igood2 == 0) {
	goto L180;
    }
/*     VIOLATED CONSTRAINTS.  PICK MINIMUM VIOLATION. */
    *ibest = 1;
    if (*cv1 >= *cv2) {
	*ibest = 2;
    }
    goto L190;
L180:
/*     NO CONSTRAINT VIOLATION.  PICK MINIMUM F. */
    *ibest = 1;
    if (*f2 <= *f1) {
	*ibest = 2;
    }
L190:
    ii = 1;
/*     ------------------------------------------------------------------ */
/*     IF CV2 IS GREATER THAN CV1, SET MOVE LIMITS TO A2. */
/*     PROGRAM MOD-FEB, 1981, GV. */
/*     ------------------------------------------------------------------ */
    if (*cv2 <= *cv1) {
	goto L195;
    }
    *alpln = *a2;
    *alpnc = *a2;
    *alpca = *a2;
L195:
    if (cnmn1_1.ncon == 0) {
	goto L230;
    }
/*     ------------------------------------------------------------------ */
/*     *****                 2 - POINT INTERPOLATION                ***** */
/*     ------------------------------------------------------------------ */
    *iii = 0;
L200:
    ++(*iii);
    c1 = g1[*iii];
    c2 = g[*iii];
    if (isc[*iii] == 0) {
	goto L210;
    }
/*     ------------------------------------------------------------------ */
/*                        LINEAR CONSTRAINT */
/*     ------------------------------------------------------------------ */
    if (c1 >= 1e-5 && c1 <= *ctbm) {
	goto L220;
    }
    cnmn07_(&ii, alp, &zro, &zro, &c1, a2, &c2, &zro, &zro);
    if (*alp <= 0.f) {
	goto L220;
    }
    if (c1 > *ctbm && *alp > *alpfes) {
	*alpfes = *alp;
    }
    if (c1 < cnmn1_1.ctl && *alp < *alpln) {
	*alpln = *alp;
    }
    goto L220;
L210:
/*     ------------------------------------------------------------------ */
/*                     NON-LINEAR CONSTRAINT */
/*     ------------------------------------------------------------------ */
    if (c1 >= 1e-5 && c1 <= *ctam) {
	goto L220;
    }
    cnmn07_(&ii, alp, &zro, &zro, &c1, a2, &c2, &zro, &zro);
    if (*alp <= 0.f) {
	goto L220;
    }
    if (c1 > *ctam && *alp > *alpfes) {
	*alpfes = *alp;
    }
    if (c1 < cnmn1_1.ct && *alp < *alpnc) {
	*alpnc = *alp;
    }
L220:
    if (*iii < cnmn1_1.ncon) {
	goto L200;
    }
L230:
    if (cnmn1_1.linobj > 0 || *slope >= 0.f) {
	goto L240;
    }
/*     CALCULATE ALPHA TO MINIMIZE FUNCTION. */
    cnmn04_(&ii, alpmin, &zro, &zro, f1, slope, a2, f2, &zro, &zro, &zro, &
	    zro);
L240:
/*     ------------------------------------------------------------------ */
/*                         PROPOSED MOVE */
/*     ------------------------------------------------------------------ */
/*     MOVE AT LEAST FAR ENOUGH TO OVERCOME CONSTRAINT VIOLATIONS. */
    *a3 = *alpfes;
/*     MOVE TO MINIMIZE FUNCTION. */
    if (*alpmin > *a3) {
	*a3 = *alpmin;
    }
/*     IF A3.LE.0, SET A3 = ALPSID. */
    if (*a3 <= 0.f) {
	*a3 = *alpsid;
    }
/*     LIMIT MOVE TO NEW CONSTRAINT ENCOUNTER. */
    if (*a3 > *alpnc) {
	*a3 = *alpnc;
    }
    if (*a3 > *alpln) {
	*a3 = *alpln;
    }
/*     MAKE A3 NON-ZERO. */
    if (*a3 <= 1e-20) {
	*a3 = 1e-20;
    }
/*     IF A3=A2=ALPSID AND F2 IS BEST, GO INVOKE SIDE CONSTRAINT */
/*     MODIFICATION. */
    alpb = 1.f - *a2 / *a3;
    alpa = 1.f - *alpsid / *a3;
    jbest = 0;
    if (abs(alpb) < 1e-10 && abs(alpa) < 1e-10) {
	jbest = 1;
    }
    if (jbest == 1 && *ibest == 2) {
	goto L20;
    }
/*     SIDE CONSTRAINT CHECK NOT SATISFIED. */
    if (cnmn1_1.ncon == 0) {
	goto L260;
    }
/*     STORE CONSTRAINT VALUES IN G2. */
    i__1 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g2[i__] = g[i__];
/* L250: */
    }
L260:
/*     IF A3=A2, SET A3=.9*A2. */
    if (abs(alpb) < 1e-10) {
	*a3 = *a2 * .9f;
    }
/*     MOVE AT LEAST .01*A2. */
    if (*a3 < *a2 * .01f) {
	*a3 = *a2 * .01f;
    }
/*     LIMIT MOVE TO 5.*A2. */
    if (*a3 > *a2 * 5.f) {
	*a3 = *a2 * 5.f;
    }
/*     LIMIT MOVE TO ALPSID. */
    if (*a3 > *alpsid) {
	*a3 = *alpsid;
    }
/*     MOVE A DISTANCE A3*S. */
    *alp = *a3 - *a2;
    *alptot += *alp;
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] += *alp * s[i__];
/* L270: */
    }
    if (cnmn1_1.iprint < 5) {
	goto L300;
    }
    s_wsfe(&io___244);
    e_wsfe();
    s_wsfe(&io___245);
    do_fio(&c__1, (char *)&(*a3), (ftnlen)sizeof(doublereal));
    e_wsfe();
    if (cnmn1_1.nscal == 0) {
	goto L290;
    }
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L280: */
	g[i__] = scal[i__] * x[i__];
    }
    s_wsfe(&io___246);
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&g[i__], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    goto L300;
L290:
    s_wsfe(&io___247);
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
L300:
/*     ------------------------------------------------------------------ */
/*              UPDATE FUNCTION AND CONSTRAINT VALUES */
/*     ------------------------------------------------------------------ */
    ++ncal[1];
    *jgoto = 2;
    return 0;
L310:
    *f3 = cnmn1_1.obj;
    if (cnmn1_1.iprint >= 5) {
	s_wsfe(&io___248);
	do_fio(&c__1, (char *)&(*f3), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (cnmn1_1.iprint < 5 || cnmn1_1.ncon == 0) {
	goto L320;
    }
    s_wsfe(&io___249);
    e_wsfe();
    s_wsfe(&io___250);
    i__1 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&g[i__], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
L320:
/*     ------------------------------------------------------------------ */
/*       CALCULATE MAXIMUM CONSTRAINT VIOLATION AND PICK BEST DESIGN */
/*     ------------------------------------------------------------------ */
    *cv3 = 0.f;
    *igood3 = 0;
    nvc1 = 0;
    if (cnmn1_1.ncon == 0) {
	goto L340;
    }
    i__1 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cc = *ctam;
	if (isc[i__] > 0) {
	    cc = *ctbm;
	}
	c1 = g[i__] - cc;
	if (c1 > *cv3) {
	    *cv3 = c1;
	}
	if (c1 > 0.f) {
	    ++nvc1;
	}
/* L330: */
    }
    if (*cv3 > 0.f) {
	*igood3 = 1;
    }
L340:
/*     DETERMINE BEST DESIGN. */
    if (*ibest == 2) {
	goto L360;
    }
/*     CHOOSE BETWEEN F1 AND F3. */
    if (*igood1 == 0 && *igood3 == 0) {
	goto L350;
    }
    if (*cv1 >= *cv3) {
	*ibest = 3;
    }
    goto L380;
L350:
    if (*f3 <= *f1) {
	*ibest = 3;
    }
    goto L380;
L360:
/*     CHOOSE BETWEEN F2 AND F3. */
    if (*igood2 == 0 && *igood3 == 0) {
	goto L370;
    }
    if (*cv2 >= *cv3) {
	*ibest = 3;
    }
    goto L380;
L370:
    if (*f3 <= *f2) {
	*ibest = 3;
    }
L380:
    *alp = *a3;
    cnmn1_1.obj = *f3;
/*     IF F3 VIOLATES FEWER CONSTRAINTS THAN F1 RETURN. */
    if (nvc1 < *nvc) {
	goto L710;
    }
/*     IF OBJECTIVE AND ALL CONSTRAINTS ARE LINEAR, RETURN. */
    if (cnmn1_1.linobj != 0 && *nlnc == cnmn1_1.ncon) {
	goto L710;
    }
/*     IF A3 = ALPLN AND F3 IS BOTH GOOD AND BEST RETURN. */
    alpb = 1.f - *alpln / *a3;
    if (abs(alpb) < 1e-20 && *ibest == 3 && *igood3 == 0) {
	goto L710;
    }
/*     IF A3 = ALPSID AND F3 IS BEST, GO INVOKE SIDE CONSTRAINT */
/*     MODIFICATION. */
    alpa = 1.f - *alpsid / *a3;
    if (abs(alpa) < 1e-20 && *ibest == 3) {
	goto L20;
    }
/*     ------------------------------------------------------------------ */
/*     **********            3 - POINT INTERPOLATION            ********* */
/*     ------------------------------------------------------------------ */
    *alpnc = *alpsid;
    *alpca = *alpsid;
    *alpfes = -1.f;
    *alpmin = -1.f;
/*     ------------------------------------------------------------------ */
/*     IF A3 IS GREATER THAN A2 AND CV3 IS GREATER THAN CV2, SET */
/*     MOVE LIMITS TO A3.  PROGRAM MOD-FEB, 1981, GV. */
/*     ------------------------------------------------------------------ */
    if (*a3 <= *a2 || *cv3 <= *cv2) {
	goto L285;
    }
    *alpln = *a3;
    *alpnc = *a3;
    *alpca = *a3;
L285:
    if (cnmn1_1.ncon == 0) {
	goto L440;
    }
    *iii = 0;
L390:
    ++(*iii);
    c1 = g1[*iii];
    c2 = g2[*iii];
    c3 = g[*iii];
    if (isc[*iii] == 0) {
	goto L400;
    }
/*     ------------------------------------------------------------------ */
/*     LINEAR CONSTRAINT.  FIND ALPFES ONLY.  ALPLN SAME AS BEFORE. */
/*     ------------------------------------------------------------------ */
    if (c1 <= *ctbm) {
	goto L430;
    }
    ii = 1;
    cnmn07_(&ii, alp, &zro, &zro, &c1, a3, &c3, &zro, &zro);
    if (*alp > *alpfes) {
	*alpfes = *alp;
    }
    goto L430;
L400:
/*     ------------------------------------------------------------------ */
/*                     NON-LINEAR CONSTRAINT */
/*     ------------------------------------------------------------------ */
    ii = 2;
    cnmn07_(&ii, alp, &zro, &zro, &c1, a2, &c2, a3, &c3);
    if (*alp <= zro) {
	goto L430;
    }
    if (c1 >= cnmn1_1.ct && c1 <= 0.f) {
	goto L410;
    }
    if (c1 > *ctam || c1 < 0.f) {
	goto L420;
    }
/*     ALP IS MINIMUM MOVE.  UPDATE FOR NEXT CONSTRAINT ENCOUNTER. */
L410:
    alpa = *alp;
    cnmn07_(&ii, alp, &alpa, &zro, &c1, a2, &c2, a3, &c3);
    if (*alp < *alpca && *alp >= alpa) {
	*alpca = *alp;
    }
    goto L430;
L420:
    if (*alp > *alpfes && c1 > *ctam) {
	*alpfes = *alp;
    }
    if (*alp < *alpnc && c1 < 0.f) {
	*alpnc = *alp;
    }
L430:
    if (*iii < cnmn1_1.ncon) {
	goto L390;
    }
L440:
    if (cnmn1_1.linobj > 0 || *slope > 0.f) {
	goto L450;
    }
/*     ------------------------------------------------------------------ */
/*              CALCULATE ALPHA TO MINIMIZE FUNCTION */
/*     ------------------------------------------------------------------ */
    ii = 3;
    if (*a2 > *a3 && (*igood2 == 0 && *ibest == 2)) {
	ii = 2;
    }
    cnmn04_(&ii, alpmin, &zro, &zro, f1, slope, a2, f2, a3, f3, &zro, &zro);
L450:
/*     ------------------------------------------------------------------ */
/*                       PROPOSED MOVE */
/*     ------------------------------------------------------------------ */
/*     MOVE AT LEAST ENOUGH TO OVERCOME CONSTRAINT VIOLATIONS. */
    *a4 = *alpfes;
/*     MOVE TO MINIMIZE FUNCTION. */
    if (*alpmin > *a4) {
	*a4 = *alpmin;
    }
/*     IF A4.LE.0, SET A4 = ALPSID. */
    if (*a4 <= 0.f) {
	*a4 = *alpsid;
    }
/*     LIMIT MOVE TO NEW CONSTRAINT ENCOUNTER. */
    if (*a4 > *alpln) {
	*a4 = *alpln;
    }
    if (*a4 > *alpnc) {
	*a4 = *alpnc;
    }
/*     LIMIT MOVE TO RE-ENCOUNTER CURRENTLY ACTIVE CONSTRAINT. */
    if (*a4 > *alpca) {
	*a4 = *alpca;
    }
/*     LIMIT A4 TO 5.*A3. */
    if (*a4 > *a3 * 5.f) {
	*a4 = *a3 * 5.f;
    }
/*     UPDATE DESIGN. */
    if (*ibest != 3 || cnmn1_1.ncon == 0) {
	goto L470;
    }
/*     STORE CONSTRAINT VALUES IN G2.  F3 IS BEST.  F2 IS NOT. */
    i__1 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g2[i__] = g[i__];
/* L460: */
    }
L470:
/*     IF A4=A3 AND IGOOD1=0 AND IGOOD3=1, SET A4=.9*A3. */
    *alp = *a4 - *a3;
    if (*igood1 == 0 && *igood3 == 1 && abs(*alp) < 1e-20) {
	*a4 = *a3 * .9f;
    }
/*     ------------------------------------------------------------------ */
/*                   MOVE A DISTANCE A4*S */
/*     ------------------------------------------------------------------ */
    *alp = *a4 - *a3;
    *alptot += *alp;
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] += *alp * s[i__];
/* L480: */
    }
    if (cnmn1_1.iprint < 5) {
	goto L510;
    }
    s_wsfe(&io___252);
    e_wsfe();
    s_wsfe(&io___253);
    do_fio(&c__1, (char *)&(*a4), (ftnlen)sizeof(doublereal));
    e_wsfe();
    if (cnmn1_1.nscal == 0) {
	goto L500;
    }
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L490: */
	g[i__] = scal[i__] * x[i__];
    }
    s_wsfe(&io___254);
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&g[i__], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    goto L510;
L500:
    s_wsfe(&io___255);
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
L510:
/*     ------------------------------------------------------------------ */
/*              UPDATE FUNCTION AND CONSTRAINT VALUES */
/*     ------------------------------------------------------------------ */
    ++ncal[1];
    *jgoto = 3;
    return 0;
L520:
    f4 = cnmn1_1.obj;
    if (cnmn1_1.iprint >= 5) {
	s_wsfe(&io___257);
	do_fio(&c__1, (char *)&f4, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (cnmn1_1.iprint < 5 || cnmn1_1.ncon == 0) {
	goto L530;
    }
    s_wsfe(&io___258);
    e_wsfe();
    s_wsfe(&io___259);
    i__1 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&g[i__], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
L530:
/*     DETERMINE ACCAPTABILITY OF F4. */
    *igood4 = 0;
    *cv4 = 0.f;
    if (cnmn1_1.ncon == 0) {
	goto L550;
    }
    i__1 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cc = *ctam;
	if (isc[i__] > 0) {
	    cc = *ctbm;
	}
	c1 = g[i__] - cc;
	if (c1 > *cv4) {
	    *cv4 = c1;
	}
/* L540: */
    }
    if (*cv4 > 0.f) {
	*igood4 = 1;
    }
L550:
    *alp = *a4;
    cnmn1_1.obj = f4;
/*     ------------------------------------------------------------------ */
/*                     DETERMINE BEST DESIGN */
/*     ------------------------------------------------------------------ */
    switch (*ibest) {
	case 1:  goto L560;
	case 2:  goto L610;
	case 3:  goto L660;
    }
L560:
/*     CHOOSE BETWEEN F1 AND F4. */
    if (*igood1 == 0 && *igood4 == 0) {
	goto L570;
    }
    if (*cv1 > *cv4) {
	goto L710;
    }
    goto L580;
L570:
    if (f4 <= *f1) {
	goto L710;
    }
L580:
/*     F1 IS BEST. */
    *alptot -= *a4;
    cnmn1_1.obj = *f1;
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] -= *a4 * s[i__];
/* L590: */
    }
    if (cnmn1_1.ncon == 0) {
	goto L710;
    }
    i__1 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] = g1[i__];
/* L600: */
    }
    goto L710;
L610:
/*     CHOOSE BETWEEN F2 AND F4. */
    if (*igood2 == 0 && *igood4 == 0) {
	goto L620;
    }
    if (*cv2 > *cv4) {
	goto L710;
    }
    goto L630;
L620:
    if (f4 <= *f2) {
	goto L710;
    }
L630:
/*     F2 IS BEST. */
    cnmn1_1.obj = *f2;
    *a2 = *a4 - *a2;
    *alptot -= *a2;
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] -= *a2 * s[i__];
/* L640: */
    }
    if (cnmn1_1.ncon == 0) {
	goto L710;
    }
    i__1 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] = g2[i__];
/* L650: */
    }
    goto L710;
L660:
/*     CHOOSE BETWEEN F3 AND F4. */
    if (*igood3 == 0 && *igood4 == 0) {
	goto L670;
    }
    if (*cv3 > *cv4) {
	goto L710;
    }
    goto L680;
L670:
    if (f4 <= *f3) {
	goto L710;
    }
L680:
/*     F3 IS BEST. */
    cnmn1_1.obj = *f3;
    *a3 = *a4 - *a3;
    *alptot -= *a3;
    i__1 = cnmn1_1.ndv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] -= *a3 * s[i__];
/* L690: */
    }
    if (cnmn1_1.ncon == 0) {
	goto L710;
    }
    i__1 = cnmn1_1.ncon;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] = g2[i__];
/* L700: */
    }
L710:
    *alp = *alptot;
    if (cnmn1_1.iprint >= 5) {
	s_wsfe(&io___260);
	e_wsfe();
    }
    *jgoto = 0;
    return 0;
/*     ------------------------------------------------------------------ */
/*                                  FORMATS */
/*     ------------------------------------------------------------------ */


} /* cnmn06_ */

/* ----- CNMN07 */
/* Subroutine */ int cnmn07_(integer *ii, doublereal *xbar, doublereal *eps, 
	doublereal *x1, doublereal *y1, doublereal *x2, doublereal *y2, 
	doublereal *x3, doublereal *y3)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal aa, bb, cc;
    static integer jj;
    static doublereal x21, x31, dy, x32, qq, yy, xb2, bac, xbar1;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*      revision history */
/*     double precision version for workstations */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     ROUTINE TO FIND FIRST XBAR.GE.EPS CORRESPONDING TO A REAL ZERO */
/*     OF A ONE-DIMENSIONAL FUNCTION BY POLYNOMIEL INTERPOLATION. */
/*     BY G. N. VANDERPLAATS                          APRIL, 1972. */
/*     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF. */
/*     II = CALCULATION CONTROL. */
/*          1:  2-POINT LINEAR INTERPOLATION, GIVEN X1, Y1, X2 AND Y2. */
/*          2:  3-POINT QUADRATIC INTERPOLATION, GIVEN X1, Y1, X2, Y2, */
/*              X3 AND Y3. */
/*     EPS MAY BE NEGATIVE. */
/*     IF REQUIRED ZERO ON Y DOES NOT EXITS, OR THE FUNCTION IS */
/*     ILL-CONDITIONED, XBAR = EPS-1.0 WILL BE RETURNED AS AN ERROR */
/*     INDICATOR. */
/*     IF DESIRED INTERPOLATION IS ILL-CONDITIONED, A LOWER ORDER */
/*     INTERPOLATION, CONSISTANT WITH INPUT DATA, WILL BE ATTEMPTED AND */
/*     II WILL BE CHANGED ACCORDINGLY. */
    xbar1 = *eps - 1.f;
    *xbar = xbar1;
    jj = 0;
    x21 = *x2 - *x1;
    if (abs(x21) < 1e-20) {
	return 0;
    }
    if (*ii == 2) {
	goto L30;
    }

L10:
/*     ------------------------------------------------------------------ */
/*                  II=1: 2-POINT LINEAR INTERPOLATION */
/*     ------------------------------------------------------------------ */
    *ii = 1;
    yy = *y1 * *y2;
    if (jj == 0 || yy < 0.f) {
	goto L20;
    }
/*     INTERPOLATE BETWEEN X2 AND X3. */
    dy = *y3 - *y2;
    if (abs(dy) < 1e-20) {
	goto L20;
    }
    *xbar = *x2 + *y2 * (*x2 - *x3) / dy;
    if (*xbar < *eps) {
	*xbar = xbar1;
    }
    return 0;
L20:
    dy = *y2 - *y1;
/*     INTERPOLATE BETWEEN X1 AND X2. */
    if (abs(dy) < 1e-20) {
	return 0;
    }
    *xbar = *x1 + *y1 * (*x1 - *x2) / dy;
    if (*xbar < *eps) {
	*xbar = xbar1;
    }
    return 0;
L30:
/*     ------------------------------------------------------------------ */
/*                 II=2: 3-POINT QUADRATIC INTERPOLATION */
/*     ------------------------------------------------------------------ */
    jj = 1;
    x31 = *x3 - *x1;
    x32 = *x3 - *x2;
    qq = x21 * x31 * x32;
    if (abs(qq) < 1e-20) {
	return 0;
    }
    aa = (*y1 * x32 - *y2 * x31 + *y3 * x21) / qq;
    if (abs(aa) < 1e-20) {
	goto L10;
    }
    bb = (*y2 - *y1) / x21 - aa * (*x1 + *x2);
    cc = *y1 - *x1 * (aa * *x1 + bb);
    bac = bb * bb - aa * 4.f * cc;
    if (bac < 0.f) {
	goto L10;
    }
    bac = sqrt(bac);
    aa = .5f / aa;
    *xbar = aa * (bac - bb);
    xb2 = -aa * (bac + bb);
    if (*xbar < *eps) {
	*xbar = xb2;
    }
    if (xb2 < *xbar && xb2 > *eps) {
	*xbar = xb2;
    }
    if (*xbar < *eps) {
	*xbar = xbar1;
    }
    return 0;
} /* cnmn07_ */

/* ----- CNMN08 */
/* Subroutine */ int cnmn08_(integer *ndb, integer *ner, doublereal *c__, 
	integer *ms1, doublereal *b, integer *n3, integer *n4, integer *n5)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal c1;
    static integer m2;
    static doublereal bb, cb, bi;
    static integer jj, kk;
    static doublereal bb1, eps;
    static integer ichk, nmax, iter1;
    static doublereal cbmin, cbmax;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*      revision history */
/*     double precision version for workstations */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     ROUTINE TO SOLVE SPECIAL LINEAR PROBLEM FOR IMPOSING S-TRANSPOSE */
/*     TIMES S.LE.1 BOUNDS IN THE MODIFIED METHOD OF FEASIBLE DIRECTIONS. */
/*     BY G. N. VANDERPLAATS                             APRIL, 1972. */
/*     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF. */
/*     REF.  'STRUCTURAL OPTIMIZATION BY METHODS OF FEASIBLE DIRECTIONS', */
/*     G. N. VANDERPLAATS AND F. MOSES, JOURNAL OF COMPUTERS */
/*     AND STRUCTURES, VOL 3, PP 739-755, 1973. */
/*     FORM OF L. P. IS BX=C WHERE 1ST NDB COMPONENTS OF X CONTAIN VECTOR */
/*     U AND LAST NDB COMPONENTS CONTAIN VECTOR V.  CONSTRAINTS ARE */
/*     U.GE.0, V.GE.0, AND U-TRANSPOSE TIMES V = 0. */
/*     NER = ERROR FLAG.  IF NER.NE.0 ON RETURN, PROCESS HAS NOT */
/*     CONVERGED IN 5*NDB ITERATIONS. */
/*     VECTOR MS1 IDENTIFIES THE SET OF BASIC VARIABLES. */
/*     ------------------------------------------------------------------ */
/*     CHOOSE INITIAL BASIC VARIABLES AS V, AND INITIALIZE VECTOR MS1 */
/*     ------------------------------------------------------------------ */
    /* Parameter adjustments */
    b_dim1 = *n3;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --c__;
    --ms1;

    /* Function Body */
    *ner = 1;
    m2 = *ndb << 1;
/*     CALCULATE CBMIN AND EPS AND INITIALIZE MS1. */
    eps = -1e10;
    cbmin = 0.f;
    i__1 = *ndb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	bi = b[i__ + i__ * b_dim1];
	cbmax = 0.f;
	if (bi < -1e-6) {
	    cbmax = c__[i__] / bi;
	}
	if (bi > eps) {
	    eps = bi;
	}
	if (cbmax > cbmin) {
	    cbmin = cbmax;
	}
/* L10: */
	ms1[i__] = 0;
    }
    eps *= 1e-4f;
/*     IF (EPS.LT.-1.0E-10) EPS=-1.0E-10 */

/*  E-10 CHANGED TO E-03 ON 1/12/81 */

    if (eps < -.001) {
	eps = -.001;
    }
    if (eps > -1e-4f) {
	eps = -1e-4f;
    }
    cbmin *= 1e-6;
/*     IF (CBMIN.LT.1.0D-10) CBMIN=1.0D-10 */

/*  E-10 CHANGED TO E-05 ON 1/12/81 */

    if (cbmin < 1e-5) {
	cbmin = 1e-5;
    }
    iter1 = 0;
    nmax = *ndb * 5;
/*     ------------------------------------------------------------------ */
/*     **********             BEGIN NEW ITERATION              ********** */
/*     ------------------------------------------------------------------ */
L20:
    ++iter1;
    if (iter1 > nmax) {
	return 0;
    }
/*     FIND MAX. C(I)/B(I,I) FOR I=1,NDB. */
    cbmax = cbmin * .9f;
    ichk = 0;
    i__1 = *ndb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c1 = c__[i__];
	bi = b[i__ + i__ * b_dim1];
/*     IF (BI.GT.EPS.OR.C1.GT.0.) GO TO 30 */
	if (bi > eps || c1 > -1e-5) {
	    goto L30;
	}

/*  0. CHANGED TO -1.0E-05 ON 1/12/81 */

	cb = c1 / bi;
	if (cb <= cbmax) {
	    goto L30;
	}
	ichk = i__;
	cbmax = cb;
L30:
	;
    }
    if (cbmax < cbmin) {
	goto L70;
    }
    if (ichk == 0) {
	goto L70;
    }
/*     UPDATE VECTOR MS1. */
    jj = ichk;
    if (ms1[jj] == 0) {
	jj = ichk + *ndb;
    }
    kk = jj + *ndb;
    if (kk > m2) {
	kk = jj - *ndb;
    }
    ms1[kk] = ichk;
    ms1[jj] = 0;
/*     ------------------------------------------------------------------ */
/*                     PIVOT OF B(ICHK,ICHK) */
/*     ------------------------------------------------------------------ */
    bb = 1.f / b[ichk + ichk * b_dim1];
    i__1 = *ndb;
    for (j = 1; j <= i__1; ++j) {
/* L40: */
	b[ichk + j * b_dim1] = bb * b[ichk + j * b_dim1];
    }
    c__[ichk] = cbmax;
    b[ichk + ichk * b_dim1] = bb;
/*     ELIMINATE COEFICIENTS ON VARIABLE ENTERING BASIS AND STORE */
/*     COEFICIENTS ON VARIABLE LEAVING BASIS IN THEIR PLACE. */
    i__1 = *ndb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ == ichk) {
	    goto L60;
	}
	bb1 = b[i__ + ichk * b_dim1];
	b[i__ + ichk * b_dim1] = 0.f;
	i__2 = *ndb;
	for (j = 1; j <= i__2; ++j) {
/* L50: */
	    b[i__ + j * b_dim1] -= bb1 * b[ichk + j * b_dim1];
	}
	c__[i__] -= bb1 * cbmax;
L60:
	;
    }
    goto L20;
L70:
    *ner = 0;
/*     ------------------------------------------------------------------ */
/*     STORE ONLY COMPONENTS OF U-VECTOR IN 'C'.  USE B(I,1) FOR */
/*     TEMPORARY STORAGE */
/*     ------------------------------------------------------------------ */
    i__1 = *ndb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b[i__ + b_dim1] = c__[i__];
/* L80: */
    }
    i__1 = *ndb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c__[i__] = 0.f;
	j = ms1[i__];
	if (j > 0) {
	    c__[i__] = b[j + b_dim1];
	}
	if (c__[i__] < 0.f) {
	    c__[i__] = 0.f;
	}
/* L90: */
    }
    return 0;
} /* cnmn08_ */

