//'NEWUOA Scilab Interface' (interface with scilab of the NEWUOA Library)
//The 'NEWUOA Library' is an unconstrainted optimization algorithm without derivtatives	
// written by MJD Powell
//
//Copyright (C) 2005 GUILBERT
//See the 'LICENCE' File

#include "stack-c.h"

#ifndef CPOWELL_H
#define CPOWELL_H

// From f2c.h

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (doublereal)min(a,b)
#define dmax(a,b) (doublereal)max(a,b)
#define bit_test(a,b)	((a) >> (b) & 1)
#define bit_clear(a,b)	((a) & ~((uinteger)1 << (b)))
#define bit_set(a,b)	((a) |  ((uinteger)1 << (b)))

/****************************************************************/
/*Function written by M.J.D. Powell and translated into C by f2c*/
/****************************************************************/

extern int newuoa(double *, double *, long int *, long int *, double *, 
	    double *, double *, long int *, long int *, 
		void (*)(long int *, double *, double *));

int newuob_(double *, long int *, long int *, double *, 
	    double *, double *, long int *, long int *, void (*) (long int *, double *, double *),
		double *, double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, long int *, double *, double *, double *);

int biglag_(long int *, long int *, double *, 
	    double *, double *, double *, long int *, long int *, 
	    long int *, double *, double *, double *, double *,
	     double *, double *, double *, double *);
	
int bigden_(long int *, long int *, double *, double *, double *, 
	    double *, long int *, long int *, long int *, long int *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *);

int update_(long int *, long int *, 
	    double *, double *, long int *, long int *, double *, 
	    double *, long int *, double *);

int trsapp_(long int *, long int *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *);

#endif
