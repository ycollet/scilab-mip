//'NEWUOA Scilab Interface' (interface with scilab of the NEWUOA Library)
//The 'NEWUOA Library' is an unconstrainted optimization algorithm without derivtatives	
// written by MJD Powell
//
//Copyright (C) 2005 GUILBERT
//See the 'LICENCE' File

#include "stack-c.h"
#include "machine.h"
#include "sciprint.h"
//#include "sun/link.h"
#include <setjmp.h>

#ifndef CINTPOWELL_H
#define CINTPOWELL_H

/******* Please do not edit ****************/

#define OK 1
#define FAIL 0
#define MAXNAME 32
#define a_chain 10
#define a_function 13

#if defined(__STDC__)
#define ARGS_cpowellobj long int *, double *, double *
typedef void (*cpowellobjf)(ARGS_cpowellobj);
#else
#define ARGS_cpowellobj 
typedef void (*cpowellobjf)();
#endif 


/** Dealing with scilab interface error**/

#define ExecSciFunction(n,mx,nx,lx,name) \
  if(! C2F(scifunction)((c_local=n,&c_local),mx,nx,lx))\
{ sciprint("Error in function %s\r\n",name);  longjmp(cpowellenv,-1); }


/** like CreateVar but with return for void functions */

#define CreateVar1(n,ct,mx,nx,lx) if(! C2F(createvar)((c_local=n,&c_local),ct,mx,nx,lx, 1L))\
        { return ;  }

typedef void (*voidf)();

typedef struct {
  char *name;
  voidf f;
} FTAB;

/***********************/
/*Definition of FTables*/
/***********************/
FTAB FTab_cpowellobj[] ={{(char *) 0, (voidf) 0}};

/***************************************************
 * data for interface 
 ***************************************************/

static char buf[MAXNAME];
static  jmp_buf cpowellenv; 
static int sci_obj, lhs_obj, rhs_obj;
cpowellobjf objective;


/*General functions for scilab interface*/
static voidf SetFunction(char *name, int *rep, FTAB *table);  
static int SearchComp(FTAB *Ftab, char *op, void (**realop) ( ));  
static void Emptyfunc(void) {} ;
extern int C2F(gettype)(int *);
extern int C2F(scifunction)(int *,int *,int *,int *);
extern int C2F(dcopy)(int *,double *,int *,double *,int *);
static void sciobj (ARGS_cpowellobj);
static voidf GetFunctionPtr(char *,int,FTAB *,voidf,int *,int*,int*);

#endif
