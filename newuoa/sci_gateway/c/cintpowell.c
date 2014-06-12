//'NEWUOA Scilab Interface' (Interface with scilab of the NEWUOA Library)
// The 'NEWUOA Library' is an unconstrainted optimization algorithm without derivtatives	
// written by MJD Powell
//
//Copyright (C) 2005 GUILBERT
//See the 'LICENCE' File

#include "cpowell.h"
#include "cintpowell.h"

extern int C2F(error)(int *n);
extern int SearchInDynLinks(char *op, void (**realop) ());

/********************************************************
 * cintpowell : interface for newuoa                    *
 *              An optimization toolbox by MJD Powell   *
 *              Optimization without derivatives        *
 *              Interface For Scilab                    *
 *                        By Matthieu GUILBERT Staubli  *              
 ********************************************************/

int cintpowell(char * fname) 
{ 
  //int returned_from_longjump ;
  int un=1;
  int sW=0;
  long int n, npt, iprint, maxfun;
  double  rhobeg, rhoend;
  
  static int l1, m1, n1,l2, m2, n2,l3, m3, n3,l4, m4, n4,l5, m5, n5,
			 l6, m6, n6,l7, m7, n7, func,w;
	
  /*    Define minls=1, maxlhs, minrhs, maxrhs   */
  static int minlhs=1, minrhs=7, maxlhs=1, maxrhs=7;

  
  /*   Check rhs and lhs   */  
  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;
  

 /*extern int newuoa(long int *n, long int *npt, double *x, 
	double *rhobeg, double *rhoend, long int *iprint, long int *
	maxfun, void (*obj)(long int *, double *, double *))*/
  
  /* RhsVar: newuoa(npt,x0,rhobeg,rhoend,iprint,maxfun,func) */
  /*                 1, 2,   3,      4,     5,    6,     7   */

  
  /*   Variable 1 (n)  pt */
  GetRhsVar(1, "i", &m1, &n1, &l1);
  npt = *istk(l1);
  
  if (m1!=1 || n1!=1)
  {
	sciprint("Error : first argument must be a scalar");
	return 0;
  }

  /*   Variable 2 (x)   */
  GetRhsVar(2, "d", &m2, &n2, &l2);
  n = m2*n2;

 	if (n==0)
	{
		sciprint("Error : the size of x must at least be 1");
		return 0;
	}

  /*   Variable 3*/
  GetRhsVar(3, "d", &m3, &n3, &l3);
  rhobeg = *stk(l3);

	if (m3!=1 || n3!=1)
	{
		sciprint("Error : Third argument must be a scalar");
		return 0;
	}

  /*   Variable 4*/
  GetRhsVar(4, "d", &m4, &n4, &l4);
  rhoend = *stk(l4);
  if (m4!=1 || n4!=1)
  {
	sciprint("Error : Fourth argument must be a scalar");
	return 0;
  }

  /*   Variable 5*/
  GetRhsVar(5, "i", &m5, &n5, &l5);
  iprint = *istk(l5);

  if (m5!=1 || n5!=1 || iprint>3 || iprint < 0)
  {
	sciprint("Error : IPRINT argument must be a scalar and equal to 0 or 1 or 2 or 3");
	return 0;
  }


  /*   Variable 6*/
  GetRhsVar(6, "i", &m6, &n6, &l6);
  maxfun = *istk(l6);

  if (m6!=1 || n6!=1)
  {
	sciprint("Error : Sixth argument must be a scalar");
	return 0;
  }

  /*   Variable 7 (objective function)   */

  objective =  (cpowellobjf) GetFunctionPtr("cpowell",7,FTab_cpowellobj,(voidf) sciobj,&sci_obj,&lhs_obj,&rhs_obj);
   if ( objective == (cpowellobjf) 0 ) 
     {
       sciprint("Error : Last argument must be a pointer to a scilab function");
       return 0;
     }
	 
  /*  Internal variables of f!*/

  //CreateVar1(8, "d", &un, &un, &func);
  CreateVar(8, "d", &un, &un, &func);

  sW = ((npt+13)*(npt+n)+3*(n)*(n+3)/2+10);
  //CreateVar1(9, "d", &sW, &un, &w);
  CreateVar(9, "d", &sW, &un, &w);
 	
  newuoa(    stk(w),
	     stk(func),
	     &n, 
	     &npt, 
	     stk(l2), 
	     &rhobeg, 
	     &rhoend, 
	     &iprint, 
	     &maxfun, 
	     objective);

  /* LhsVar: x = newuoa (...) */

  LhsVar(1) = 2;
  return 0;
  
}

/*Function in order to call the scilab function*/

static void sciobj(n,x,f)
     long int *n;
     double *x,*f;    
{
  int scilab_n, scilab_x;
  int un=1;
  int nb = (int) *n;

  CreateVar1(10,"d",&un,&un, &scilab_n);  
  stk(scilab_n)[0] = (double) (*n);

  CreateVar1(11,"d",&nb,&un, &scilab_x);  
  C2F(dcopy)(&nb,x,&un,stk(scilab_x),&un);

  ExecSciFunction(10, &sci_obj, &lhs_obj, &rhs_obj,"sciobj");
  
  /* One output at position of first input (16) */
  /* scilab_fj=scilab_j */
  
  *f=stk(scilab_n)[0];
}

voidf GetFunctionPtr(name,n,Table,scifun,ifunc,lhs,rhs) 
	char *name;
	int n,*lhs,*rhs,*ifunc;
	FTAB *Table;
	voidf scifun;
{
  int type,rep,mm,nn;
  voidf f;
  type=VarType(n);
  switch ( type) 
    {
    case a_chain : 
      GetRhsVar(n, "c", &mm, &nn, ifunc);
      f = SetFunction(cstk(*ifunc),&rep,Table);
      if ( rep == 1 )
	{
	  Error(999);
	  return (voidf) 0;
	}
      return f ;
    case  a_function : 
      GetRhsVar(n, "f", lhs,rhs, ifunc);
      return (voidf) scifun ;
    default: 
      sciprint("Wrong parameter in %s ! (number %d)\r\n",name,n);
      Error(999);
      return (voidf) 0 ;
    }
}

/*******************************************
 * General function 
 *******************************************/

static voidf SetFunction(name,rep,table) 
     char *name;
     int *rep;
     FTAB table[];
{
  void (*loc)();
  char *s;
  strncpy(buf,name,MAXNAME);
  s=buf ; while ( *s != ' ' && *s != '\0') { s++;};
  *s= '\0';
  if ( SearchComp(table,buf,&loc) == OK) 
    {
      *rep = 0;
      return(loc);
    }
  if ( SearchInDynLinks(buf,&loc) >= 0 )
    {
      *rep = 0;
      return(loc);
    }
  loc = Emptyfunc;
  *rep = 1;
  sciprint(" Function %s not found\r\n",name);
  return(loc);
}


static int SearchComp(Ftab,op,realop) 
     FTAB Ftab[];
     void (**realop)();
     char *op;
{
  int i=0;
  while ( Ftab[i].name != (char *) 0) 
     {
       int j;
       j = strcmp(op,Ftab[i].name);
       if ( j == 0 )
         {
           *realop = Ftab[i].f;
	   return(OK);
	 }
       else
         { 
	   if ( j <= 0)
             {
               /* sciprint("\nUnknow function <%s>\r\n",op); */
               return(FAIL);
             }
	   else i++;
         }
     }
  /* sciprint("\n Unknow function <%s>\r\n",op); */
  return(FAIL);
}
