extern "C" {
#include <stack-c.h>
#include <sciprint.h>
#include <Scierror.h>
}

#include <stack-c.h>
#include <api_scilab.h>

#include <call_function.hpp>

#include <math.h>
#include <stdio.h>
#include <string.h>

#include <vector>

using namespace std;

#define OK   1
#define FAIL 0

//#define DEBUG 1

#define _(A) (A)

typedef void (*call_f_objf)(double *, double *, int, void *);

extern voidf AddFunctionInTable (char *name, int *rep, FTAB *table);  

//
// General functions for scilab interface 
//

// SearchInDynLinks: a scilab function which tries to find a C/C++/Fortran function loaded via link
// SearchComp:       a function which tries to find a C/C++/Fortran function stored in the FTAB table (ours is empty)
// SetFunction:      this function looks inside the table or inside Scilab for dynamically loaded functions
// sciobj:           a wrapper function which has the same prototype as a dynamically loaded C/C++/Fortran function and which is used
//                   to call a Scilab script function like a C/C++/Fortran function
// Emptyfunc:        an empty function which is used when no dynamically functions has been found
// GetFunctionPtr:   the main function. Get the parameter on the stack. If it's a Scilab function, then it returns a pointer to
//                   sciobj. Otherwise, returns a pointer to the dynamically loaded function


extern "C" int SearchInDynLinks(char *op, void (**realop) ());
static int     SearchComp(FTAB *Ftab, char *op, void (**realop) ( ));  
static voidf   SetFunction(char *name, int *rep, FTAB *table);  
int            sciobj (double * x, double * f, int n_size_x, void * param);
static void    Emptyfunc(void) {} ;

// This function returns a pointer to a function or a Scilab pointer to a function
// The choice is performed on the type of the given Scilab parameter
voidf GetFunctionPtr(char * name, int n, FTAB Table[], voidf scifun, int * ifunc, int * lhs, int * rhs, int * is_list, vector<param_fobj> & param_list) 
{
  int type, rep, * var_addr = NULL, var_type, var_list_type;
  int nb_elems, * list_item_addr = NULL, i;
  struct param_fobj parameter_tmp;
  voidf f;
  char * func_name = NULL;
  SciErr _SciErr;

  _SciErr = getVarAddressFromPosition(pvApiCtx, n, &var_addr);
  _SciErr = getVarType(pvApiCtx, var_addr, &var_type);

  switch(var_type) 
    {
      // The function is defined as a string
    case sci_strings:
      getAllocatedSingleString(pvApiCtx, var_addr, &func_name);
      f = SetFunction(func_name, &rep, Table);

      if (rep==1)
	{
	  Scierror(999,"SetFunction problem\n");
	  f = (voidf)0;
	}

      freeAllocatedSingleString(func_name);
      break;

      // The function is a Scilab macro
    case sci_c_function: 
      GetRhsVar(n, EXTERNAL_DATATYPE, lhs, rhs, ifunc);
      f = (voidf)scifun ;
      break;

      // The function is a list
    case sci_list: 
      if (ElementType(n, 1)!=sci_c_function)
	{
	  Scierror(999,"First element of the list must be a Scilab compiled function\n");
	  f = (voidf)0;
	  break;
	}
      GetListRhsVar(n, 1, EXTERNAL_DATATYPE, lhs, rhs, ifunc);
      f = (voidf)scifun ;

      // Now get the parameters      
      _SciErr = getListItemNumber(pvApiCtx, var_addr, &nb_elems);
      for(i=1;i<nb_elems; i++)
	{
	  _SciErr = getListItemAddress(pvApiCtx, var_addr, i+1, &list_item_addr);
	  _SciErr = getVarType(pvApiCtx, list_item_addr, &var_list_type);
	  parameter_tmp.param_addr = list_item_addr;
	  parameter_tmp.param_type = var_list_type;
	  param_list.push_back(parameter_tmp);
	}
      *is_list = 1;
      break;
      
    default:
      Scierror(999,"Wrong parameter in %s ! (number %d)\n",name,n);
      f = (voidf)0;
    }

  return f;
}

// This function searches in the FTAB or in Scilab for corresponding function
voidf SetFunction(char * name, int * rep, FTAB table[]) 
{
  voidf loc;
  char * s = NULL;
  char buf[csiz];

  strncpy(buf,name,csiz);
  s = buf;
  while((*s!=' ')&&(*s != '\0')) {s++;};
  *s =  '\0';

  if (SearchComp(table,buf,&loc)==OK) 
    {
      *rep = 0;
      return(loc);
    }

  if (SearchInDynLinks(buf,&loc)>=0)
    {
      *rep = 0;
      return(loc);
    }

  loc = Emptyfunc;
  *rep = 1;

  sciprint(" Function %s not found\n",name);

  return(loc);
}


// This function search in FTAB (here, we will use FTab_call_f) for the corresponding name of the function
int SearchComp(FTAB Ftab[], char * op, void (**realop)()) 
{
  int i=0;

  while(Ftab[i].name!=(char *)0) 
    {
      int j;

#ifdef DEBUG
      printf("SearchComp: Ftab[%d] = %s\n",i, Ftab[i].name);
#endif
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
	      return(FAIL);
	    }
	  else i++;
	}
    }

  return(FAIL);
}

