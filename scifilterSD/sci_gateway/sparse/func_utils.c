#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <setjmp.h>

#include <stack-c.h>
#include <sciprint.h>
#include <Scierror.h>
#include <sci_types.h>
#include <MALLOC.h>

#include <api_scilab.h>
#include <api_parameters.h>
#include <freeArrayOfString.h>

#define OK   1
#define FAIL 0

#define GLCPD_ERROR_RETURN_NULL if(_SciErr.iErr)	\
    {							\
      printError(&_SciErr, 0);				\
      return NULL;					\
    }

/******************************
 * Work around for the        *
 * AddFunctionInTable problem *
 ******************************/

typedef void (*voidf)();
typedef struct {
  char *name;
  voidf f;
} FTAB;

/*******************************************
 * Tools for linked functions manipulation *
 *******************************************/

void  Emptyfunc(void) {};

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

  return(loc);
}

// This function search in FTAB (here, we will use FTab_cobyla_function) for the corresponding name of the function
int SearchComp(FTAB Ftab[], char * op, void (**realop)()) 
{
  int i=0;

  while(Ftab[i].name!=(char *)0) 
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
	      return(FAIL);
	    }
	  else i++;
	}
    }

  return(FAIL);
}

// This function returns a pointer to a function or a Scilab pointer to a function
// The choice is performed on the type of the given Scilab parameter
voidf GetFunctionPtr(char * name, int n, FTAB Table[], voidf scifun, int * ifunc, int * lhs, int * rhs) 
{
  int type, rep, m_tmp, n_tmp, i;
  int * tmp_addr = NULL;
  int * pi_len = NULL;
  char ** pst_strings = NULL;
  voidf f;
  SciErr _SciErr;

  _SciErr = getVarAddressFromPosition(pvApiCtx, n, &tmp_addr); GLCPD_ERROR_RETURN_NULL;
  _SciErr = getVarType(pvApiCtx, tmp_addr, &type); GLCPD_ERROR_RETURN_NULL;

  switch(type) 
    {
    case sci_strings: 
      _SciErr = getMatrixOfString(pvApiCtx, tmp_addr, &n_tmp, &m_tmp, NULL, NULL); GLCPD_ERROR_RETURN_NULL;
      pi_len = (int *)MALLOC(n_tmp*m_tmp*sizeof(int));
      pst_strings = (char **)MALLOC(n_tmp*m_tmp*sizeof(char *));
      _SciErr = getMatrixOfString(pvApiCtx, tmp_addr, &n_tmp, &m_tmp, pi_len, NULL); GLCPD_ERROR_RETURN_NULL;
      for(i=0;i<n_tmp*m_tmp;i++) pst_strings[i] = (char *)MALLOC((pi_len[i]+1)*sizeof(char));
      _SciErr = getMatrixOfString(pvApiCtx, tmp_addr, &n_tmp, &m_tmp, pi_len, pst_strings); GLCPD_ERROR_RETURN_NULL;

      f = SetFunction(pst_strings[0], &rep, Table);

      if (pst_strings) freeArrayOfString(pst_strings, n_tmp*m_tmp);
      if (pi_len) FREE(pi_len);

      if (rep==1)
	{
	  Scierror(999,"Function not found is %s\n", name);
	  return (voidf)0;
	}

      return f ;

    case sci_c_function: 
      GetRhsVar(n, "f", lhs, rhs, ifunc);
      return (voidf)scifun ;
      
    default: 
      Scierror(999,"Wrong parameter in %s ! (number %d)\r\n",name,n);
      return (voidf)0;
    }
}
