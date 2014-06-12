#ifndef FUNC_UTILS_H
#define FUNC_UTILS_H

#define OK   1
#define FAIL 0

/******************************
 * Work around for the        *
 * AddFunctionInTable problem *
 ******************************/

typedef void (*voidf)();
typedef struct {
  char *name;
  voidf f;
} FTAB;

extern voidf AddFunctionInTable(char *name, int *rep, FTAB *table);  

/******************************************
 * General functions for scilab interface *
 ******************************************/

// SearchInDynLinks: a scilab function which tries to find a C/C++/Fortran function loaded via link
// SearchComp: a function which tries to find a C/C++/Fortran function stored in the FTAB table (ours is empty)
// SetFunction: this function looks inside the table or inside Scilab for dynamically loaded functions
// sciobj: a wrapper function which has the same prototype as a dynamically loaded C/C++/Fortran function and which is used
//         to call a Scilab script function like a C/C++/Fortran function
// Emptyfunc: an empty function which is used when no dynamically functions has been found
// GetFunctionPtr: the main function. Get the parameter on the stack. If it's a Scilab function, then it returns a pointer to
//                 sciobj. Otherwise, returns a pointer to the dynamically loaded function

extern int   SearchInDynLinks(char *op, void (**realop) ());
int   SearchComp(FTAB *Ftab, char *op, void (**realop) ( ));  
voidf SetFunction(char *name, int *rep, FTAB *table);  
void  Emptyfunc(void);
voidf GetFunctionPtr(char *, int, FTAB *, voidf, int *, int*, int*);

#endif
