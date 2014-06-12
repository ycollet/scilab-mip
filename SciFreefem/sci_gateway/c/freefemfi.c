#include <string.h>

#include <api_scilab.h>
#include <MALLOC.h>
#include <Scierror.h>
#include <stack-c.h>

void C2F(delete_lb_lobj)();
void C2F(delete_lexp)();
void C2F(add_lobj)();
void C2F(delete_ptr)(void * ptr);

void C2F(delete_scilabana)(int * ierr);
void C2F(freefem_code)(char * Message, int * err);
void C2F(ff_problem)(char * Message, int * flag, int * err);
void C2F(get_ff_result)(char * s, int * lhs, 
			double ** tNode, int * mNode, int * nNode, 
			double ** tFunc, int * mFunc, int * nFunc, 
			int ** tTriangle, int * mTriangle, int * nTriangle,
			char * Message,int * err);
void C2F(put_scilab_mesh)(char * s,double * tNode, int * mNode, int * tTriangle, int * mTriangle, char * Message, int * err);
void C2F(put_scilab_border_1)(char * nameBorder,double * tBorder, int * nbData,int * normal, char * Message,int * err);
void C2F(put_scilab_border_2)(char * nameBorder,char * equation, double * start, double * stop,	int * ng,int * nrm, char * Message, int * err);
void C2F(build_scilab_mesh)(char * nameMesh, char * Message,int * err);
void C2F(get_matrix)(double ** tMatrix, int * sztMatrix, int ** tjlow,int ** tjhigh, int * sztVect, int ** szblock, char * Message, int * err);

int int_ff_exec(char * fname)
{
  int minrhs = 1, maxrhs = 1;
  int minlhs = 0, maxlhs = 1;
  int m_message, n_message, * message_addr = NULL;
  int err = 0;
  char * message = NULL;
  SciErr _sciErr;

  CheckRhs(minrhs, maxrhs);
  CheckLhs(minlhs, maxlhs);

  _sciErr = getVarAddressFromPosition(pvApiCtx, 1, &message_addr);
  _sciErr = getVarDimension(pvApiCtx, message_addr, &n_message, &m_message);

  if (n_message*m_message!=1)
    {
      Scierror(999,"%s: Error, a single string expected.\n", fname);
      return 0;
    }

  getAllocatedSingleString(pvApiCtx, message_addr, &message);

  C2F(freefem_code)(message,&err);

  if (err)
    {
      Scierror(999,"%s: problem when calling freefem_code: %s.\n",fname, message);
      freeAllocatedSingleString(message);
      return 0;
    }

  LhsVar(1) = 0;

  freeAllocatedSingleString(message);

  return 0;
}

int int_ff_problem(char * fname)
{
  int minrhs = 1, maxrhs = 2;
  int minlhs = 0, maxlhs = 1;
  int * message_addr = NULL, m_message, n_message;
  int * choice_addr = NULL, m_choice, n_choice;
  char * message = NULL;
  int defval = 0, err = 0, choice;
  double tmpdbl;
  SciErr _sciErr;

  CheckRhs(minrhs, maxrhs);
  CheckLhs(minlhs, maxrhs);

  _sciErr = getVarAddressFromPosition(pvApiCtx, 1, &message_addr);
  _sciErr = getVarDimension(pvApiCtx, message_addr, &n_message, &m_message);

  if (n_message*m_message!=1)
    {
      Scierror(999,"%s: Error, a single string expected.\n", fname);
      return 0;
    }

  getAllocatedSingleString(pvApiCtx, message_addr, &message);

  if (Rhs == 1)
    {
      defval = 0;
      C2F(ff_problem)(message,&defval,&err);
    }
  else 
    {
      if (Rhs == 2)
	{
	  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &choice_addr);
	  getScalarDouble(pvApiCtx, choice_addr, &tmpdbl);
	  choice = (int)tmpdbl;
	  
	  C2F(ff_problem)(message,&choice,&err);
	}
      else
	{
	  err = 1;
	}
    }

  if (err)
    {
      Scierror(999,"%s: problem when calling ff_problem - too many arguments.\n",fname);
      return 0;
    }

  LhsVar(1) = 0;

  return 0;
}

int int_ff_end(char * fname)
{
  int minrhs = 0, maxrhs = 0;
  int minlhs = 0, maxlhs = 1;
  int result = 0;

  CheckLhs(minlhs, maxlhs);
  CheckRhs(minrhs, maxrhs);

  C2F(delete_scilabana)(&result);

  if (result)
    {
      Scierror(999, "%s: Syntax analyser already destroyed \n", fname);
      return 0;
    }

  createScalarDouble(pvApiCtx, 1, (double)result);

  LhsVar(1) = 1;

  return 0;
}

int int_get_ff_result(char * fname)
{
  int mNode,     nNode;     double * lNode     = NULL;
  int mTriangle, nTriangle; int    * lTriangle = NULL;
  int mFunc,     nFunc;     double * lFunc     = NULL;
  char * str = NULL, * param = NULL, buffer[4096];
  int m_param, n_param, * param_addr = NULL, i = 0;
  int minrhs = 0, maxrhs = 1;
  int minlhs = 1, maxlhs = 3;
  int lhs = Lhs, err = 0;
  int * tmpIntPtr = NULL;
  double * lTriangleDbl = NULL;
  SciErr _sciErr;

  CheckRhs(minrhs, maxrhs);
  CheckLhs(minlhs, maxlhs);

  if (Lhs == 2)
    {
      if (Rhs == 1)
	{
	  Scierror(999, "%s: No input parameter needed in this configuration.\n", fname);
	  return 0;
	}
    }
  else
    {
      if (Rhs == 0)
	{
	  Scierror(999, "%s: You need to specify on which data we get the result.\n", fname);
	  return 0;
	}
      
      _sciErr = getVarAddressFromPosition(pvApiCtx, 1, &param_addr);
      _sciErr = getVarDimension(pvApiCtx, param_addr, &n_param, &m_param);

      if (n_param*m_param!=1)
	{
	  Scierror(999,"%s: Error, a single string expected.\n", fname);
	  return 0;
	}

      getAllocatedSingleString(pvApiCtx, param_addr, &param);
    }
  
  buffer[0] = ' ';

  if (Lhs != 2)
    {
      C2F(get_ff_result)(param, &lhs, 
			 &lNode, &mNode, &nNode,
			 &lFunc, &mFunc, &nFunc,
			 &lTriangle, &mTriangle, &nTriangle,
			 buffer, &err);
    }
  else
    {
      str = (char *)MALLOC(3*sizeof(char)); str[0] = '0';
      C2F(get_ff_result)(str, &lhs,
			 &lNode, &mNode, &nNode,
			 &lFunc, &mFunc, &nFunc,
			 &lTriangle, &mTriangle, &nTriangle,
			 buffer, &err);
      free(str);
    }
  
  freeAllocatedSingleString(param);

  if (err) 
    {
      Scierror(999, "%s: problem when calling get_ff_results.\n", fname);
      return 0;
    }
  
  if (Lhs != 1)
    {
      createMatrixOfDouble(pvApiCtx, Rhs + 1, mNode, nNode, lNode);
      C2F(delete_ptr)(lNode);
      lTriangleDbl = (double *)MALLOC((mTriangle*nTriangle)*sizeof(double));
      for(i = 0; i<mTriangle*nTriangle; i++) lTriangleDbl[i] = (double)lTriangle[i];
      
      _sciErr = createMatrixOfDouble(pvApiCtx, Rhs + 2, mTriangle, nTriangle, lTriangleDbl);
      free(lTriangle);
      FREE(lTriangleDbl);
    }
  
  if (Lhs != 2)
    {
      _sciErr = createMatrixOfDouble(pvApiCtx, Rhs + Lhs, mFunc, nFunc, lFunc);
      C2F(delete_ptr)(lFunc);
    }
  
  if (Lhs == 1)
    {
      LhsVar(1) = Rhs + 1;
    }
  else if (Lhs == 2)
    {
      LhsVar(1) = Rhs + 1;
      LhsVar(2) = Rhs + 2;
    }
  else if (Lhs == 3) 
    {
      LhsVar(1) = Rhs + 1;
      LhsVar(2) = Rhs + 2;
      LhsVar(3) = Rhs + 3;
    }
  
  return 0;
}

int int_put_mesh(char * fname)
{
  char bordername[20], * fct = NULL;
  int minrhs = 2, maxrhs = 2;
  int minlhs = 1, maxlhs = 1;
  int m_list, n_list, * list_addr = NULL, * item_list_addr = NULL, * item_item_list_addr = NULL;
  int m_name_mesh, n_name_mesh, * name_mesh_addr = NULL;
  int m_matrix_dbl, n_matrix_dbl, nb_points = 0;
  int m_string_list, n_string_list;
  int m_param, n_param;
  char * name_mesh = NULL, buffer[4096];
  char ** string_list = NULL;
  int nb_item_list = 0, type = 0, i = 0, err = 0;
  double tmpDbl = 0.0, tmpDbl1 = 0.0, tmpDbl2 = 0.0, * matrix_dbl = NULL;
  int tmpInt1 = 0, tmpInt2 = 0;
  SciErr _sciErr;

  // Arg 1: a list
  // Arg 2: border name
  // Example:
  // bord = tlist(['border';'a';'b';'c';'d';'e';'f'],...
  //               list('x = t; y = 0',0,1,6,1),...
  //               list('x = 1; y = t',0,0.5,4, 1),...
  //               list('x = 1 - t; y = 0.5',0,0.5,4,1),...
  //               list('x = 0.5; y = t',0.5,1,4,1),...
  //               list('x = 1 - t; y = 1',0.5,1,4,1),...
  //               list('x = 0; y = 1 - t',0,1,6,1));
  //
  // buildMesh(bord,'th'); //mesh building
  //
  // Other example:
  //
  // //Definition of the border
  // deff('[x,y] = f1(t)','x = 0; y = 1 - t');
  //
  // t = 0:1/6:1;
  // A = [zeros(t); 1-t]';
  //
  // bord = tlist(['border';'a';'b';'c';'d';'e';'f'],...
  //               list('x = t; y = 0',0,1,6,1),...
  //               list('x = 1; y = t',0,0.5,4, 1),...
  //               list('x = 1 - t; y = 0.5',0,0.5,4,1),...
  //               list('x = 0.5; y = t',0.5,1,4,1),...
  //               list('x = 1 - t; y = 1',0.5,1,4,1),...
  //               list(A, 1));
  //
  // buildMesh(bord,'th');

  CheckRhs(minrhs, maxrhs);
  CheckLhs(minlhs, maxlhs);

  _sciErr = getVarAddressFromPosition(pvApiCtx, 1, &list_addr);
  _sciErr = getListItemNumber(pvApiCtx, list_addr, &nb_item_list);

  _sciErr = getListItemAddress(pvApiCtx, list_addr, 1, &item_list_addr);
  getAllocatedMatrixOfString(pvApiCtx, item_list_addr, &m_string_list, &n_string_list, &string_list);

  C2F(delete_lb_lobj)();
  C2F(delete_lexp)();

  for(i=2; i<=nb_item_list; i++)
    {
      strcpy(bordername, string_list[i-1]);

      // The first element is a matrix of string, we skip it (i=2), now we get a list
      _sciErr = getListInList(pvApiCtx, list_addr, i, &item_list_addr);
      // We get the first item to check if it's a value (type = 1) or a string (type = 10)
      _sciErr = getListItemAddress(pvApiCtx, item_list_addr, 1, &item_item_list_addr);
      _sciErr = getVarType(pvApiCtx, item_item_list_addr, &type);

      if (type == 1) // real of complex matrix
	{
	  // Case where we have a list of values
	  // deff('[x,y] = f1(t)','x = 0; y = 1 - t');
	  // t = 0:1/6:1;
	  // A = [zeros(t), 1-t]';
	  // list(A, 1));

	  // the first real matrix
	  _sciErr = getListItemAddress(pvApiCtx, item_list_addr, 1, &item_item_list_addr);
	  _sciErr = getMatrixOfDouble(pvApiCtx, item_item_list_addr, &m_matrix_dbl, &n_matrix_dbl, &matrix_dbl);
	  nb_points = m_matrix_dbl*n_matrix_dbl / 2;

	  // the second value (the first integer value)
	  _sciErr = getListItemAddress(pvApiCtx, item_list_addr, 2, &item_item_list_addr);
	  getScalarDouble(pvApiCtx, item_item_list_addr, &tmpDbl);
	  tmpInt1 = (int)tmpDbl;

	  buffer[0] = ' ';
	  C2F(put_scilab_border_1)(        bordername,       matrix_dbl,       &nb_points,   &tmpInt1,     buffer,        &err);
	  // void C2F(put_scilab_border_1)(char * nameBorder,double * tBorder, int * nbPoint,int * normal, char * Message,int * err);
	  
	  if (err)
	    {
	      Scierror(999,"%s: Error while putting a border defined by coordinates.\n", fname);
	      freeAllocatedMatrixOfString(m_string_list, n_string_list, string_list);
	      return 0;
	    }
	}
      else if (type == 10) // matrix of string
	{
	  // Case where the freefem syntax is given:
	  // list('x = 1 - t; y = 0.5',0,0.5,4,1)

	  // We get fct = 'x = 1 - t; y = 0.5'
	  _sciErr = getListItemAddress(pvApiCtx, item_list_addr, 1, &item_item_list_addr);
	  _sciErr = getVarDimension(pvApiCtx, item_item_list_addr, &n_param, &m_param);
	  
	  if (n_param*m_param!=1)
	    {
	      Scierror(999,"%s: Error, a single string expected.\n", fname);
	      freeAllocatedMatrixOfString(m_string_list, n_string_list, string_list);
	      return 0;
	    }

	  getAllocatedSingleString(pvApiCtx, item_item_list_addr, &fct);
	  
	  // the first real value
	  _sciErr = getListItemAddress(pvApiCtx, item_list_addr, 2, &item_item_list_addr);
	  getScalarDouble(pvApiCtx, item_item_list_addr, &tmpDbl1);
	  
	  // the second real value
	  _sciErr = getListItemAddress(pvApiCtx, item_list_addr, 3, &item_item_list_addr);
	  getScalarDouble(pvApiCtx, item_item_list_addr, &tmpDbl2);
	  
	  // the third value (the first integer value)
	  _sciErr = getListItemAddress(pvApiCtx, item_list_addr, 4, &item_item_list_addr);
	  getScalarDouble(pvApiCtx, item_item_list_addr, &tmpDbl);
	  tmpInt1 = (int)tmpDbl;

	  // the fourth value (the second integer value)
	  _sciErr = getListItemAddress(pvApiCtx, item_list_addr, 5, &item_item_list_addr);
	  getScalarDouble(pvApiCtx, item_item_list_addr, &tmpDbl);
	  tmpInt2 = (int)tmpDbl;
	  
	  buffer[0] = ' ';
	  
	  C2F(put_scilab_border_2)(        bordername,        fct,             &tmpDbl1,       &tmpDbl2,      &tmpInt1, &tmpInt2,  buffer,         &err);
	  // void C2F(put_scilab_border_2)(char * nameBorder, char * equation, double * start, double * stop, int * ng, int * nrm, char * Message, int * err);

	  if (err)
	    {
	      Scierror(999,"%s: Error while putting a border defined by a freefem expression: %s.\n", fname, fct);
	      freeAllocatedMatrixOfString(m_string_list, n_string_list, string_list);
	      freeAllocatedSingleString(fct);
	      return 0;
	    }

	  freeAllocatedSingleString(fct);
	}
      else
	{
	  Scierror(999,"%s: Error: a matrix + int or a freefem expression requested.\n", fname);
	  freeAllocatedMatrixOfString(m_string_list, n_string_list, string_list);
	  return 0;
	}
    }
  
  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &name_mesh_addr);
  _sciErr = getVarDimension(pvApiCtx, item_item_list_addr, &n_param, &m_param);
  
  if (n_param*m_param!=1)
    {
      Scierror(999,"%s: Error, a single string expected.\n", fname);
      freeAllocatedMatrixOfString(m_string_list, n_string_list, string_list);
      return 0;
    }

  getAllocatedSingleString(pvApiCtx, name_mesh_addr, &name_mesh);

  buffer[0] = ' ';
  C2F(build_scilab_mesh)(       name_mesh,       buffer,        &err);
  //void C2F(build_scilab_mesh)(char * nameMesh, char * Message,int * err)

  freeAllocatedSingleString(name_mesh);

  if (err)
    {
      Scierror(999, "%s: Error while building mesh.\n", fname);
      return 0;
    }

  return 0;
}

int int_put_new_mesh(char * fname)
{
  char bordername[20], * fct = NULL;
  int minrhs = 2, maxrhs = 2;
  int minlhs = 1, maxlhs = 1;
  int m_list, n_list, * list_addr = NULL, * item_list_addr = NULL, * item_item_list_addr = NULL;
  int m_name_mesh, n_name_mesh, * name_mesh_addr = NULL;
  int m_matrix_dbl, n_matrix_dbl, nb_points = 0;
  int m_string_list, n_string_list;
  char * name_mesh = NULL, buffer[4096];
  char ** string_list = NULL;
  int nb_item_list = 0, type = 0, i = 0, err = 0;
  double tmpDbl = 0.0, tmpDbl1 = 0.0, tmpDbl2 = 0.0, * matrix_dbl = NULL;
  int tmpInt1 = 0, tmpInt2 = 0;
  SciErr _sciErr;

  CheckRhs(minrhs, maxrhs);
  CheckLhs(minlhs, maxlhs);

  // Arg 1: a list
  // Arg 2: border name
  // Example of use:
  // //Definition of the border
  // deff('[x,y] = f1(t)','x = 0; y = 1 - t');
  //
  // t = 0:1/6:1;
  // A = [zeros(t); 1-t]';
  //
  // bord = tlist(['border';'a';'b';'c';'d';'e';'f'],...
  //               list('x = t; y = 0',0,1,6,1),...
  //               list('x = 1; y = t',0,0.5,4, 1),...
  //               list('x = 1 - t; y = 0.5',0,0.5,4,1),...
  //               list('x = 0.5; y = t',0.5,1,4,1),...
  //               list('x = 1 - t; y = 1',0.5,1,4,1),...
  //               list(A, 1));
  //
  // ...
  //
  // // Update mesh
  // deff('[x,y] = f1(t)','x = (0.5/(4-i))*cos(t); y = 0.5+0.5*sin(t)');
  // bord('f')(1) = evalf(f1,%pi/2:%pi/(5*i):3*%pi/2);
  // UpdateMesh(bord,'th');

  _sciErr = getVarAddressFromPosition(pvApiCtx, 1, &list_addr);
  _sciErr = getListItemNumber(pvApiCtx, list_addr, &nb_item_list);

  _sciErr = getListItemAddress(pvApiCtx, list_addr, 1, &item_list_addr);
  getAllocatedMatrixOfString(pvApiCtx, item_list_addr, &m_string_list, &n_string_list, &string_list);
      
  C2F(delete_lb_lobj)();

  for(i=2; i<=nb_item_list; i++)
    {
      strcpy(bordername, string_list[i-1]);

      // The first element is a matrix of string, we skip it (i=2), now we get a list
      _sciErr = getListInList(pvApiCtx, list_addr, i, &item_list_addr);
      // We get the first item to check if it's a value (type = 1) or a string (type = 10)
      _sciErr = getListItemAddress(pvApiCtx, item_list_addr, 1, &item_item_list_addr);
      _sciErr = getVarType(pvApiCtx, item_item_list_addr, &type);
      
      if (type == 1) // real or complex matrix
	{
	  // the first real matrix
	  _sciErr = getListItemAddress(pvApiCtx, item_list_addr, 1, &item_item_list_addr);
	  _sciErr = getMatrixOfDouble(pvApiCtx, item_item_list_addr, &m_matrix_dbl, &n_matrix_dbl, &matrix_dbl);
	  nb_points = m_matrix_dbl*n_matrix_dbl / 2;

	  // the second value (the first integer value)
	  _sciErr = getListItemAddress(pvApiCtx, item_list_addr, 2, &item_item_list_addr);
	  getScalarDouble(pvApiCtx, item_item_list_addr, &tmpDbl);
	  tmpInt1 = (int)tmpDbl;

	  buffer[0] = ' ';
	  C2F(put_scilab_border_1)(        bordername,       matrix_dbl,       &nb_points,   &tmpInt1,     buffer,        &err);
	  // void C2F(put_scilab_border_1)(char * nameBorder,double * tBorder, int * nbPoint,int * normal, char * Message,int * err);

	  if (err)
	    {
	      Scierror(999,"%s: Error while putting a border defined by coordinates.\n", fname);
	      freeAllocatedMatrixOfString(m_string_list, n_string_list, string_list);
	      return 0;
	    }
	}
      else if (type == 10) // matrix of string
	{
	  // Case where we give the freefem syntax
	  C2F(add_lobj)();
	}
      else
	{
	  Scierror(999,"%s: Error: a matrix + int or a freefem expression requested.\n", fname);
	  freeAllocatedMatrixOfString(m_string_list, n_string_list, string_list);
	  return 0;
	}
    }

  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &name_mesh_addr);
  getAllocatedSingleString(pvApiCtx, name_mesh_addr, &name_mesh);

  buffer[0] = ' ';
  C2F(build_scilab_mesh)(        name_mesh,      buffer,        &err);
  // void C2F(build_scilab_mesh)(char * nameMesh,char * Message,int * err)

  freeAllocatedSingleString(name_mesh);
  freeAllocatedMatrixOfString(m_string_list, n_string_list, string_list);

  if (err)
    {
      Scierror(999, "%s: Error while updating mesh.\n", fname);
      return 0;
    }

  return 0;
}

int int_get_matrix(char * fname)
{
  double * lMatrix = NULL;
  int * ljlow = NULL,  * ljhigh = NULL, * lblock = NULL, mMatrix = 0;
  int minrhs = 0, maxrhs = 0;
  int minlhs = 4, maxlhs = 4;
  int err = 0, lsize, i;
  char buffer[4096];
  double * tmpDblPtr = NULL, tmpDbl = 0.0;
  SciErr _sciErr;

  CheckRhs(minrhs, maxrhs);
  CheckLhs(minlhs, maxlhs);

  C2F(get_matrix)(&lMatrix, &mMatrix,
		  &ljlow, &ljhigh, &lsize,
		  &lblock, buffer, &err);

  if (err)
    {
      Scierror(999, "%s: ???\n", fname);
      return 0;
    }

  _sciErr = createMatrixOfDouble(pvApiCtx, 1, mMatrix, 1, lMatrix);
  C2F(delete_ptr)(lMatrix);

  tmpDblPtr = (double *)MALLOC(lsize*sizeof(double));
  for(i=0; i<lsize; i++) tmpDblPtr[i] = ljlow[i];
  _sciErr = createMatrixOfDouble(pvApiCtx, 2, lsize, 1, tmpDblPtr);
  C2F(delete_ptr)(ljlow);
  FREE(tmpDblPtr);

  tmpDblPtr = (double *)MALLOC(lsize*sizeof(double));
  for(i=0; i<lsize; i++) tmpDblPtr[i] = ljhigh[i];
  _sciErr = createMatrixOfDouble(pvApiCtx, 3, lsize, 1, tmpDblPtr);
  C2F(delete_ptr)(ljhigh);
  FREE(tmpDblPtr);

  tmpDbl = (double)lblock[0];
  createScalarDouble(pvApiCtx, 4, tmpDbl);
  C2F(delete_ptr)(lblock);
  
  LhsVar(1) = 1;
  LhsVar(2) = 2;
  LhsVar(3) = 3;
  LhsVar(4) = 4;

  return 0;
}
