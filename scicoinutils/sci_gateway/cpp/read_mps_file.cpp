//////////////////////////////////////////////////////////////
// read_mps_file, read_mps_file_mp: Tools to read mps files //
//////////////////////////////////////////////////////////////

//  Copyright (C) 2008-2010 Yann Collette.
//
//  read_mps_file is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2, or (at your option) any
//  later version.
//
//  read_mps_file is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with SciCoinOR; see the file COPYING.  If not, write to the Free
//  Software Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#include <stdio.h>

extern "C"
{
#include <stack-c.h>
#include <sciprint.h>
#include <Scierror.h>
#include <MALLOC.h>
#include <freeArrayOfString.h>
}

#include <api_scilab.h>

#include <CoinMpsIO.hpp>

//#define DEBUG 1
//#define DBGPRINTF printf
#define DBGPRINTF sciprint

#include <helper.hpp>

// read in all necessary elements for retrieving the LP/MILP
extern "C" int read_mps_file(char * fname)
{  
  int * filename_addr = NULL, * type_addr = NULL, * verbosity_addr = NULL, * list_addr = NULL;
  int nb_constr;
  int nb_obj_var;
  int nb_val_constr_mat;
  int nb_int_var;
  int i;
  double status, verbosity, type, tmp_dbl;
  static char * ListLabels[] = {"plist","nb_constr","nb_obj_var","nb_val_constr_mat","nb_int_var"};
  CoinMpsIO        mpsReader;
  DerivedHandler * printer = NULL;
  char * filename = NULL;
  SciErr _SciErr;

#ifdef DEBUG
  DBGPRINTF("DEBUG: number of parameters = %d\n", Rhs);
#endif

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &filename_addr); SCICOINOR_ERROR;
  getAllocatedSingleString(pvApiCtx, filename_addr, &filename);

  if (Rhs>=2)
    {
      _SciErr = getVarAddressFromPosition(pvApiCtx, 2, &type_addr); SCICOINOR_ERROR;
      if (isEmptyMatrix(pvApiCtx, type_addr))
	{
	  type = 0;
	}
      else
	{
	  getScalarDouble(pvApiCtx, type_addr, &type);
	}
      // 0: mps
      // 1: gms  (GAMS)
      // 2: gmpl (AMPL)
    }
  else
    {
      type = 0; // mps
    }
  if (Rhs>=3)
    {
      _SciErr = getVarAddressFromPosition(pvApiCtx, 3, &verbosity_addr); SCICOINOR_ERROR;
      if (isEmptyMatrix(pvApiCtx, verbosity_addr))
	{
	  verbosity = 0;
	}
      else
	{
	  getScalarDouble(pvApiCtx, type_addr, &verbosity);
	}
    }
  else
    {
      verbosity = 0;
    }

#ifdef DEBUG
  DBGPRINTF("DEBUG: filename = %s, verbosity = %d\n", filename, verbosity);
#endif

  // Turn on/off Terminal Output
  if(verbosity!=0) 
    {
      printer = new DerivedHandler(); // assumed open	
      printer->setLogLevel((int)verbosity);		 
    }

  switch((int)type)
    {
    case 1:
      if (verbosity) sciprint("read_mps_file: reading a .gms file\n");
      status = mpsReader.readGms(filename);
      break;
    case 2:
      if (verbosity) sciprint("read_mps_file: reading a .gmpl file\n");
      status = mpsReader.readGMPL(filename);
      break;
    default:
      if (verbosity) sciprint("read_mps_file: reading a .mps file\n");
      status = mpsReader.readMps(filename);
      break;
    }

  if (status) 
    {
      Scierror(999, "%s: Reading file %c failed\n", fname, filename);
      freeAllocatedSingleString(filename);
      return 0;
    }

  // retrieve number of constraints
  nb_constr = mpsReader.getNumRows();
  
  // retrieve number of objective variables
  nb_obj_var = mpsReader.getNumCols();
  
  // retrieve number of non-zero elements in constraint matrix
  nb_val_constr_mat = mpsReader.getNumElements();
  
  // retrieve number of integer variables
  status = 0;
  for(i=0;i<nb_obj_var;i++)
    {
      if (mpsReader.isInteger(i)) status++;
    }
  nb_int_var = status;

#ifdef DEBUG
  DBGPRINTF("DEBUG: nb_constr         = %d\n", nb_constr);
  DBGPRINTF("DEBUG: nb_obj_var        = %d\n", nb_obj_var);
  DBGPRINTF("DEBUG: nb_val_constr_mat = %d\n", nb_val_constr_mat);
  DBGPRINTF("DEBUG: nb_int_var        = %d\n", nb_int_var);
#endif

  //////////////////////////////////////////////
  // Creation  of the output scilab variables //
  //////////////////////////////////////////////

  // Now create scilab variables which will be stored in the list via a call to CreateListVarFrom

  _SciErr = createMList(pvApiCtx, Rhs+1, 5, &list_addr); SCICOINOR_ERROR;
  _SciErr = createMatrixOfStringInList(pvApiCtx, Rhs+1, list_addr, 1, 1, 5, ListLabels); SCICOINOR_ERROR;
  tmp_dbl = (double)nb_constr;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, list_addr, 2, 1, 1, &tmp_dbl); SCICOINOR_ERROR;
  tmp_dbl = (double)nb_obj_var;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, list_addr, 3, 1, 1, &tmp_dbl); SCICOINOR_ERROR;
  tmp_dbl = (double)nb_val_constr_mat;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, list_addr, 4, 1, 1, &tmp_dbl); SCICOINOR_ERROR;
  tmp_dbl = (double)nb_int_var;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, list_addr, 5, 1, 1, &tmp_dbl); SCICOINOR_ERROR;

  LhsVar(1) = Rhs+1;

  freeAllocatedSingleString(filename);

  return 0;
}

// retrieve all missing values of LP/MILP
extern "C" int read_mps_file_mp(char * fname)
{
  int nb_constr;
  int nb_obj_var;
  int nb_val_constr_mat;
  int index_start;

  int * filename_addr = NULL, * type_addr = NULL, * verbosity_addr = NULL;

  char ** RowNames = NULL;
  char ** ColNames = NULL;
  static char * ListLabels [] = {"plist","constr_mat","constr_sense","obj_coeff",
				 "rhs","lhs","bounds_lower","bounds_upper","obj_var_is_int",
				 "pb_name","col_name","row_name"};

  CoinMpsIO        mpsReader;
  DerivedHandler * printer = NULL;

  static SciSparse * ConstrMat = NULL;
  int i, j;
  double status, verbosity, type;
  double * obj_var_is_int = NULL;
  char * filename = NULL;
  SciErr _SciErr;

#ifdef DEBUG
  DBGPRINTF("DEBUG: number of parameters = %d\n", Rhs);
#endif

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &filename_addr); SCICOINOR_ERROR;
  getAllocatedSingleString(pvApiCtx, filename_addr, &filename);

  if (Rhs>=2) 
    {
      _SciErr = getVarAddressFromPosition(pvApiCtx, 2, &type_addr); SCICOINOR_ERROR;
      if (isEmptyMatrix(pvApiCtx, type_addr))
	{
	  type = 0;
	}
      else
	{
	  getScalarDouble(pvApiCtx, type_addr, &type);
	}
    }
  else
    {
      type = 0;
    }
  if (Rhs>=3) 
    {
      _SciErr = getVarAddressFromPosition(pvApiCtx, 3, &verbosity_addr); SCICOINOR_ERROR;
      if (isEmptyMatrix(pvApiCtx, verbosity_addr))
	{
	  verbosity = 0;
	}
      else
	{
	  getScalarDouble(pvApiCtx, type_addr, &verbosity);
	}
    }
  else
    {
      verbosity = 0;
    }

  // Turn on/off Terminal Output
  if(verbosity!=0) 
    {
      printer = new DerivedHandler(); // assumed open	
      printer->setLogLevel((int)verbosity);		 
    }

  switch((int)type)
    {
    case 1:
      if (verbosity) sciprint("read_mps_file: reading a .gms file\n");
      status = mpsReader.readGms(filename);
      break;
    case 2:
      if (verbosity) sciprint("read_mps_file: reading a .gmpl file\n");
      status = mpsReader.readGMPL(filename);
      break;
    default:
      if (verbosity) sciprint("read_mps_file: reading a .mps file\n");
      status = mpsReader.readMps(filename);
      break;
    }

  if (status) 
    {
      Scierror(999, "%s: Reading file %c failed\n", fname, filename);
      freeAllocatedSingleString(filename);
      return 0;
    }

  // retrieve number of constraints
  nb_constr = mpsReader.getNumRows();
  
  // retrieve number of objective variables
  nb_obj_var = mpsReader.getNumCols();
  
  // retrieve number of non-zero elements in constraint matrix
  nb_val_constr_mat = mpsReader.getNumElements();

  ///////////////////////////////////////////
  // Create of the output scilab variables //
  ///////////////////////////////////////////

  // First, allocate a new sparse matrix
  ConstrMat      = (SciSparse *)MALLOC(1*sizeof(SciSparse));

  ConstrMat->n    = nb_obj_var;
  ConstrMat->m    = nb_constr;
  ConstrMat->it   = 0;
  ConstrMat->nel  = nb_val_constr_mat;
  ConstrMat->mnel = (int *)MALLOC(nb_constr*sizeof(int));
  ConstrMat->icol = (int *)MALLOC(nb_val_constr_mat*sizeof(int));
  ConstrMat->R    = (double *)MALLOC(nb_val_constr_mat*sizeof(double));

  if ((ConstrMat==(SciSparse *)0) || (ConstrMat->mnel==NULL) || (ConstrMat->R==NULL))
    {
      Scierror(999, "%s: error while allocating the sparse\n",fname);
      freeAllocatedSingleString(filename);
      return 0;
    }

  // Copy the sense of constraintes
  // 'L' - smaller than - <=
  // 'E' - equality     - =
  // 'G' - greater than - >=
  // 'R' - Range        - <= + >=
  // 'N' - Free         - no constraints

  // Copy Row name and column names
  RowNames = (char **)MALLOC(sizeof(char *) * mpsReader.getNumRows());
  for(i=0;i<mpsReader.getNumRows();i++)
    {
      RowNames[i] = (char *)MALLOC(sizeof(char) * (strlen(mpsReader.rowName(i))+1));
      strncpy(RowNames[i],mpsReader.rowName(i),strlen(mpsReader.rowName(i))+1);
    }

  ColNames = (char **)MALLOC(sizeof(char *) * mpsReader.getNumCols());
  for(i=0;i<mpsReader.getNumCols();i++)
    {
      ColNames[i] = (char *)MALLOC(sizeof(char) * (strlen(mpsReader.columnName(i))+1));
      strncpy(ColNames[i],mpsReader.columnName(i),strlen(mpsReader.columnName(i))+1);
    }
  
  // retrieve column specific data (values, bounds and type)
  obj_var_is_int = (double *)MALLOC(nb_obj_var*sizeof(double));

  for (i=0; i<nb_obj_var; i++) 
    {
      // set to TRUE if objective variable is integer or binary  
      if (mpsReader.isInteger(i)) *(obj_var_is_int+i) = 1.0;
      else                        *(obj_var_is_int+i) = 0.0;
    }
  
  // the constraint matrix
  for (i=0; i<mpsReader.getMatrixByRow()->getSizeVectorLengths(); i++) 
    {
      ConstrMat->mnel[i] = mpsReader.getMatrixByRow()->getVectorLengths()[i];
      index_start = mpsReader.getMatrixByRow()->getVectorStarts()[i];
#ifdef DEBUG
      DBGPRINTF("mat[%d] : mnel = %d index_start = %d\n",i,ConstrMat->mnel[i],index_start);
#endif
      for(j=0;j<ConstrMat->mnel[i];j++)
	{
	  ConstrMat->icol[index_start+j] = mpsReader.getMatrixByRow()->getIndices()[index_start+j] + 1;
	  ConstrMat->R[index_start+j]    = mpsReader.getMatrixByRow()->getElements()[index_start+j];
#ifdef DEBUG
	  DBGPRINTF("mat[%d][%d] : icol = %d R = %f\n",i,j,ConstrMat->icol[index_start+j],ConstrMat->R[index_start+j]);
#endif
	}
    }

#ifdef DEBUG
  DBGPRINTF("DEBUG: nb_constr  = %d\n", nb_constr);
  DBGPRINTF("DEBUG: nb_obj_var = %d\n", nb_obj_var);
  DBGPRINTF("DEBUG: nb_val_constr_mat = %d\n", nb_val_constr_mat);

  DBGPRINTF("DEBUG: ConstrMat->it  = %d\n", ConstrMat->it);
  DBGPRINTF("DEBUG: ConstrMat->m   = %d\n", ConstrMat->m);
  DBGPRINTF("DEBUG: ConstrMat->n   = %d\n", ConstrMat->n);
  DBGPRINTF("DEBUG: ConstrMat->nel = %d\n", ConstrMat->nel);
  DBGPRINTF("Exiting ...\n");
#endif

  // Create the scilab arrays to store the parameters of the problem
  int * list_addr = NULL;
  char * tmp_char = NULL;

  _SciErr = createMList(pvApiCtx, Rhs+1, 12, &list_addr); SCICOINOR_ERROR;
  _SciErr = createMatrixOfStringInList(pvApiCtx, Rhs+1, list_addr, 1, 1, 12, ListLabels); SCICOINOR_ERROR;
  _SciErr = createSparseMatrixInList(pvApiCtx, Rhs+1, list_addr, 2, 
				     ConstrMat->m, 
				     ConstrMat->n, 
				     ConstrMat->nel,  
				     ConstrMat->mnel, 
				     ConstrMat->icol,
				     ConstrMat->R); SCICOINOR_ERROR;
  tmp_char = (char *)mpsReader.getRowSense();
  _SciErr = createMatrixOfStringInList(pvApiCtx, Rhs+1, list_addr, 3, 1, 1, &tmp_char); SCICOINOR_ERROR;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, list_addr, 4, 1, nb_obj_var, (double *)mpsReader.getObjCoefficients()); SCICOINOR_ERROR;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, list_addr, 5, 1, nb_constr,  (double *)mpsReader.getRowUpper()); SCICOINOR_ERROR;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, list_addr, 6, 1, nb_constr,  (double *)mpsReader.getRowLower()); SCICOINOR_ERROR;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, list_addr, 7, 1, nb_obj_var, (double *)mpsReader.getColLower()); SCICOINOR_ERROR;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, list_addr, 8, 1, nb_obj_var, (double *)mpsReader.getColUpper()); SCICOINOR_ERROR;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, list_addr, 9, 1, nb_obj_var, obj_var_is_int); SCICOINOR_ERROR;
  tmp_char = (char *)mpsReader.getProblemName();
  _SciErr = createMatrixOfStringInList(pvApiCtx, Rhs+1, list_addr, 10, 1, 1, &tmp_char); SCICOINOR_ERROR;
  _SciErr = createMatrixOfStringInList(pvApiCtx, Rhs+1, list_addr, 11, 1, nb_obj_var, ColNames); SCICOINOR_ERROR;
  _SciErr = createMatrixOfStringInList(pvApiCtx, Rhs+1, list_addr, 12, 1, nb_constr,  RowNames); SCICOINOR_ERROR;

  LhsVar(1) = Rhs+1;

  if (ConstrMat->mnel) free(ConstrMat->mnel);
  if (ConstrMat->icol) free(ConstrMat->icol);
  if (ConstrMat->R)    free(ConstrMat->R);
  if (ConstrMat)       free(ConstrMat);

  freeArrayOfString(RowNames,mpsReader.getNumRows());
  freeArrayOfString(ColNames,mpsReader.getNumCols());

  if (printer) delete printer;

  freeAllocatedSingleString(filename);
  FREE(obj_var_is_int);

  return 0;
}
