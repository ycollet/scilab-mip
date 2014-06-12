////////////////////////////////////////////////////////////
// read_lp_file, read_lp_file_mp: Tools to read mps files //
////////////////////////////////////////////////////////////

//  Copyright (C) 2008-2010 Yann Collette.
//
//  read_lp_file is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2, or (at your option) any
//  later version.
//
//  This part of code is distributed with the FURTHER condition that it 
//  can be compiled and linked with the Matlab libraries and it can be 
//  used within the Matlab environment.
//
//  read_lp_file is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with SciCoinOR; see the file COPYING.  If not, write to the Free
//  Software Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#include <stdio.h>
#include <fstream>

extern "C"
{
#include <stack-c.h>
#include <sciprint.h>
#include <Scierror.h>
#include <MALLOC.h>
#include <freeArrayOfString.h>
}

#include <api_scilab.h>

#include <CoinPackedMatrix.hpp>
#include <CoinLpIO.hpp>

//#define DEBUG 1
//#define DBGPRINTF printf
#define DBGPRINTF sciprint

#include <helper.hpp>

using namespace std;

// Add columnName(int), rowName(int), getProblemName(), getObjName()

// read in all necessary elements for retrieving the LP/MILP
extern "C" int read_lp_file(char * fname)
{  
  int * filename_addr = NULL;          // char - Input
  int nb_constr;
  int nb_obj_var;
  int nb_val_constr_mat;
  int nb_int_var;
  int status, i;
  static char * ListLabels[] = {"plist","nb_constr","nb_obj_var","nb_val_constr_mat","nb_int_var"};
  CoinLpIO      LpReader;
  char * filename = NULL;
  SciErr _SciErr;

#ifdef DEBUG
  DBGPRINTF("DEBUG: number of parameters = %d\n", Rhs);
#endif

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &filename_addr); SCICOINOR_ERROR;
  getAllocatedSingleString(pvApiCtx, filename_addr, &filename);

#ifdef DEBUG
  DBGPRINTF("DEBUG: filename = %s\n", filename);
#endif

  ifstream infile;
  infile.open(filename,ifstream::in);
  if (!infile.is_open())
    {
      sciprint("read_lp_file: error while opening file %s\n", filename);
      freeAllocatedSingleString(filename);
      return 0;
    }
  infile.close();

  LpReader.readLp(filename);

  // retrieve number of constraints
  nb_constr = LpReader.getNumRows();
  
  // retrieve number of objective variables
  nb_obj_var = LpReader.getNumCols();
  
  // retrieve number of non-zero elements in constraint matrix
  nb_val_constr_mat = LpReader.getNumElements();
  
  // retrieve number of integer variables
  status = 0;
  for(i=0;i<nb_obj_var;i++)
    {
      if (LpReader.isInteger(i)) status++;
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

  int * list_addr = NULL;
  double tmp_dbl;

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
extern "C" int read_lp_file_mp(char * fname)
{
  int * filename_addr;  // char - Input

  int nb_constr;
  int nb_obj_var;
  int nb_val_constr_mat;
  int index_start;

  char ** RowNames = NULL;
  char ** ColNames = NULL;
  char * filename  = NULL;
  double * obj_var_is_int = NULL;

  static char * ListLabels [] = {"plist","constr_mat","constr_sense","obj_coeff",
				 "rhs","lhs","bounds_lower","bounds_upper","obj_var_is_int",
				 "pb_name","col_name","row_name"};

  CoinLpIO           LpReader;
  static SciSparse * ConstrMat = NULL;
  int i, j;
  SciErr _SciErr;

#ifdef DEBUG
  DBGPRINTF("DEBUG: number of parameters = %d\n", Rhs);
#endif

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &filename_addr); SCICOINOR_ERROR;
  getAllocatedSingleString(pvApiCtx, filename_addr, &filename);

  ifstream infile;
  infile.open(filename,ifstream::in);
  if (!infile.is_open())
    {
      sciprint("read_lp_file_mp: error while opening file %s\n", filename);
      freeAllocatedSingleString(filename);
      return 0;
    }
  infile.close();

  LpReader.readLp(filename);

  // retrieve number of constraints
  nb_constr = LpReader.getNumRows();
  
  // retrieve number of objective variables
  nb_obj_var = LpReader.getNumCols();
  
  // retrieve number of non-zero elements in constraint matrix
  nb_val_constr_mat = LpReader.getNumElements();

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
      return 0;
    }
  
  // Copy the sense of constraintes
  // 'L' - smaller than - <=
  // 'E' - equality     - =
  // 'G' - greater than - >=
  // 'R' - Range        - <= + >=
  // 'N' - Free         - no constraints

  // Copy Row name and column names
  RowNames = (char **)MALLOC(sizeof(char *) * LpReader.getNumRows());
  for(i=0;i<LpReader.getNumRows();i++)
    {
      RowNames[i] = (char *)MALLOC(sizeof(char) * (strlen(LpReader.rowName(i))+1));
      strncpy(RowNames[i],LpReader.rowName(i),strlen(LpReader.rowName(i))+1);
    }

  ColNames = (char **)MALLOC(sizeof(char *) * LpReader.getNumCols());
  for(i=0;i<LpReader.getNumCols();i++)
    {
      ColNames[i] = (char *)MALLOC(sizeof(char) * (strlen(LpReader.columnName(i))+1));
      strncpy(ColNames[i],LpReader.columnName(i),strlen(LpReader.columnName(i))+1);
    }

#ifdef DEBUG
  sciprint("constr_sense = %s\n", LpReader.getRowSense());
#endif
  
  // retrieve column specific data (values, bounds and type)
  obj_var_is_int = (double *)MALLOC(nb_obj_var*sizeof(double));

  for (i=0; i<nb_obj_var; i++) 
    {
      // set to TRUE if objective variable is integer or binary  
      if (LpReader.isInteger(i)) *(obj_var_is_int+i) = 1;
      else                       *(obj_var_is_int+i) = 0;
    }
      
  // the constraint matrix
  for (i=0; i<LpReader.getMatrixByRow()->getSizeVectorLengths(); i++) 
    {
      ConstrMat->mnel[i] = LpReader.getMatrixByRow()->getVectorLengths()[i];
      index_start = LpReader.getMatrixByRow()->getVectorStarts()[i];
#ifdef DEBUG
      DBGPRINTF("mat[%d] : mnel = %d index_start = %d\n",i,ConstrMat->mnel[i],index_start);
#endif
      for(j=0;j<ConstrMat->mnel[i];j++)
	{
	  ConstrMat->icol[index_start+j] = LpReader.getMatrixByRow()->getIndices()[index_start+j] + 1;
	  ConstrMat->R[index_start+j]    = LpReader.getMatrixByRow()->getElements()[index_start+j];
#ifdef DEBUG
	  DBGPRINTF("mat[%d][%d] : icol = %d R = %f\n",i,j,ConstrMat->icol[index_start+j],ConstrMat->R[index_start+j]);
#endif
	}
    }

#ifdef DEBUG
  DBGPRINTF("DEBUG: nb_constr  = %d\n", nb_constr);
  DBGPRINTF("DEBUG: nb_obj_var = %d\n", nb_obj_var);
  DBGPRINTF("DEBUG: nb_val_constr_mat = %d\n", nb_val_constr_mat);

  DBGPRINTF("DEBUG: constrmat_it  = %d\n", constrmat_it);
  DBGPRINTF("DEBUG: constrmat_m   = %d\n", constrmat_m);
  DBGPRINTF("DEBUG: constrmat_n   = %d\n", constrmat_n);
  DBGPRINTF("DEBUG: constrmat_nel = %d\n", constrmat_nel);

  DBGPRINTF("DEBUG: ConstrMat->it  = %d\n", ConstrMat->it);
  DBGPRINTF("DEBUG: ConstrMat->m   = %d\n", ConstrMat->m);
  DBGPRINTF("DEBUG: ConstrMat->n   = %d\n", ConstrMat->n);
  DBGPRINTF("DEBUG: ConstrMat->nel = %d\n", ConstrMat->nel);
  DBGPRINTF("Exiting ...\n");
#endif

  // Now create the mlist of type plist
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
  tmp_char = (char *)LpReader.getRowSense();
  _SciErr = createMatrixOfStringInList(pvApiCtx, Rhs+1, list_addr, 3, 1, 1, &tmp_char); SCICOINOR_ERROR;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, list_addr, 4, 1, nb_obj_var, (double *)LpReader.getObjCoefficients()); SCICOINOR_ERROR;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, list_addr, 5, 1, nb_constr,  (double *)LpReader.getRowUpper()); SCICOINOR_ERROR;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, list_addr, 6, 1, nb_constr,  (double *)LpReader.getRowLower()); SCICOINOR_ERROR;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, list_addr, 7, 1, nb_obj_var, (double *)LpReader.getColLower()); SCICOINOR_ERROR;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, list_addr, 8, 1, nb_obj_var, (double *)LpReader.getColUpper()); SCICOINOR_ERROR;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, list_addr, 9, 1, nb_obj_var, obj_var_is_int); SCICOINOR_ERROR;
  tmp_char = (char *)LpReader.getProblemName();
  _SciErr = createMatrixOfStringInList(pvApiCtx, Rhs+1, list_addr, 10, 1, 1, &tmp_char); SCICOINOR_ERROR;
  _SciErr = createMatrixOfStringInList(pvApiCtx, Rhs+1, list_addr, 11, 1, nb_obj_var, ColNames); SCICOINOR_ERROR;
  _SciErr = createMatrixOfStringInList(pvApiCtx, Rhs+1, list_addr, 12, 1, nb_constr,  RowNames); SCICOINOR_ERROR;

  LhsVar(1) = Rhs+1;

  if (ConstrMat->mnel) FREE(ConstrMat->mnel);
  if (ConstrMat->icol) FREE(ConstrMat->icol);
  if (ConstrMat->R)    FREE(ConstrMat->R);
  if (ConstrMat)       FREE(ConstrMat);

  freeArrayOfString(RowNames,LpReader.getNumRows());
  freeArrayOfString(ColNames,LpReader.getNumCols());

  freeAllocatedSingleString(filename);
  FREE(obj_var_is_int);

  return 0;
}
