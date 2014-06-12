//////////////////////////////////////////////
// write_mps_file: Tools to write mps files //
//////////////////////////////////////////////

//  Copyright (C) 2008-2010 Yann Collette.
//
//  write_mps_file is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2, or (at your option) any
//  later version.
//
//  This part of code is distributed with the FURTHER condition that it 
//  can be compiled and linked with the Matlab libraries and it can be 
//  used within the Matlab environment.
//
//  write_mps_file is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with SciCoinOR; see the file COPYING.  If not, write to the Free
//  Software Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#include <stdio.h>
#include <string>
#include <vector>
#include <sstream>

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
#include <CoinPackedMatrix.hpp>

#define DEBUG 1
//#define DBGPRINTF printf
#define DBGPRINTF sciprint

#include <helper.hpp>

using namespace std;

#define C_IN        1
#define A_IN        2
#define LHS_IN      3
#define RHS_IN      4
#define UB_IN       5
#define LB_IN       6
#define BTYPE_IN    7
#define VARTYPE_IN  8
#define PB_NAME_IN  9
#define COL_NAME_IN 10
#define ROW_NAME_IN 11
#define FILENAME_IN 12
#define VERBOSE_IN  13
#define STATUS_OUT  14

// status = write_mps_file(c,a,lhs,rhs,ub,lb,btype,vartype,filename);
// btype is not yet used. But I keep this parameter for the future
extern "C" int write_mps_file(char * fname)
{
  int m_c,        n_c;
  int m_a,        n_a;
  int m_lhs,      n_lhs;
  int m_rhs,      n_rhs;
  int m_ub,       n_ub;
  int m_lb,       n_lb;
  int m_col_name, n_col_name;
  int m_row_name, n_row_name;
  vector<string>  RowNames, ColNames;
  stringstream    tmpstring;

  CoinMpsIO        mpsWriter;
  DerivedHandler * printer = NULL;
  CoinPackedMatrix A_matrix;
  SciSparse        S_A;
  char ** pColNames = NULL;
  char ** pRowNames = NULL;
  int i, j, status;
  int nrows, ncols, count, type;
  int * c_addr = NULL, * lhs_addr = NULL, * rhs_addr = NULL, * ub_addr = NULL, * lb_addr = NULL;
  int * btype_addr = NULL, * vartype_addr = NULL, * pb_name_addr = NULL, * a_addr = NULL;
  int * col_name_addr = NULL, * row_name_addr = NULL, * filename_addr = NULL, * verbose_addr = NULL;
  double * c = NULL, * lhs = NULL, * rhs = NULL, * ub = NULL, * lb = NULL, * a = NULL;
  double verbose = 0;
  char * btype = NULL, * vartype = NULL, * pb_name = NULL, * filename = NULL;
  SciErr _SciErr;

  _SciErr = getVarAddressFromPosition(pvApiCtx, C_IN, &c_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, c_addr, &n_c, &m_c, &c); SCICOINOR_ERROR;
  _SciErr = getVarAddressFromPosition(pvApiCtx, LHS_IN, &lhs_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, lhs_addr, &n_lhs, &m_lhs, &lhs); SCICOINOR_ERROR;
  _SciErr = getVarAddressFromPosition(pvApiCtx, RHS_IN, &rhs_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, rhs_addr, &n_rhs, &m_rhs, &rhs); SCICOINOR_ERROR;
  _SciErr = getVarAddressFromPosition(pvApiCtx, UB_IN, &ub_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, ub_addr, &n_ub, &m_ub, &ub); SCICOINOR_ERROR;
  _SciErr = getVarAddressFromPosition(pvApiCtx, LB_IN, &lb_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, lb_addr, &n_lb, &m_lb, &lb); SCICOINOR_ERROR;

  _SciErr = getVarAddressFromPosition(pvApiCtx, BTYPE_IN, &btype_addr); SCICOINOR_ERROR;
  getAllocatedSingleString(pvApiCtx, btype_addr, &btype);
  _SciErr = getVarAddressFromPosition(pvApiCtx, VARTYPE_IN, &vartype_addr); SCICOINOR_ERROR;
  getAllocatedSingleString(pvApiCtx, vartype_addr, &vartype);
  _SciErr = getVarAddressFromPosition(pvApiCtx, PB_NAME_IN, &pb_name_addr); SCICOINOR_ERROR;
  getAllocatedSingleString(pvApiCtx, pb_name_addr, &pb_name);

  _SciErr = getVarAddressFromPosition(pvApiCtx, COL_NAME_IN, &col_name_addr); SCICOINOR_ERROR;
  getAllocatedMatrixOfString(pvApiCtx, col_name_addr, &n_col_name, &m_col_name, &pColNames);
  _SciErr = getVarAddressFromPosition(pvApiCtx, ROW_NAME_IN, &row_name_addr); SCICOINOR_ERROR;
  getAllocatedMatrixOfString(pvApiCtx, row_name_addr, &n_row_name, &m_row_name, &pRowNames);

  _SciErr = getVarAddressFromPosition(pvApiCtx, FILENAME_IN, &filename_addr); SCICOINOR_ERROR;
  getAllocatedSingleString(pvApiCtx, filename_addr, &filename);

  _SciErr = getVarAddressFromPosition(pvApiCtx, VERBOSE_IN, &verbose_addr); SCICOINOR_ERROR;
  getScalarDouble(pvApiCtx, verbose_addr, &verbose);

  nrows = n_lhs * m_lhs;
  if (nrows==0) nrows = n_rhs*m_rhs;
  ncols = n_c * m_c;

#ifdef DEBUG  
  DBGPRINTF("rowname n = %d m = %d\n", n_row_name, m_row_name);
#endif

  if (n_row_name * m_row_name==0)
    {
      RowNames.resize(nrows);
      for(i=0;i<nrows;i++)
	{
	  tmpstring.str("");
	  tmpstring << "R" << i;
	  RowNames[i] = tmpstring.str();
	}
    }
  else
    {
      RowNames.resize(nrows);
      for(i=0;i<nrows;i++)
	{
	  RowNames[i] = string(pRowNames[i]);
	}
    }

#ifdef DEBUG
  DBGPRINTF("colname n = %d m = %d\n", n_col_name, m_col_name);
#endif

  if (n_col_name * m_col_name==0)
    {
      ColNames.resize(ncols);
      for(i=0;i<ncols;i++)
	{
	  tmpstring.str("");
	  tmpstring << "R" << i;
	  ColNames[i] = tmpstring.str();
	}
    }
  else
    {
      ColNames.resize(ncols);
      for(i=0;i<ncols;i++)
	{
	  ColNames[i] = string(pColNames[i]);
	}
    }

  //////////////////
  // The A matrix //
  //////////////////

  _SciErr = getVarAddressFromPosition(pvApiCtx, A_IN, &a_addr); SCICOINOR_ERROR;
  _SciErr = getVarType(pvApiCtx, a_addr, &type); SCICOINOR_ERROR;

  if(type!=sci_sparse)
    {
      _SciErr = getMatrixOfDouble(pvApiCtx, a_addr, &n_a, &m_a, &a); SCICOINOR_ERROR;
      A_matrix.setDimensions(nrows,ncols);

      if (a==NULL) 
	{
	  Scierror(999,"%s: invalid value of matrix a\n",fname);

	  freeAllocatedSingleString(btype);
	  freeAllocatedSingleString(vartype);
	  freeAllocatedSingleString(pb_name);
	  freeAllocatedSingleString(filename);
	  if (pRowNames) freeAllocatedMatrixOfString(m_row_name, n_row_name, pRowNames);
	  if (pColNames) freeAllocatedMatrixOfString(m_col_name, n_col_name, pColNames);

	  return 0;
	}
      
      for(i=0; i<m_a; i++)
	{
	  for(j=0; j<n_a; j++)
	    {
	      if (*(a+i+j*m_a) != 0) A_matrix.modifyCoefficient(i,j,*(a+i+j*m_a));
	    }
	}
    }
  else
    {
      getAllocatedSparseMatrix(pvApiCtx, a_addr, &S_A.m, &S_A.n, &S_A.nel, &S_A.mnel, &S_A.icol, &S_A.R);
      A_matrix.setDimensions(nrows,ncols);
      
      count = 0;
      for(i=0;i<S_A.m;i++)
	{
	  if (S_A.mnel[i]!=0) 
	    {
	      for(j=0;j<S_A.mnel[i];j++)
		{
		  count++;
		  A_matrix.modifyCoefficient(i,S_A.icol[count-1]-1,S_A.R[count-1]);
		}
	    }
	}

      freeAllocatedSparseMatrix(S_A.mnel, S_A.icol, S_A.R);
    }

  if (verbose)
    {
      printer = new DerivedHandler();
      printer->setLogLevel((int)verbose);
      mpsWriter.passInMessageHandler(printer);
    }


  // void setMpsData (const CoinPackedMatrix &m, const double infinity,
  //                  const double *collb, const double *colub, 
  //                  const double *obj, const char *integrality,
  //                  const double *rowlb, const double *rowub,
  //                  const std::vector< std::string > &colnames, const std::vector< std::string > &rownames)

  mpsWriter.setMpsData(A_matrix, mpsWriter.getInfinity(),
		       lb, ub, 
		       c, vartype,
		       lhs, rhs,
		       ColNames, RowNames);

  mpsWriter.setProblemName(pb_name);

  status = mpsWriter.writeMps(filename);

  createScalarDouble(pvApiCtx, STATUS_OUT, (double)status);

  LhsVar(1) = STATUS_OUT;

  freeAllocatedSingleString(btype);
  freeAllocatedSingleString(vartype);
  freeAllocatedSingleString(pb_name);
  freeAllocatedSingleString(filename);

  if (pRowNames) freeAllocatedMatrixOfString(m_row_name, n_row_name, pRowNames);
  if (pColNames) freeAllocatedMatrixOfString(m_col_name, n_col_name, pColNames);

  return 0;
}
