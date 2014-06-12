//////////////////////////////////////////////////////////////////////////////////
// scicsdp: A scilab interface to the CSDP library for semidefinite programming //
//////////////////////////////////////////////////////////////////////////////////

//  Copyright (C) 2009-2010 Yann Collette.
//
//  SCICSDP is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2, or (at your option) any
//  later version.
//
//  This part of code is distributed with the FURTHER condition that it 
//  can be compiled and linked with the Matlab libraries and it can be 
//  used within the Matlab environment.
//
//  SCICSDP is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License
//  If not, write to the Free Software Foundation, 
//  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#include <cstdio>
#include <cstring>
#include <exception>

#include <api_scilab.h>

extern "C"
{
  // Scilab
#include <stack-c.h>
#include <sciprint.h>
#include <Scierror.h>
#include <MALLOC.h>
#include <api_parameters.h>
  // CSDP
#include <declarations.h>
}

#define CSDP_ERROR if(sciErr.iErr)		\
    {                                           \
      printError(&sciErr, 0);			\
      return sciErr.iErr;			\
    }

//#define DEBUG 1

using namespace std;

// Input arguments
#define C_IN       1
#define A_IN       2
#define B_IN       3
#define PARAM_IN   4
#define LAST_PARAM PARAM_IN

// Output arguments
#define X_OUT      5
#define Y_OUT      6
#define Z_OUT      7
#define F_OUT      8
#define STATUS_OUT 9
#define EXTRA_OUT  10

extern "C" int scicsdp(char * fname)
{
  // The problem and solution data.
  struct blockmatrix C;
  double *b = NULL;
  struct constraintmatrix *constraints = NULL;

  // Storage for the initial and final solutions.
  struct blockmatrix X,Z;
  double *y = NULL;
  double pobj,dobj;

  // blockptr will be used to point to blocks in constraint matrices.
  struct sparseblock *blockptr = NULL;
  struct paramstruc csdp_params;

  // A return code for the call to easy_sdp().
  int ret;

  // Get pointers to input
  // The variable A:
  // a list of list of matrix (sparse or not)
  // length(a) = the number of constraints
  // length(a(1)) = the number of blocks in constraints 1
  // The variable C:
  // a list of matrix (sparse or not)
  // length(c) = the number of blocks in matrix C
  int nb_vars, nb_constr;
  int i, j, k, l, ii, Index = 0, Index_sparse = 0, Log = 0;

  int c_constr_matr_nb_rows, c_constr_matr_nb_cols, c_block_list_nb_items;
  int * c_block_list_address = NULL, * c_constr_matr_address = NULL;
  double * c_constr_matr = NULL;
  int a_constr_matr_nb_rows, a_constr_matr_nb_cols, a_constr_matr_nb_items;
  int a_constr_list_nb_items, a_block_list_nb_items, * a_constr_matr_nb_items_row = NULL, * a_constr_matr_col_pos = NULL;
  int * a_constr_list_address = NULL, * a_block_list_address = NULL, * a_constr_matr_address = NULL;
  int m_size, issparse = 0;
  double * a_constr_matr = NULL;
  int b_vect_nb_rows, b_vect_nb_cols;
  int * b_vect_address = NULL;
  double * b_vect = NULL;
  int     tmp_int, tmp_res;
  double  tmp_double;
  char *  tmp_char = NULL;
  char *  writesdpa_filename = NULL;
  int     printlevel = 0;
  int * x_out_address = NULL;
  int * z_out_address = NULL;
  int f_out_nb_rows      = 1, f_out_nb_cols      = 1;
  int status_out_nb_rows = 1, status_out_nb_cols = 1;
  double * x_out = NULL, * y_out = NULL, * z_out = NULL, * f_out = NULL, * status_out = NULL;
  int * param_in_addr = NULL;
  int nb_entries = 0;
  // SDP parameters
  int type;
  SciErr sciErr;

  if (Rhs<LAST_PARAM) 
    {
      Scierror(999,"%s: %d inputs required in call to %s. Bug in csdp.sci ?...\n", LAST_PARAM, fname, fname);
      return 0;
    }

  //////////////////
  // The A matrix //
  //////////////////

  sciErr = getVarAddressFromPosition(pvApiCtx,A_IN,&a_constr_list_address); CSDP_ERROR;
  sciErr = getVarType(pvApiCtx, a_constr_list_address, &type); CSDP_ERROR;

  if (type!=sci_list)
    {
      Scierror(999,"%s: %the A matrix must be a list\n", fname);
      return 0;
    }

  sciErr = getListItemNumber(pvApiCtx,a_constr_list_address, &a_constr_list_nb_items); CSDP_ERROR;

  nb_constr = a_constr_list_nb_items;

  // The next major step is to setup the two constraint matrices A1 and A2. Again, because C indexing starts with 0, we have to allocate space for
  // one more constraint.  constraints[0] is not used.

  constraints = (struct constraintmatrix *)MALLOC((a_constr_list_nb_items+1)*sizeof(struct constraintmatrix));
  if (constraints==NULL)
    {
      Scierror(999,"%s: Failed to allocate storage for constraints!\n",fname);
      return 0;
    }

  for(i=0;i<a_constr_list_nb_items;i++)
    {
      // For each constraints
      sciErr = getListItemAddress(pvApiCtx,a_constr_list_address, i+1, &a_block_list_address); CSDP_ERROR;
      sciErr = getListItemNumber(pvApiCtx,a_block_list_address, &a_block_list_nb_items); CSDP_ERROR;

      constraints[i+1].blocks = NULL;

      //for(j=0;j<a_block_list_nb_items;j++)
      for(j=a_block_list_nb_items-1;j>=0;j--)
	{
	  // For each blocks, Setup the constraint matrix.  
	  sciErr = getListItemAddress(pvApiCtx,a_block_list_address, j+1, &a_constr_matr_address); CSDP_ERROR;
	  sciErr = getVarType(pvApiCtx, a_constr_matr_address, &type); CSDP_ERROR;
	  if (type==sci_sparse)
	    {
	      sciErr = getSparseMatrix(pvApiCtx,a_constr_matr_address, &a_constr_matr_nb_rows, &a_constr_matr_nb_cols, 
				       &a_constr_matr_nb_items, &a_constr_matr_nb_items_row, &a_constr_matr_col_pos, &a_constr_matr); CSDP_ERROR;
	      issparse = 1;
	      nb_entries = a_constr_matr_nb_items;
	    }
	  else
	    {
	      sciErr = getMatrixOfDouble(pvApiCtx,a_constr_matr_address, &a_constr_matr_nb_rows, &a_constr_matr_nb_cols, &a_constr_matr); CSDP_ERROR;
	      issparse = 0;

	      nb_entries = 0;
	      for(k=0; k<a_constr_matr_nb_rows*a_constr_matr_nb_cols; k++) 
		{
		  if (a_constr_matr[k]!=0.0) nb_entries++;
		}
	    }	      
#ifdef DEBUG
	  printf("DEBUG: matrix A - constraint %d - reading block %d: m = %d n = %d\n", i+1, j+1, a_constr_matr_nb_rows, a_constr_matr_nb_cols);
	  printf("DEBUG: nb_entries = %d\n", nb_entries);
#endif
	  // If the block matrix is empty, then, we skipped it.
	  if (nb_entries==0) 
	    {
#ifdef DEBUG
	      printf("DEBUG: skipping block %d from constraint %d.\n", j+1, i+1);
#endif
	      continue;
	    }

	  if (a_constr_matr_nb_rows != a_constr_matr_nb_cols)
	    {
	      Scierror(999,"%s: matrix A - the bloc %d is not square\n",fname, j+1);
	      if (constraints) FREE(constraints);
	      return 0;
	    }
	  
	  // Allocate space for block
	  blockptr = (struct sparseblock *)MALLOC(sizeof(struct sparseblock));
	  if (blockptr==NULL)
	    {
	      Scierror(999,"%s: Allocation of constraint block failed - step 1!\n",fname);
	      return 0;
	    }
	  
	  // For each block of constraint k
	  // Initialize block
	  blockptr->blocknum      = j+1;
	  blockptr->blocksize     = a_constr_matr_nb_rows;
	  blockptr->constraintnum = i+1;
	  blockptr->next          = NULL;
	  blockptr->nextbyblock   = NULL;
	  blockptr->numentries    = nb_entries;

#ifdef DEBUG
	  printf("DEBUG: blockptr: blocknum = %d, blocksize = %d, constraintnum = %d\n", blockptr->blocknum, blockptr->blocksize, blockptr->constraintnum); 
#endif

	  m_size = nb_entries;
	  blockptr->issparse = issparse;
	  
	  blockptr->entries = (double *)MALLOC((m_size+1)*sizeof(double));
	  if (blockptr->entries==NULL)
	    {
	      Scierror(999,"%s: Allocation of constraint block failed - step 2!\n",fname);
	      if (blockptr)    FREE(blockptr);
	      if (constraints) FREE(constraints);
	      return 0;
	    }
	  for(k=0; k<m_size+1; k++) blockptr->entries[k] = 0.0;
	  
	  blockptr->iindices = (int *) MALLOC((m_size+1)*sizeof(int));
	  if (blockptr->iindices==NULL)
	    {
	      Scierror(999,"%s: Allocation of constraint block failed - step 3!\n",fname);
	      if (blockptr->entries) FREE(blockptr->entries);
	      if (blockptr)          FREE(blockptr);
	      if (constraints)       FREE(constraints);
	      return 0;
	    }
	  for(k=0; k<m_size+1; k++) blockptr->iindices[k] = 0;
	  
	  blockptr->jindices = (int *) MALLOC((m_size+1)*sizeof(int));
	  if (blockptr->jindices==NULL)
	    {
	      Scierror(999,"%s: Allocation of constraint block failed - step 4!\n",fname);
	      if (blockptr->iindices) FREE(blockptr->iindices);
	      if (blockptr->entries)  FREE(blockptr->entries);
	      if (blockptr)           FREE(blockptr);
	      if (constraints)        FREE(constraints);
	      return 0;
	    }
	  for(k=0; k<m_size+1; k++) blockptr->jindices[k] = 0;
	  
	  // Now fill the blocks of the constraints

	  Index = 1;

	  if (issparse)
	    {
	      Index_sparse = 0;
	      for(k=0;k<a_constr_matr_nb_rows;k++)
		{
		  for(l=0;l<a_constr_matr_nb_items_row[k];l++)
		    {
		      // If a value is null, we skipped it
		      if (a_constr_matr[Index_sparse]==0.0) 
			{
#ifdef DEBUG
			  printf("DEBUG: skipping value %d - %d\n", k+1, a_constr_matr_col_pos[Index_sparse]);
#endif
			  Index_sparse++;
			  continue;
			}

		      // Put the entries into the first block.
		      blockptr->iindices[Index] = k+1;
		      blockptr->jindices[Index] = a_constr_matr_col_pos[Index_sparse];
		      blockptr->entries[Index]  = a_constr_matr[Index_sparse];
#ifdef DEBUG
		      printf("DEBUG: Index = %d, A[%d][%d] = %f\n",Index,blockptr->iindices[Index],blockptr->jindices[Index],blockptr->entries[Index]);
#endif
		      Index_sparse++;
		      Index++;
		    }
		}
	    }
	  else
	    {
	      for(k=0;k<blockptr->blocksize;k++)
		{
		  for(l=0;l<blockptr->blocksize;l++)
		    {
		      // If a value is null, we skipped it
		      if (a_constr_matr[k + l*blockptr->blocksize]==0.0) 
			{
#ifdef DEBUG
			  printf("DEBUG: skipping value %d - %d\n", k+1, l+1);
#endif
			  continue;
			}

		      // Put the entries into the first block.
		      blockptr->iindices[Index] = k+1;
		      blockptr->jindices[Index] = l+1;
		      blockptr->entries[Index]  = a_constr_matr[k + l*blockptr->blocksize];
#ifdef DEBUG
		      printf("DEBUG: Index = %d, A[%d][%d] = %f\n",Index, blockptr->iindices[Index],blockptr->jindices[Index],blockptr->entries[Index]);
#endif
		      Index++;
		    }
		}
	    }

#ifdef DEBUG
	  printf("DEBUG: j = %d / %d\n", j, a_block_list_nb_items);
#endif
	  // Insert block  into the linked list of A1 blocks.  
	  blockptr->next = constraints[i+1].blocks; // VALGRIND: Invalid read of size 8
	  constraints[i+1].blocks = blockptr; // VALGRIND: Invalid read of size 8
	}
#ifdef DEBUG
      printf("DEBUG: i = %d / %d\n", i, a_constr_list_nb_items);
#endif
    }
  
  //////////////////
  // The C matrix //
  //////////////////

  sciErr = getVarAddressFromPosition(pvApiCtx,C_IN,&c_block_list_address); CSDP_ERROR;
  sciErr = getVarType(pvApiCtx, c_block_list_address, &type); CSDP_ERROR;
  if (type!=sci_list)
    {
      Scierror(999,"%s: %the C matrix must be a list\n", fname);
      return 0;
    }

  sciErr = getListItemNumber(pvApiCtx,c_block_list_address, &c_block_list_nb_items); CSDP_ERROR;

  // First, allocate storage for the C matrix.  We have three blocks, but because C starts arrays with index 0, we have to allocate space for
  // four blocks- we'll waste the 0th block.  Notice that we check to make sure that the malloc succeeded.

  C.nblocks = c_block_list_nb_items;
  C.blocks  = (struct blockrec *)MALLOC((C.nblocks+1)*sizeof(struct blockrec));
  
#ifdef DEBUG
  printf("DEBUG: C matrix - %d items in the list\n", c_block_list_nb_items);
#endif

  if (C.blocks == NULL)
    {
      Scierror(999,"Couldn't allocate storage for C!\n");
      if (blockptr->jindices) FREE(blockptr->jindices);
      if (blockptr->iindices) FREE(blockptr->iindices);
      if (blockptr->entries)  FREE(blockptr->entries);
      if (blockptr)           FREE(blockptr);
      if (constraints)        FREE(constraints);
      return 0;
    }

  nb_vars = 0;

  for(i=0;i<c_block_list_nb_items;i++)
    {
      sciErr = getListItemAddress(pvApiCtx,c_block_list_address, i+1, &c_constr_matr_address); CSDP_ERROR;

#ifdef DEBUG
      printf("DEBUG: C matrix: type item %d\n", i+1);
#endif

      sciErr = getMatrixOfDouble(pvApiCtx,c_constr_matr_address, &c_constr_matr_nb_rows, &c_constr_matr_nb_cols, &c_constr_matr); CSDP_ERROR;

      if ((c_constr_matr_nb_rows != c_constr_matr_nb_cols) &&
	  !((c_constr_matr_nb_rows == 1) || (c_constr_matr_nb_cols == 1)))
	{
	  Scierror(999,"%s: matrix C - the bloc %d is not square or is not a vector\n",fname, i+1);
	  if (C.blocks)           FREE(C.blocks);
	  if (blockptr->jindices) FREE(blockptr->jindices);
	  if (blockptr->iindices) FREE(blockptr->iindices);
	  if (blockptr->entries)  FREE(blockptr->entries);
	  if (blockptr)           FREE(blockptr);
	  if (constraints)        FREE(constraints);
	  return 0;
	}

      // Setup the block.
      m_size = c_constr_matr_nb_rows * c_constr_matr_nb_rows;

      if ((c_constr_matr_nb_rows == 1) || (c_constr_matr_nb_cols == 1))
	{
	  C.blocks[i+1].blocksize     = m_size;
	  C.blocks[i+1].blockcategory = DIAG;
	  C.blocks[i+1].data.vec      = (double *)MALLOC((m_size + 1)* sizeof(double));
	}
      else
	{
	  C.blocks[i+1].blocksize     = c_constr_matr_nb_rows;
	  C.blocks[i+1].blockcategory = MATRIX;
	  C.blocks[i+1].data.mat      = (double *)MALLOC((m_size + 1)* sizeof(double));
	}
      
      nb_vars += c_constr_matr_nb_rows;

#ifdef DEBUG
      printf("DEBUG: matrix C: reading bloc %d: m = %d n = %d - size = %d\n", i+1, c_constr_matr_nb_rows, c_constr_matr_nb_cols, C.blocks[i+1].blocksize);
#endif
      
      if (C.blocks[i+1].blockcategory==MATRIX)
	{
	  if (C.blocks[i+1].data.mat == NULL)
	    {
	      Scierror(999,"%s: Couldn't allocate storage for the %dth block.\n",fname,C.blocks[i+1].blocksize);
	      if (C.blocks)           FREE(C.blocks);
	      if (blockptr->jindices) FREE(blockptr->jindices);
	      if (blockptr->iindices) FREE(blockptr->iindices);
	      if (blockptr->entries)  FREE(blockptr->entries);
	      if (blockptr)           FREE(blockptr);
	      if (constraints)        FREE(constraints);
	      return 0;
	    }
	}
      else
	{
	  if (C.blocks[i+1].data.vec == NULL)
	    {
	      Scierror(999,"%s: Couldn't allocate storage for the %dth block.\n",fname,C.blocks[i+1].blocksize);
	      if (C.blocks)           FREE(C.blocks);
	      if (blockptr->jindices) FREE(blockptr->jindices);
	      if (blockptr->iindices) FREE(blockptr->iindices);
	      if (blockptr->entries)  FREE(blockptr->entries);
	      if (blockptr)           FREE(blockptr);
	      if (constraints)        FREE(constraints);
	      return 0;
	    }
	}

      if (c_constr_matr==NULL) 
	{
	  Scierror(999,"%s: invalid value of matrix a\n",fname);
	  for(ii=0;ii<i+1;ii++)
	    {
	      if (C.blocks[ii+1].data.mat) FREE(C.blocks[ii+1].data.mat);
	    }
	  if (C.blocks)           FREE(C.blocks);
	  if (blockptr->jindices) FREE(blockptr->jindices);
	  if (blockptr->iindices) FREE(blockptr->iindices);
	  if (blockptr->entries)  FREE(blockptr->entries);
	  if (blockptr)           FREE(blockptr);
	  if (constraints)        FREE(constraints);
	  return 0;
	}

      if (C.blocks[i+1].blockcategory==MATRIX)
	{
	  for(j=0;j<c_constr_matr_nb_rows;j++)
	    {
	      for(k=0;k<c_constr_matr_nb_cols;k++)
		{
		  // Put the entries into the block.
		  C.blocks[i+1].data.mat[ijtok(j+1,k+1,C.blocks[i+1].blocksize)] = c_constr_matr[j+k*c_constr_matr_nb_rows];
#ifdef DEBUG
		  printf("DEBUG: MATRIX - C[%d][%d] = %f\n", j+1, k+1, C.blocks[i+1].data.mat[ijtok(j+1,k+1,C.blocks[i+1].blocksize)]);
#endif
		}
	    }
	}
      else
	{
	  for(j=0;j<m_size;j++)
	    {
	      // Put the entries into the block.
	      C.blocks[i+1].data.vec[j+1] = c_constr_matr[j];
#ifdef DEBUG
	      printf("DEBUG: DIAG - C[%d] = %f\n", j+1, C.blocks[i+1].data.vec[j+1]);
#endif
	    }
	}
    }

  //////////////////
  // The B vector //
  //////////////////

#ifdef DEBUG
  printf("DEBUG: reading the b vector\n");
#endif

  // YC: check the size of the b vector

  sciErr = getVarAddressFromPosition(pvApiCtx,B_IN,&b_vect_address); CSDP_ERROR;
  sciErr = getMatrixOfDouble(pvApiCtx,b_vect_address, &b_vect_nb_rows, &b_vect_nb_cols, &b_vect); CSDP_ERROR;

  // Allocate storage for the right hand side, b.
  m_size = b_vect_nb_rows * b_vect_nb_cols;
  b = (double *)MALLOC((m_size + 1)*sizeof(double));
  if (b == NULL)
    {
      Scierror(999,"%s: Failed to allocate storage for b!\n",fname);
      for(i=0;i<c_block_list_nb_items;i++)
	{
	  if (C.blocks[i+1].data.mat) FREE(C.blocks[i+1].data.mat);
	}
      if (C.blocks)           FREE(C.blocks);
      if (blockptr->jindices) FREE(blockptr->jindices);
      if (blockptr->iindices) FREE(blockptr->iindices);
      if (blockptr->entries)  FREE(blockptr->entries);
      if (blockptr)           FREE(blockptr);
      if (constraints)        FREE(constraints);
      return 0;
    }

  // Fill in the entries in b.
  for(i=0;i<m_size;i++) b[i+1] = b_vect[i];

  /////////////////////////
  // Process the options //
  /////////////////////////

#ifdef DEBUG
  printf("DEBUG: process the options\n");
#endif

  // Initialization of the parameters
  initparams(&csdp_params, &printlevel);

  initPList(pvApiCtx, PARAM_IN, &param_in_addr);

  if (!checkPList(pvApiCtx, param_in_addr))
    {
      Scierror(999, "%s: argument n° %d is not a plist\n", fname, PARAM_IN);

      if (b) FREE(b);
      for(i=0;i<c_block_list_nb_items;i++)
	{
	  if (C.blocks[i+1].data.mat) FREE(C.blocks[i+1].data.mat);
	}
      if (C.blocks)           FREE(C.blocks);
      if (blockptr->jindices) FREE(blockptr->jindices);
      if (blockptr->iindices) FREE(blockptr->iindices);
      if (blockptr->entries)  FREE(blockptr->entries);
      if (blockptr)           FREE(blockptr);
      if (constraints)        FREE(constraints);
      
      return 0;
    }
  
  // solver option
  getIntInPList(pvApiCtx, param_in_addr, "printlevel", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) printlevel = tmp_int;

  getDoubleInPList(pvApiCtx, param_in_addr, "axtol", &tmp_double, &tmp_res, 1.0e-8, Log, CHECK_NONE);
  if (tmp_res!=-1) csdp_params.axtol = tmp_double;

  getDoubleInPList(pvApiCtx, param_in_addr, "atytol", &tmp_double, &tmp_res, 1.0e-8, Log, CHECK_NONE);
  if (tmp_res!=-1) csdp_params.atytol = tmp_double;

  getDoubleInPList(pvApiCtx, param_in_addr, "objtol", &tmp_double, &tmp_res, 1.0e-8, Log, CHECK_NONE);
  if (tmp_res!=-1) csdp_params.objtol = tmp_double;

  getDoubleInPList(pvApiCtx, param_in_addr, "pinftol", &tmp_double, &tmp_res, 1.0e8, Log, CHECK_NONE);
  if (tmp_res!=-1) csdp_params.pinftol = tmp_double;

  getDoubleInPList(pvApiCtx, param_in_addr, "dinftol", &tmp_double, &tmp_res, 1.0e8, Log, CHECK_NONE);
  if (tmp_res!=-1) csdp_params.dinftol = tmp_double;

  getIntInPList(pvApiCtx, param_in_addr, "maxiter", &tmp_int, &tmp_res, 100, Log, CHECK_NONE);
  if (tmp_res!=-1) csdp_params.maxiter = tmp_int;

  getDoubleInPList(pvApiCtx, param_in_addr, "minstepfrac", &tmp_double, &tmp_res, 0.90, Log, CHECK_NONE);
  if (tmp_res!=-1) csdp_params.minstepfrac = tmp_double;

  getDoubleInPList(pvApiCtx, param_in_addr, "maxstepfrac", &tmp_double, &tmp_res, 0.97, Log, CHECK_NONE);
  if (tmp_res!=-1) csdp_params.maxstepfrac = tmp_double;

  getDoubleInPList(pvApiCtx, param_in_addr, "minstepp", &tmp_double, &tmp_res, 1.0e-8, Log, CHECK_NONE);
  if (tmp_res!=-1) csdp_params.minstepp = tmp_double;

  getDoubleInPList(pvApiCtx, param_in_addr, "minstepd", &tmp_double, &tmp_res, 1.0e-8, Log, CHECK_NONE);
  if (tmp_res!=-1) csdp_params.minstepd = tmp_double;

  getIntInPList(pvApiCtx, param_in_addr, "usexzgap", &tmp_int, &tmp_res, 2, Log, CHECK_NONE);
  if (tmp_res!=-1) csdp_params.usexzgap = tmp_int;

  getIntInPList(pvApiCtx, param_in_addr, "tweakgap", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) csdp_params.tweakgap = tmp_int;

  getIntInPList(pvApiCtx, param_in_addr, "affine", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) csdp_params.affine = tmp_int;

  getDoubleInPList(pvApiCtx, param_in_addr, "perturbobj", &tmp_double, &tmp_res, 1.0, Log, CHECK_NONE);
  if (tmp_res!=-1) csdp_params.perturbobj = tmp_double;

  getIntInPList(pvApiCtx, param_in_addr, "fastmode", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) csdp_params.fastmode = tmp_int;

  // writesdpa option
  getStringInPList(pvApiCtx, param_in_addr, "writesdpa", &tmp_char, &tmp_res, "test.sdpa", Log, CHECK_NONE);
  if (tmp_res!=-1) 
    {
      // Write the problem out in SDPA sparse format.
      write_prob(tmp_char,nb_vars,nb_constr,C,b,constraints);
      FREE(tmp_char);
#ifdef DEBUG
      printf("DEBUG: writesdpa_filename = %s\n", writesdpa_filename);
#endif
    }

  ////////////////
  // Setup done //
  ////////////////

#ifdef DEBUG
  printf("DEBUG: setup done\n");
  printf("nb_vars = %d, nb_constr = %d\n", nb_vars, nb_constr);
#endif

  // Create an initial solution.  This allocates space for X, y, and Z, and sets initial values.
  initsoln(nb_vars,nb_constr,C,b,constraints,&X,&y,&Z);

  // Solve the problem.
#ifdef DEBUG
  printf("DEBUG: solving the problem nb_vars = %d nb_constr = %d\n",nb_vars, nb_constr);
  printf("DEBUG: parameters:\n");
  printf("DEBUG:   axtol       = %e\n", csdp_params.axtol);
  printf("DEBUG:   atytol      = %e\n", csdp_params.atytol);
  printf("DEBUG:   objtol      = %e\n", csdp_params.objtol);
  printf("DEBUG:   pinftol     = %e\n", csdp_params.pinftol);
  printf("DEBUG:   dinftol     = %e\n", csdp_params.dinftol);
  printf("DEBUG:   maxiter     = %d\n", csdp_params.maxiter);
  printf("DEBUG:   minstepfrac = %e\n", csdp_params.minstepfrac);
  printf("DEBUG:   maxstepfrac = %e\n", csdp_params.maxstepfrac);
  printf("DEBUG:   minstepp    = %e\n", csdp_params.minstepp);
  printf("DEBUG:   minstepd    = %e\n", csdp_params.minstepd);
  printf("DEBUG:   usexzgap    = %d\n", csdp_params.usexzgap);
  printf("DEBUG:   tweakgap    = %d\n", csdp_params.tweakgap);
  printf("DEBUG:   affine      = %d\n", csdp_params.affine);
  printf("DEBUG:   perturbobj  = %e\n", csdp_params.perturbobj);
  printf("DEBUG:   fastmode    = %d\n", csdp_params.fastmode);
  printf("DEBUG:   printlevel  = %d\n", printlevel);
#endif

  ret = sci_easy_sdp(nb_vars,nb_constr,C,b,constraints,0.0,&X,&y,&Z,&pobj,&dobj, csdp_params, printlevel);
  
#ifdef DEBUG
  printf("DEBUG: end of resolution - ret = %d\n", ret);
#endif
  
  ///////////////////////////////////
  // End of copy paste of easy_sdp //
  ///////////////////////////////////
  
  // Process X
  sciErr = createList(pvApiCtx,X_OUT,X.nblocks,&x_out_address); CSDP_ERROR;

#ifdef DEBUG
  printf("DEBUG: x-out: %d blocks to process\n", X.nblocks);
#endif

  for(i=0;i<X.nblocks;i++)
    {
#ifdef DEBUG
      printf("DEBUG: x-out - processing block %d - size = %d\n", i+1, X.blocks[i+1].blocksize);
      printf("DEBUG: x-out - category = %d\n", X.blocks[i+1].blockcategory);
#endif
      sciErr = allocMatrixOfDoubleInList(pvApiCtx,X_OUT, x_out_address, i+1, X.blocks[i+1].blocksize, X.blocks[i+1].blocksize, &x_out); CSDP_ERROR; // VALGRIND: conditional jump depends on uninitialied value
      m_size = X.blocks[i+1].blocksize * X.blocks[i+1].blocksize;
      for(j=0;j<X.blocks[i+1].blocksize;j++) // VALGRIND: Conditional jump depends on uninitialized values
	{
	  for(k=0;k<X.blocks[i+1].blocksize;k++) // VALGRIND: idem
	    {
	      x_out[j + k*X.blocks[i+1].blocksize] = X.blocks[i+1].data.mat[ijtok(j+1,k+1,X.blocks[i+1].blocksize)]; // VALGRIND: Use of uninitialized value of size 8
	    }
	}
    }

  // Process Y
  sciErr = allocMatrixOfDouble(pvApiCtx,Y_OUT, nb_constr, 1, &y_out); CSDP_ERROR; // VALGRIND: conditional jum depends on uninitialized value
  memcpy(y_out,y,nb_constr*sizeof(double)); // VALGRIND: idem

  // Process Z
  sciErr = createList(pvApiCtx,Z_OUT,Z.nblocks,&z_out_address); CSDP_ERROR; // VALGRIND: idem

#ifdef DEBUG
  printf("DEBUG: z-out: %d blocks to process\n", Z.nblocks);
#endif

  for(i=0;i<Z.nblocks;i++)
    {
      sciErr = allocMatrixOfDoubleInList(pvApiCtx,Z_OUT, z_out_address, i+1, Z.blocks[i+1].blocksize, Z.blocks[i+1].blocksize, &z_out); CSDP_ERROR; // VALGRIND: use of uninitialized value of size 8

#ifdef DEBUG
      printf("DEBUG: z-out - processing block %d - size = %d\n", i+1, Z.blocks[i+1].blocksize);
      printf("DEBUG: z-out - caterogy = %d\n", Z.blocks[i+1].blockcategory);
#endif

      m_size = Z.blocks[i+1].blocksize * Z.blocks[i+1].blocksize;
      for(j=0;j<Z.blocks[i+1].blocksize;j++)
	{
	  for(k=0;k<Z.blocks[i+1].blocksize;k++)
	    {
	      z_out[j + k*Z.blocks[i+1].blocksize] = Z.blocks[i+1].data.mat[ijtok(j+1,k+1,X.blocks[i+1].blocksize)];
	    }
	}
    }
  
#ifdef DEBUG
  printf("DEBUG: now processing F_OUT and STATUS_OUT\n");
#endif

  sciErr = allocMatrixOfDouble(pvApiCtx, F_OUT,      f_out_nb_rows,      f_out_nb_cols,      &f_out);  CSDP_ERROR;
  sciErr = allocMatrixOfDouble(pvApiCtx, STATUS_OUT, status_out_nb_rows, status_out_nb_cols, &status_out);  CSDP_ERROR;

  *status_out = ret;
  *f_out      = (dobj+pobj)/2;

  LhsVar(1) = X_OUT;
  LhsVar(2) = Y_OUT;
  LhsVar(3) = Z_OUT;
  LhsVar(4) = F_OUT;
  LhsVar(5) = STATUS_OUT;

  // Free storage allocated for the problem and return.
  free_prob(nb_vars,nb_constr,C,b,constraints,X,y,Z);

  return 0;
}
