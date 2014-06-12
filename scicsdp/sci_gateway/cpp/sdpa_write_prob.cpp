/////////////////////////////////////////////////////////////////////////////////////////
// sdpa_read_prob: A scilab interface to the CSDP library for semidefinite programming //
/////////////////////////////////////////////////////////////////////////////////////////

//  Copyright (C) 2009-2010 Yann Collette.
//
//  sdpa_read_prob is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2, or (at your option) any
//  later version.
//
//  This part of code is distributed with the FURTHER condition that it 
//  can be compiled and linked with the Matlab libraries and it can be 
//  used within the Matlab environment.
//
//  sdpa_read_prob is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License
//  If not, write to the Free Software Foundation, 
//  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#include <string.h>
#include <stdio.h>

#define NOSHORTS 1

extern "C" {
#include <stack-c.h>
#include <Scierror.h>
#include <sciprint.h>
#include <declarations.h>
#include <MALLOC.h>
}

#include <api_scilab.h>

//#define DEBUG 1

#define FILENAME_IN 1
#define C_IN        2
#define A_IN        3
#define B_IN        4
#define STATUS_OUT  5

#define CSDP_ERROR if(sciErr.iErr)		\
    {                                           \
      printError(&sciErr, 0);			\
      return sciErr.iErr;			\
    }

extern "C" int sci_sdpa_write_prob(char * fname)
{
  // The problem and solution data.
  int nb_vars, nb_constr;
  struct blockmatrix C;
  struct constraintmatrix * constraints = NULL;
  struct sparseblock * blockptr = NULL;
  int * a_constr_list_address = NULL, a_constr_list_nb_items, * a_block_list_address = NULL; 
  int a_block_list_nb_items, * a_constr_matr_address = NULL, a_constr_matr_nb_items;
  int * a_constr_matr_col_pos = NULL, * a_constr_matr_nb_items_row = NULL;
  int a_constr_matr_nb_rows, a_constr_matr_nb_cols;
  double * a_constr_matr = NULL;
  int issparse;
  int * c_block_list_address = NULL, c_block_list_nb_items, * c_constr_matr_address = NULL;
  int c_constr_matr_nb_rows, c_constr_matr_nb_cols, m_size, ret;
  double * c_constr_matr = NULL;
  int * b_vect_address = NULL, b_vect_nb_rows, b_vect_nb_cols;
  double * b_vect = NULL, * b = NULL;
  double * status_out = NULL;
  int minrhs = 4, maxrhs = 4;
  int minlhs = 0, maxlhs = 1;
  char ** filename = NULL;
  int * filename_address = NULL, filename_rows, filename_cols, * filename_length = NULL;
  int i, j, k, ii, l, Index, Index_sparse, type, nb_entries;
  SciErr sciErr;

  CheckRhs(minrhs, maxrhs);
  CheckLhs(minlhs, maxlhs);

  ///////////////////////
  // Read the filename //
  ///////////////////////

  sciErr = getVarAddressFromPosition(pvApiCtx,FILENAME_IN,&filename_address); CSDP_ERROR;
  sciErr = getMatrixOfString(pvApiCtx, filename_address, &filename_rows, &filename_cols, NULL, NULL); CSDP_ERROR;
  if (filename_rows*filename_cols!=1)
    {
      Scierror(999,"%s: an unique string is requested for the filename\n",fname);
      return 0;
    }
  filename_length = (int *)MALLOC(filename_rows*filename_cols*sizeof(int));
  sciErr = getMatrixOfString(pvApiCtx,filename_address, &filename_rows, &filename_cols, filename_length, NULL); CSDP_ERROR;
  filename = (char **)MALLOC(filename_rows*filename_cols*sizeof(char *));
  for(i=0;i<filename_rows*filename_cols;i++)
    {
      filename[i] = (char *)MALLOC((filename_length[i]+1)*sizeof(char));
    }
  sciErr = getMatrixOfString(pvApiCtx, filename_address, &filename_rows, &filename_cols, filename_length, filename); CSDP_ERROR;

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
#ifdef DEBUG
	  printf("DEBUG: block %d / %d\n", j, a_block_list_nb_items);
#endif
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
	  printf("DEBUG: matrix A - constraint %d - block %d: size m = %d n = %d\n", i+1, j+1, a_constr_matr_nb_rows, a_constr_matr_nb_cols);
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
	      Scierror(999,"%s: matrix A - the block %d is not square\n",fname, j+1);
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
	  printf("DEBUG: blockptr: blocknum = %d, blocksize = %d, constraintnum = %d nbentries = %d\n", blockptr->blocknum, blockptr->blocksize, blockptr->constraintnum,blockptr->numentries); 
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
		      printf("DEBUG: Constraint = %d, Block = %d -Index = %d - A[%d][%d] = %f\n", i+1, j+1, Index, blockptr->iindices[Index],blockptr->jindices[Index],blockptr->entries[Index]);
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
		      printf("DEBUG: Constraint = %d, Block = %d - Index = %d - A[%d][%d] = %f\n", i+1, j+1, Index, blockptr->iindices[Index],blockptr->jindices[Index],blockptr->entries[Index]);
#endif
		      Index++;
		    }
		}
	    }

#ifdef DEBUG
	  printf("DEBUG: block = %d / %d\n", j, a_block_list_nb_items);
#endif
	  // Insert block  into the linked list of A1 blocks.  
	  blockptr->next = constraints[i+1].blocks; // VALGRIND: Invalid read of size 8
	  constraints[i+1].blocks = blockptr; // VALGRIND: Invalid read of size 8
	}
#ifdef DEBUG
      printf("DEBUG: constraint = %d / %d\n", i, a_constr_list_nb_items);
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

      if ((c_constr_matr_nb_rows != c_constr_matr_nb_cols) || 
	  (c_constr_matr_nb_rows != 1) ||
	  (c_constr_matr_nb_cols != 1))
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
	  C.blocks[i+1].blockcategory = DIAG;
	  C.blocks[i+1].blocksize     = m_size;
	  C.blocks[i+1].data.vec      = (double *)MALLOC((m_size + 1)* sizeof(double));
	}
      else
	{
	  C.blocks[i+1].blockcategory = MATRIX;
	  C.blocks[i+1].blocksize     = c_constr_matr_nb_rows;
	  C.blocks[i+1].data.mat      = (double *)MALLOC((m_size + 1)* sizeof(double));
	}

      nb_vars += c_constr_matr_nb_rows;

#ifdef DEBUG
      printf("DEBUG: matrix C: reading bloc %d: m = %d n = %d - size = %d\n", i+1, c_constr_matr_nb_rows, c_constr_matr_nb_cols, C.blocks[i+1].blocksize);
#endif
      
      if (C.blocks[i+1].blockcategory == MATRIX)
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

      if (C.blocks[i+1].blockcategory == DIAG)
	{
	  for(j=0;j<m_size;j++)
	    {
	      C.blocks[i+1].data.vec[j+1] = c_constr_matr[j];
#ifdef DEBUG
	      printf("DEBUG: DIAG - C[%d] = %f\n", j+1 C.blocks[i+1].data.vec[j+1]);
#endif
	    }
	}
      else
	{
	  for(j=0;j<c_constr_matr_nb_rows;j++)
	    {
	      for(k=0;k<c_constr_matr_nb_cols;k++)
		{
		  // Put the entries into the block.
		  C.blocks[i+1].data.mat[ijtok(j+1,k+1,C.blocks[i+1].blocksize)] = c_constr_matr[j+k*c_constr_matr_nb_rows];
#ifdef DEBUG
		  printf("DEBUG: C[%d][%d] = %f\n", j+1, k+1, C.blocks[i+1].data.mat[ijtok(j+1,k+1,C.blocks[i+1].blocksize)]);
#endif
		}
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

  ///////////////////////
  // Write the problem //
  ///////////////////////

#ifdef DEBUG
  for (l=1; l<=C.nblocks; l++)
    {
      for (i=1; i<=C.blocks[l].blocksize; i++)
	for (j=i; j<=C.blocks[l].blocksize; j++)
	  {
	    if (C.blocks[l].data.mat[ijtok(i,j,C.blocks[l].blocksize)] != 0.0 )
	      printf("DEBUG: C - 0 %d %d %d %.18e \n",l,i,j,C.blocks[l].data.mat[ijtok(i,j,C.blocks[l].blocksize)]);
	  }
    }

  struct sparseblock * p = NULL;

  for(i=1; i<=nb_constr; i++)
    {
      printf("DEBUG: constraint %d\n",i);
      printf("block num     = %d\n", constraints[i].blocks->blocknum);
      printf("blocksize     = %d\n", constraints[i].blocks->blocksize);
      printf("constraintnum = %d\n", constraints[i].blocks->constraintnum);
      printf("next          = %x\n", constraints[i].blocks->next);
      printf("numentries    = %d\n", constraints[i].blocks->numentries);

      p = constraints[i].blocks;
      while (p!=NULL)
	{
	  for (j=1; j<=p->numentries; j++)
	    {
	      printf("DEBUG: %d %d %d %d %.18e \n",i,p->blocknum,
		     p->iindices[j],
		     p->jindices[j],
		     p->entries[j]);
	    };
	  p = p->next;
	};
    }

  printf("DEBUG: nb_vars = %d, nb_constr = %d\n", nb_vars, nb_constr);
#endif

  ret = write_prob(filename[0], nb_vars, nb_constr, C, b, constraints);

  //////////////////////////////
  // Create the status output //
  //////////////////////////////

  sciErr = allocMatrixOfDouble(pvApiCtx, STATUS_OUT, 1, 1, &status_out); CSDP_ERROR;

  *status_out = (double)ret;

  LhsVar(1) = STATUS_OUT;

  return 0;
}
