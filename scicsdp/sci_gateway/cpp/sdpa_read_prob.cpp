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

#define NOSHORTS 1

extern "C" {
#include <stack-c.h>
#include <Scierror.h>
#include <sciprint.h>
#include <declarations.h>
#include <MALLOC.h>
#include <freeArrayOfString.h>
}

#include <api_scilab.h>

//#define DEBUG 1

#define FILENAME_IN   1
#define PRINTLEVEL_IN 2
#define C_OUT         3
#define A_OUT         4
#define B_OUT         5
#define STATUS_OUT    6

#define CSDP_ERROR if(sciErr.iErr)		\
    {                                           \
      printError(&sciErr, 0);			\
      return sciErr.iErr;			\
    }

extern "C" int sci_sdpa_read_prob(char * fname)
{
  int nb_vars, nb_constr, i, j, k, Index;
  struct blockmatrix pC;
  double * pa = NULL;
  struct constraintmatrix * pconstraints = NULL;
  struct sparseblock * tmp_block = NULL;
  int printlevel = 0, ret = 0, res = 0;
  char ** filename = NULL;
  int * filename_address = NULL, filename_rows, filename_cols, * filename_length = NULL;
  int * printlevel_address = NULL, printlevel_rows, printlevel_cols;
  double * printlevel_value = NULL, * c_out = NULL, * b_out = NULL;
  int * c_out_address = NULL, * a_out_address = NULL;
  int * constraint_address = NULL;
  int constraint_block_rows, constraint_block_cols, constraint_nb_item;
  double * block_value = NULL;
  int minrhs = 1, maxrhs = 2;
  int minlhs = 4, maxlhs = 4;
  SciErr sciErr;

  CheckRhs(minrhs,maxrhs);
  CheckLhs(minlhs,maxlhs);

  ///////////////////////
  // Read the filename //
  ///////////////////////

  sciErr = getVarAddressFromPosition(pvApiCtx,FILENAME_IN,&filename_address); CSDP_ERROR;
  sciErr = getMatrixOfString(pvApiCtx,filename_address, &filename_rows, &filename_cols, NULL, NULL); CSDP_ERROR;
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
  sciErr = getMatrixOfString(pvApiCtx,filename_address, &filename_rows, &filename_cols, filename_length, filename); CSDP_ERROR;

  /////////////////////
  // Read printlevel //
  /////////////////////

  if (Rhs==2)
    {
      sciErr = getVarAddressFromPosition(pvApiCtx,PRINTLEVEL_IN, &printlevel_address); CSDP_ERROR;
      sciErr = getMatrixOfDouble(pvApiCtx,printlevel_address, &printlevel_rows, &printlevel_cols, &printlevel_value); CSDP_ERROR;
      printlevel = (int)*printlevel_value;
    }
  else
    {
      printlevel = 0;
    }

  ///////////////////////////
  // Read the sdpa problem //
  ///////////////////////////

  ret =  read_prob(filename[0], &nb_vars, &nb_constr, &pC, &pa, &pconstraints, printlevel);

#ifdef DEBUG
  printf("DEBUG: read file %s, return code = %d nb_vars = %d nb_constr = %d\n", filename[0], ret, nb_vars, nb_constr);
#endif

  ///////////////
  // Process C //
  ///////////////

  sciErr = createList(pvApiCtx,C_OUT,pC.nblocks,&c_out_address); CSDP_ERROR;

#ifdef DEBUG
  printf("DEBUG: process c: %d blocks\n", pC.nblocks);
#endif

  for(i=0;i<pC.nblocks;i++)
    {
#ifdef DEBUG
      printf("DEBUG: processing block %d: blockssize = %d blockcat = %d\n", i+1, pC.blocks[i+1].blocksize, pC.blocks[i+1].blockcategory);
#endif
      if (pC.blocks[i+1].blockcategory==MATRIX)
	{
	  sciErr = allocMatrixOfDoubleInList(pvApiCtx, C_OUT, c_out_address, i+1, pC.blocks[i+1].blocksize, pC.blocks[i+1].blocksize, &c_out); CSDP_ERROR;
	  
	  for(j=0;j<pC.blocks[i+1].blocksize;j++)
	    {
	      for(k=0;k<pC.blocks[i+1].blocksize;k++)
		{
		  c_out[j + k*pC.blocks[i+1].blocksize] = pC.blocks[i+1].data.mat[ijtok(j+1,k+1,pC.blocks[i+1].blocksize)];
#ifdef DEBUG
		  printf("DEBUG: MATRIX - c_out_m[%d][%d] = %f - Index = %d\n",j,k,c_out[j + k*pC.blocks[i+1].blocksize], ijtok(j+1,k+1,pC.blocks[i+1].blocksize));
#endif
		}
	    }
	}
      else if (pC.blocks[i+1].blockcategory==DIAG)
	{
	  sciErr = allocMatrixOfDoubleInList(pvApiCtx, C_OUT, c_out_address, i+1, pC.blocks[i+1].blocksize, 1, &c_out); CSDP_ERROR;
	  
	  for(j=0;j<pC.blocks[i+1].blocksize;j++)
	    {
	      c_out[j] = pC.blocks[i+1].data.vec[j+1];
#ifdef DEBUG
	      printf("DEBUG: DIAG - c_out_d[%d] = %f\n",j,c_out[j]);
#endif
	    }	  
	}
      else
	{
	  Scierror(999,"%s: wrong blockcat type PACKEDMATRIX\n",fname);
	  return 0;
	}
    }

  ///////////////
  // Process A //
  ///////////////

#ifdef DEBUG
  printf("DEBUG: process a - nb_constr = %d\n",nb_constr);
#endif
  sciErr = createList(pvApiCtx,A_OUT,nb_constr,&a_out_address); CSDP_ERROR;

  for(i=0;i<nb_constr;i++)
    {
#ifdef DEBUG
      printf("DEBUG: processing constraint %d\n", i+1);
#endif
      // constraint_nb_item: nb blocks for constraint i+1
      Index = 0;
      tmp_block = pconstraints[i+1].blocks;
      while (tmp_block)
	{
	  tmp_block = tmp_block->next;
	  Index++;
	}
      constraint_nb_item = Index;

      sciErr = createListInList(pvApiCtx,A_OUT, a_out_address, i+1, constraint_nb_item, &constraint_address); CSDP_ERROR;

      Index = 0;
      tmp_block = pconstraints[i+1].blocks;
      for(j=0;j<constraint_nb_item;j++)
	{
#ifdef DEBUG
	  printf("DEBUG: processing block %d: nb_block = %d blocksize = %d numentries %d\n", j+1, constraint_nb_item, tmp_block->blocksize, tmp_block->numentries);
	  printf("DEBUG:                      blocknum = %d constraintnum = %d issparse = %d\n", tmp_block->blocknum, tmp_block->constraintnum, tmp_block->issparse);
#endif
	  constraint_block_rows = tmp_block->blocksize;
	  constraint_block_cols = tmp_block->blocksize;
#ifdef DEBUG
	  int ii;
	  for(ii=0;ii<tmp_block->numentries;ii++)
	    {
	      printf("entries[%d] = %f\n",ii+1,tmp_block->entries[ii+1]);
	    }
#endif
	  sciErr = allocMatrixOfDoubleInList(pvApiCtx, A_OUT, constraint_address, j+1, constraint_block_rows, constraint_block_cols, &block_value); CSDP_ERROR;
	  if (res)
	    {
	      Scierror(999,"%s: error while allocating a matrix A on the stack\n", fname);
	      return 0;
	    }
	  for(k=0;k<constraint_block_rows*constraint_block_cols;k++) block_value[k] = 0.0;

	  for(k=0;k<tmp_block->numentries;k++)
	    {
	      block_value[(tmp_block->iindices[k+1]-1) + (tmp_block->jindices[k+1]-1)*tmp_block->blocksize] = tmp_block->entries[k+1];
#ifdef DEBUG
	      printf("a_out[%d][%d] = %f\n",tmp_block->iindices[k+1],tmp_block->jindices[k+1],
			block_value[(tmp_block->iindices[k+1]-1)+(tmp_block->jindices[k+1]-1)*tmp_block->blocksize]);
#endif
	    }
	  tmp_block = tmp_block->next;
	}
    }

  ///////////////
  // Process B //
  ///////////////

#ifdef DEBUG
  printf("DEBUG: process B\n");
#endif
  sciErr = allocMatrixOfDouble(pvApiCtx, B_OUT, nb_constr, 1, &b_out); CSDP_ERROR;
  if (res)
    {
      Scierror(999,"%s: error while allocating a vector on the stack\n", fname);
      return 0;
    }
  for(i=0;i<nb_constr;i++) 
    {
      b_out[i] = pa[i+1];
#ifdef DEBUG
      printf("a[%d] = %f\n", i, b_out[i]);
#endif
    }

  ////////////////////
  // Process status //
  ////////////////////

#ifdef DEBUG
  printf("DEBUG: process status\n");
#endif
  sciErr = allocMatrixOfDouble(pvApiCtx, STATUS_OUT, 1, 1, &b_out); CSDP_ERROR;
  *b_out = (double)ret;
#ifdef DEBUG
  printf("DEBUG: status = %d,%d\n",*b_out,ret);
#endif

  LhsVar(1) = C_OUT;
  LhsVar(2) = A_OUT;
  LhsVar(3) = B_OUT;
  LhsVar(4) = STATUS_OUT;

  return 0;
}
